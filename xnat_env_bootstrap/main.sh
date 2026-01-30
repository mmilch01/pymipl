#!/bin/bash

#1 inputs session resource name with runtime configuration

#2 download, import and activate the specified micromamba environment from project resource

#3 read/interpret the runtime configuration from the session resource;

#4 either execute scripts directory from runtime confi(if they are stored verbatim), 
# or generate scripts and run them in this environment.


set -eu

exit_with_error() { echo "ERROR: $*" 1>&2; usage 1>&2; exit 2; }

usage() {
  echo "usage: main.sh <project> <subject> <experiment> <runtime_resource> [options]"
  echo ""
  echo "required positional:"
  echo "  project           XNAT project"
  echo "  subject           XNAT subject label"
  echo "  experiment        XNAT experiment label"
  echo "  runtime_resource  XNAT session resource containing runtime config (yaml/scripts)"
  echo ""
  echo "required options:"
  echo "  -microenv <proj_resource>   XNAT project resource file with micromamba env tarball (e.g. ENVS/myenv.tar.gz)"
  echo ""
  echo "repo source (exactly one required):"
  echo "  -repo_git <url@sha>         public git repo URL with commit SHA (required); cloned to /opt/packages/user/alg_repo"
  echo "  -repo_env_dir <path>        path inside extracted env that contains repo; symlinked to /opt/packages/user/alg_repo"
  echo ""
  echo "optional:"
  echo "  -host        <XNAT_HOST>        [default: \$XNAT_HOST]"
  echo "  -user        <XNAT_USER>        [default: \$XNAT_USER]"
  echo "  -pass        <XNAT_PASS>        [default: \$XNAT_PASS]"
  echo "  -job    <string>                job name. will look for \$job.yaml or \$job.sh."
  echo "  -input_mount <path>             [default: /input]"
  echo "  -pymipl_dir  <path>             [default: /opt/packages/pymipl]"
  #option for nnUnet as an alternative to proj_resource, will be added later
}

cmdline="$0 $*"

if [[ $# -lt 4 ]]; then usage; exit -1; fi
PROJECT="$1"; SUBJECT="$2"; EXPERIMENT="$3"; RUNTIME_RESOURCE="$4"; shift 4

# ---- defaults ----
XNAT_HOST="${XNAT_HOST:-}"
XNAT_USER="${XNAT_USER:-}"
XNAT_PASS="${XNAT_PASS:-}"
MICROENV_RESOURCE=""
JOB=""
REPO_GIT=""
REPO_ENV_DIR=""
INPUT_MOUNT="/input"
PYMIPL_DIR="/opt/packages/pymipl"

# ---- parse options ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    -host) XNAT_HOST="${2:-}"; shift 2 ;;
    -user) XNAT_USER="${2:-}"; shift 2 ;;
    -pass) XNAT_PASS="${2:-}"; shift 2 ;;
    -microenv) MICROENV_RESOURCE="${2:-}"; shift 2 ;;
    -job) JOB="${2:-}"; shift 2 ;;
    -input_mount) INPUT_MOUNT="${2:-}"; shift 2 ;;
    -repo_git) REPO_GIT="${2:-}"; shift 2 ;;
    -repo_env_dir) REPO_ENV_DIR="${2:-}"; shift 2 ;;
    -pymipl_dir) PYMIPL_DIR="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) exit_with_error "Unknown option: $1" ;;
  esac
done
# ---- validate required ----
[[ -n "$MICROENV_RESOURCE" ]] || exit_with_error "Missing required: -microenv"
[[ -n "$XNAT_HOST" && -n "$XNAT_USER" && -n "$XNAT_PASS" ]] || exit_with_error "Missing XNAT credentials"
[[ -n "$JOB" ]] || exit_with_error "Must specify job name"

mkdir -p /workdir

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# Download runtime environment from -microenv and import it into microconda. 
env_loc="$tmpdir/$MICROENV_RESOURCE"
python "$PYMIPL_DIR"/xnat_workflow/sync-resource-with-xnat.py \
    --xnat_host "$XNAT_HOST"                \
    --user "$XNAT_USER"                     \
    --password "$XNAT_PASS"                 \
    --level project                         \
    --project "$PROJECT"                    \
    --remote_resource "$MICROENV_RESOURCE"  \
    --local_resource "$env_loc"             \
    --upload 0

mkdir -p /opt/packages/user

# extract runtime environment to workdir
env_repo_prefix="/opt/packages/user/env_repo"
mkdir -p "$env_repo_prefix"
rm -rf "$env_repo_prefix"/*
tar -xzf "$env_loc" -C "$env_repo_prefix"

# clone or link the main algorithm repo to standard location.
alg_repo_prefix="/opt/packages/user/alg_repo"
rm -rf "$alg_repo_prefix"

repo_sha=""
if [[ -n "$REPO_GIT" ]]; then
  #clone main algorithm repo from public git link
  case "$REPO_GIT" in
    *@*) repo_url="${REPO_GIT%@*}"; repo_sha="${REPO_GIT#*@}";;
    *) exit_with_error "-repo_git must be url@sha (SHA required)";;
  esac  
  [[ -n "$repo_url" && -n "$repo_sha" && "$repo_sha" != "$repo_url" ]] || exit_with_error "-repo_git must be url@sha (SHA required)"
  git clone "$repo_url" "$alg_repo_prefix"
  (cd "$alg_repo_prefix" && git checkout "$repo_sha")
#link main algorithm repo from the downloaded environment
elif [[ -n "$REPO_ENV_DIR" ]]; then
  repo_in_env="$REPO_ENV_DIR"
  case "$repo_in_env" in 
    /*) :;; 
    *) repo_in_env="$env_repo_prefix/$repo_in_env";; 
  esac
  [[ -d "$repo_in_env" ]] || exit_with_error "-repo_env_dir does not exist after env extraction: $repo_in_env"
  ln -s "$repo_in_env" "$alg_repo_prefix"
fi


resource_loc="$INPUT_MOUNT/RESOURCES/$RUNTIME_RESOURCE"
job_sh_remote="$resource_loc/$JOB.sh"
job_yaml_remote="$resource_loc/$JOB.yaml"
job_sh="$tmpdir/$JOB.sh"


#copy or compile job shell script
generated=0
if [ -f "$job_sh_remote" ]; then 
    cp -f "$job_sh_remote" "$job_sh"

elif [ -f "$job_yaml_remote" ]; then 
    python "$PYMIPL_DIR"/xnat_workflow/run_container_adapter.py \
    --job_yaml "$job_yaml_remote" \
    --output_batch_file "$job_sh" \
    --input_mount_path "$INPUT_MOUNT" \
    --local_workdir_path "/workdir" \
    --pymipl_dir "$PYMIPL_DIR" \
    --project "$PROJECT" \
    --main_repo_dir "$alg_repo_prefix"
    generated=1
else 
    #4. If $INPUT directory has no $job.yaml or $job.sh, exit with error
    exit_with_error "Neither job yaml $job_yaml_remote nor bash script $job_sh_remote was found"
fi
#5. Execute main.sh. That script would upload outputs back to xnat on its own.
chmod +x "$job_sh"
"$job_sh" 

# There's curretnly a catch: the script would need installed PYMIPL_DIR and pyxnat in the 
# default invironment but run the rest of the commands in downloaded environment. How to merge the two?
# Solution: the batch script must call each command using appropriate environment, not activate
# one environment and run everything in it. 

#6. write log to job runtime log directory. stdout and stderr will be uploaded via container service
ts="$(date +%m%d%H%M)"
job_log_dir="$tmpdir/$JOB-$ts"
mkdir -p "$job_log_dir"
cp -f "$job_sh" "$job_log_dir/$JOB.sh"

{
  echo "timestamp=$ts"
  echo "cmdline=$cmdline"
  echo "PROJECT=$PROJECT"
  echo "SUBJECT=$SUBJECT"
  echo "EXPERIMENT=$EXPERIMENT"
  echo "RUNTIME_RESOURCE=$RUNTIME_RESOURCE"
  echo "MICROENV_RESOURCE=$MICROENV_RESOURCE"
  echo "REPO_GIT=$REPO_GIT"
  echo "REPO_ENV_DIR=$REPO_ENV_DIR"
  echo "repo_sha=$repo_sha"
  echo "INPUT_MOUNT=$INPUT_MOUNT"
  echo "PYMIPL_DIR=$PYMIPL_DIR"
  echo "XNAT_HOST=$XNAT_HOST"
  echo "JOB=$JOB"
} > "$job_log_dir/$JOB.log"

#7, However, we need to upload $job_sh to the runtime resource if it was generated, so do that.
# Remember, resource mounts are not writable.

python "$PYMIPL_DIR"/xnat_workflow/sync-resource-with-xnat.py \
    --xnat_host "$XNAT_HOST"                \
    --user "$XNAT_USER"                     \
    --password "$XNAT_PASS"                 \
    --level experiment                      \
    --project "$PROJECT"                    \
    --subject "$SUBJECT"                    \
    --experiment "$EXPERIMENT"              \
    --remote_resource "$RUNTIME_RESOURCE"   \
    --local_resource "$job_log_dir"         \
    --upload 1

echo "Completed main functionality"