#!/bin/bash

#1 inputs session resource name with runtime configuration

#2 download, import and activate the specified micromamba environment from project resource

#3 read/interpret the runtime configuration from the session resource;

#4 either execute scripts directory from runtime confi(if they are stored verbatim), 
# or generate scripts and run them in this environment.


set -eu

exit_with_error() { echo "ERROR: $*" 1>&2; usage 1>&2; exit 2; }

usage() {
  echo "usage: main.sh <project> <subject> <experiment> <workflow_id> [options]"
  echo ""
  echo "required positional:"
  echo "  project           XNAT project"
  echo "  subject           XNAT subject label"
  echo "  experiment        XNAT experiment label"
  echo "  workflow_id  XNAT session resource containing runtime config and, optionally, outputs"
  echo ""
  echo "Mounted environment options:"
  #echo "  -microenv <proj_resource>   XNAT project resource file with micromamba env tarball (e.g. ENVS/myenv.tar.gz)"
  echo "  -user-env <local_path>      User-supplied micromamba environment path"
  echo "  -user-src <local_path>      User-supplied source code/binaries path"
  echo ""
#  echo "repo source (exactly one required):"
#  echo "  -repo_git <url@sha>         public git repo URL with commit SHA (required); cloned to /opt/packages/user/alg_repo"
#  echo "  -repo_env_dir <path>        path inside extracted env that contains repo; symlinked to /opt/packages/user/alg_repo"
  echo ""
  echo "optional:"
  echo "  -host        <XNAT_HOST>        [default: \$XNAT_HOST]"
  echo "  -user        <XNAT_USER>        [default: \$XNAT_USER]"
  echo "  -pass        <XNAT_PASS>        [default: \$XNAT_PASS]"
  echo "  -input_mount <path>             [default: /input]"
  echo "  -pymipl_dir  <path>             [default: /opt/packages/pymipl]"
  #option for nnUnet as an alternative to proj_resource, will be added later
}


cmdline="$0 $*"
echo "running command: $cmdline"


if [[ $# -lt 4 ]]; then usage; exit -1; fi
PROJECT="$1"; SUBJECT="$2"; EXPERIMENT="$3"; WORKFLOW_ID="$4"; shift 4

# ---- defaults ----
XNAT_HOST="${XNAT_HOST}"
XNAT_USER="${XNAT_USER}"
XNAT_PASS="${XNAT_PASS}"
MICROENV_RESOURCE=""
REPO_GIT=""
USER_ENV_DIR=""
USER_SRC_DIR=""
INPUT_MOUNT="/input"
PYMIPL_DIR="/opt/packages/pymipl"

# ---- parse options ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    -host) XNAT_HOST="${2}"; shift 2 ;;
    -user) XNAT_USER="${2}"; shift 2 ;;
    -pass) XNAT_PASS="${2}"; shift 2 ;;
#    -microenv) MICROENV_RESOURCE="${2}"; shift 2 ;;
    -input_mount) INPUT_MOUNT="${2}"; shift 2 ;;
#    -repo_git) REPO_GIT="${2}"; shift 2 ;;
    -user-env) USER_ENV_DIR="${2}"; shift 2 ;;
    -user-src) USER_SRC_DIR="$2"; shift 2 ;;
    -pymipl_dir) PYMIPL_DIR="${2}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) exit_with_error "Unknown option: $1" ;;
  esac
done

if [ "$REPO_GIT" == "NONE" ]; then REPO_GIT=""; fi

# ---- validate required ----
#[[ -n "$MICROENV_RESOURCE" ]] || exit_with_error "Missing required: -microenv"
[[ -n "$XNAT_HOST" && -n "$XNAT_USER" && -n "$XNAT_PASS" ]] || exit_with_error "Missing XNAT credentials"

mkdir -p /workdir

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

run_python="micromamba run -n base python" 

# # Download runtime environment from -microenv and import it into microconda. 
# env_loc="$tmpdir/$MICROENV_RESOURCE"
# cmd=($run_python "$PYMIPL_DIR"/xnat_workflow/sync-resource-with-xnat.py \
#     --xnat_host "$XNAT_HOST"                \
#     --user "$XNAT_USER"                     \
#     --password "$XNAT_PASS"                 \
#     --level project                         \
#     --project "$PROJECT"                    \
#     --remote_resource "$MICROENV_RESOURCE"  \
#     --local_resource "$env_loc"             \
#     --upload 0)
#TODO DEBUG only
#echo "${cmd[@]}"
#"${cmd[@]}"
#if (( $? )); then exit_with_error "Failed to download microenv: $MICROENV_RESOURCE"; fi

mkdir -p /opt/packages/user

# extract runtime environment to workdir
#TODO DEBUG ONLY
env_repo_prefix="/opt/packages/user/env_repo"
if [ -d "$USER_ENV_DIR" ]; then 
  rm -rf "$env_repo_prefix"
  ln -sf "$USER_ENV_DIR" "$env_repo_prefix"
fi

#mkdir -p "$env_repo_prefix"
#rm -rf "$env_repo_prefix"/*
#echo tar -xzf "$env_loc" -C "$env_repo_prefix"
#tar -xzf "$env_loc" -C "$env_repo_prefix"
#if (( $? )); then exit_with_error "Failed to extract env tarball: $env_loc"; fi

# clone or link the main algorithm repo to standard location.
alg_repo_prefix="/opt/packages/user/alg_repo"
if [ -d "$USER_SRC_DIR" ]; then 
  rm -rf "$alg_repo_prefix"
  ln -sf "$USER_SRC_DIR" "$alg_repo_prefix"
fi

# repo_sha=""
# if [[ -n "$REPO_GIT" ]]; then
#   #clone main algorithm repo from public git link
#   case "$REPO_GIT" in
#     *@*) repo_url="${REPO_GIT%@*}"; repo_sha="${REPO_GIT#*@}";;
#     *) exit_with_error "-repo_git must be url@sha (SHA required)";;
#   esac  
#   [[ -n "$repo_url" && -n "$repo_sha" && "$repo_sha" != "$repo_url" ]] || exit_with_error "-repo_git must be url@sha (SHA required)"
#   echo git clone "$repo_url" "$alg_repo_prefix"
#   git clone "$repo_url" "$alg_repo_prefix"
#   (cd "$alg_repo_prefix" && git checkout "$repo_sha")
#   if (( $? )); then exit_with_error "git checkout failed: $repo_url:$repo_sha"; fi

# #link main algorithm repo from the downloaded environment
# elif [[ -n "$USER_ENV_DIR" ]]; then
#   repo_in_env="$USER_ENV_DIR"
#   case "$repo_in_env" in
#     /*) :;;
#     *) repo_in_env="$env_repo_prefix/$repo_in_env";;
#   esac
#   [[ -d "$repo_in_env" ]] || exit_with_error "-repo_env_dir does not exist after env extraction: $repo_in_env"
#   ln -s "$repo_in_env" "$alg_repo_prefix"
#   if (( $? )); then exit_with_error "Failed to link repo: $repo_in_env -> $alg_repo_prefix"; fi
# fi

resource_loc="$INPUT_MOUNT/RESOURCES/$WORKFLOW_ID"
job_sh_remote="$resource_loc/job.sh"
job_yaml_remote="$resource_loc/job.yaml"
job_sh="$tmpdir/job.sh"


#copy or compile job shell script
generated=0
if [ -f "$job_sh_remote" ]; then
    cp -f "$job_sh_remote" "$job_sh"
    if (( $? )); then exit_with_error "Failed to copy job script: $job_sh_remote"; fi

elif [ -f "$job_yaml_remote" ]; then 
    cmd=($run_python python "$PYMIPL_DIR"/xnat_workflow/run_container_adapter.py \
    --job_yaml "$job_yaml_remote" \
    --output_batch_file "$job_sh" \
    --input_mount_path "$INPUT_MOUNT" \
    --local_workdir_path "/workdir" \
    --pymipl_dir "$PYMIPL_DIR" \
    --project "$PROJECT" \
    --main_repo_dir "$alg_repo_prefix")
    echo "${cmd[@]}"
    "${cmd[@]}"
    if (( $? )); then exit_with_error "Failed to generate job script from: $job_yaml_remote"; fi
    generated=1
else
    #4. If $INPUT directory has no $job.yaml or $job.sh, exit with error
    exit_with_error "Neither job yaml $job_yaml_remote nor bash script $job_sh_remote was found"
fi

#5. Execute main.sh. That script would upload outputs back to xnat on its own.
chmod +x "$job_sh"
echo "$job_sh"
"$job_sh"
if (( $? )); then exit_with_error "Job script failed: $job_sh"; fi

# There's curretnly a catch: the script would need installed PYMIPL_DIR and pyxnat in the 
# default invironment but run the rest of the commands in downloaded environment. How to merge the two?
# Solution: the batch script must call each command using appropriate environment, not activate
# one environment and run everything in it. 

#6. write log to job runtime log directory. stdout and stderr will be uploaded via container service
ts="$(date +%m%d%H%M)"
job_log_dir="$tmpdir/job-$ts"
mkdir -p "$job_log_dir"
cp -f "$job_sh" "$job_log_dir/job.sh"

{
  echo "timestamp=$ts"
  echo "cmdline=$cmdline"
  echo "PROJECT=$PROJECT"
  echo "SUBJECT=$SUBJECT"
  echo "EXPERIMENT=$EXPERIMENT"
  echo "WORKFLOW_ID=$WORKFLOW_ID"
  echo "MICROENV_RESOURCE=$MICROENV_RESOURCE"
  echo "REPO_GIT=$REPO_GIT"
  echo "USER_ENV_DIR=$USER_ENV_DIR"
  echo "USER_SRC_DIR=$USER_SRC_DIR"
  #echo "repo_sha=$repo_sha"
  echo "INPUT_MOUNT=$INPUT_MOUNT"
  echo "PYMIPL_DIR=$PYMIPL_DIR"
  echo "XNAT_HOST=$XNAT_HOST"
} > "$job_log_dir/job.log"

#7, However, we need to upload $job_sh to the runtime resource if it was generated, so do that.
# Remember, resource mounts are not writable.
cmd=($run_python "$PYMIPL_DIR"/xnat_workflow/sync-resource-with-xnat.py \
    --xnat_host "$XNAT_HOST"                \
    --user "$XNAT_USER"                     \
    --password "$XNAT_PASS"                 \
    --level experiment                      \
    --project "$PROJECT"                    \
    --subject "$SUBJECT"                    \
    --experiment "$EXPERIMENT"              \
    --remote_resource "$WORKFLOW_ID"   \
    --local_resource "$job_log_dir"         \
    --upload 1)
echo "${cmd[@]}"
"${cmd[@]}"
if (( $? )); then exit_with_error "Failed uploading logs: $job_sh"; fi

echo "Completed main functionality"
