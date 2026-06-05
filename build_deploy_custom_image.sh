
#!/bin/bash

exit_with_error() {
    echo "${script_name}: ERROR: $*" >&2
    exit 1
}

script_name=$(basename "$0")
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

show_usage() {
    cat <<EOF
Usage: ${script_name} <config_file>

Build and deploy custom XNAT workflow bootstrap image. 

Example config: xnat_env_bootstrap/example-custom-image.conf
Config file must define:
  LOCAL_IMAGE       local Docker image name, e.g. xnat-env:latest
  REMOTE_IMAGE      remote Docker image name with tag, e.g. mmilch01/xnat-env:latest
  DOCKER_REPO       Docker repository, default: docker.io
  LOCAL_ENV_FOLDER  folder with user environment to embed into the Docker image.
  ENV_VARS          list extra environment variables for ENV dockerfile directive, e.g. "VAR1=val1 VAR2=val2".
EOF
}

config_file="$1"
if [ "$#" -ne 1 ]; then show_usage;  exit 1; fi
if [ ! -f "$config_file" ]; then show_usage; exit 1; fi
if [ ! -d "$PYMIPL_DIR" ]; then exit_with_error "PYMIPL_DIR must be set in order to proceed with the build."; fi

source "$config_file"
DOCKER_REPO="${DOCKER_REPO:-docker.io}"

[ -n "${LOCAL_IMAGE:-}" ] || exit_with_error "LOCAL_IMAGE is not set in config file."
[ -n "${REMOTE_IMAGE:-}" ] || exit_with_error "REMOTE_IMAGE is not set in config file."
[ -d "${LOCAL_ENV_FOLDER:-}" ] || exit_with_error "LOCAL_ENV_FOLDER $LOCAL_ENV_FOLDER must be a directory."

tmp_dockerfile="$(mktemp "${TMPDIR:-/tmp}/Dockerfile.custom.XXXXXX")"
#trap 'rm -f "$tmp_dockerfile"' EXIT
echo "using temp dockerfile: $tmp_dockerfile"

if [ -n "$ENV_VARS" ]; then 
  cp "$script_dir/xnat_env_bootstrap/Dockerfile.custom" "$tmp_dockerfile"
  printf '\nENV %s\n' "$ENV_VARS" >> "$tmp_dockerfile"
fi

echo "building image"

echo docker build --tag "$LOCAL_IMAGE" --build-context user_env="$LOCAL_ENV_FOLDER" -f "$tmp_dockerfile" $PYMIPL_DIR
docker build --tag "$LOCAL_IMAGE" --build-context user_env="$LOCAL_ENV_FOLDER" -f "$tmp_dockerfile" $PYMIPL_DIR || exit_with_error "Image build failed"

echo "deploying image"
docker tag "$LOCAL_IMAGE" "$REMOTE_IMAGE" || exit_with_error "docker tag failed: $LOCAL_IMAGE -> $REMOTE_IMAGE"
docker login "$DOCKER_REPO" || exit_with_error "docker login failed for repository: $DOCKER_REPO"
docker push "$REMOTE_IMAGE" || exit_with_error "docker push failed: $REMOTE_IMAGE"

echo "Finished deployment"
