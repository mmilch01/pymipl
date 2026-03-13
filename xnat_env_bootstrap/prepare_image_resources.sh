#!/usr/bin/env bash

if [ -z "$1" ]; then 
    echo "Usage: prepare_image_resources.sh <configuration file>"
    echo ""

    echo "  Configuration file should define the following variables: "
    echo "  UPDATE_MODELS (0 or 1): whether to update models"
    echo "  ARCHIVED_MODELS_DIR: directory with archived nnunet models\n"

    echo "  LOAD_ENV    (0 or 1): whether to load user environment from server"
    echo "  SERVER      XNAT sever e.g. https://myxnat.org"
    echo "  PROJECT     Project where environment is stored, e.g MYPROJECT"
    echo "  RESOURCE    Project resource where environment is stored, e.g. ENVS"
    echo "  ENV_FILE    File name under the project resource, e.g. myenv.tar.gz"
    echo "  Note that   server access credentials must be entered by user interactively when running this script."
    exit 1
fi

set -eu
config_file="$1"
if [ ! -f "$config_file" ]; then 
    echo "configuration file $config_file does not exist"
    exit -1
fi
source "$config_file"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
#1. User environment.
# read the config file
if (( LOAD_ENV == 1 )); then
    source "$config_file"

    if [ -z "${SERVER:-}" ] || [ -z "${PROJECT:-}" ] || [ -z "${RESOURCE:-}" ] || [ -z "${ENV_FILE:-}" ]; then
        echo "Config file must define SERVER, PROJECT, RESOURCE, and ENV_FILE"
        exit 1
    fi
    # ask for user credentials for the specified SERVER
    read -r -p "user:" XNAT_USER
    read -r -s -p "password:" XNAT_PASS

    # using curl, via REST call download the environment from SERVER/PROJECT/RESOURCE/ENV_FILE, error out if not successful
    temp="$(mktemp)";  trap 'rm -f "$temp"' EXIT
    url="${SERVER}/data/projects/${PROJECT}/resources/${RESOURCE}/files/${ENV_FILE}"

    curl -fL -u "$XNAT_USER:$XNAT_PASS" "$url" -o "$temp" || {
        echo "Failed to download environment archive from $url"
        exit 1
    }

    # clear the contents of <SCRIPT_DIR>/user_env directory
    env_dir="$SCRIPT_DIR/user_env"
    rm -rf "$env_dir"; mkdir -p "$env_dir"

    # Unzip and untar if ENV_FILE extension is tar.gz or just untar if ENV_FILE extension is .tar into <SCRIPT_DIR>/user_env directory
    echo "extracting ${ENV_FILE}:"
    case "$ENV_FILE" in 
        *.tar.gz)  tar -xzf "$temp" -C "$env_dir" ;;
        *.tar)     tar -xf "$temp" -C "$env_dir" ;;
        *)         echo "Unsupported environment archive extension: $ENV_FILE";  exit 1 ;;
    esac    
fi

#2. nnunet models.
if (( UPDATE_MODELS )); then 
    TARGET_DIR="$SCRIPT_DIR/nnunet_models"

    mkdir -p "$TARGET_DIR"

    shopt -s nullglob
    for zip_file in "$ARCHIVED_MODELS_DIR"/*.zip; do
        unzip -u "$zip_file" -d "$TARGET_DIR"
    done
fi
