#!/usr/bin/env bash

function exit_with_error
{
    echo "ERROR: $1"
    exit 1
}

if [ -z "$1" ]; then 
    echo "Usage: upload-env-to-xnat.sh <configuration file>"
    echo "Prepares custom resources to upload to XNAT project,"
    echo "including 1) user micromamba environment and 2) user code base"

    echo "  Configuration file should define the following variables: "
    echo "  UPLOAD_ENV (0 or 1): whether to upload environment\n"
    echo "  LOCAL_ENV_FOLDER: local directory with micromamba environment\n"
    echo "  REMOTE_ENV_RESOURCE: XNAT project resource containing environments, defaults to ENVS\n"
    #echo "  REMOTE_ENV_FOLDER: folder under XNAT project resource which will contain the uploaded environment\n"

    echo "  UPLOAD_USER_SOURCE (0 or 1): whether to upload user environment to XNAT project"
    echo "  REMOTE_USER_SOURCE_RESOURCE: XNAT project resource that contains user sources, defaults to SRC\n"
    echo "  LOCAL_USER_SOURCE_FOLDER: local folder containing user executables/scripts/resources\n"
    #echo "  REMOTE_USER_SOURCE_FOLDER: folder under XNAT project user sources resource which will contain the uploaded user source\n"


    echo "  SERVER      XNAT sever e.g. https://myxnat.org"
    echo "  PROJECT     Project where data will be stored, e.g MYPROJECT"
    echo "  Note that   server access credentials must be entered by user interactively when running this script."
    exit 1
fi

set -eu
config_file="$1"
[ -f "$config_file" ] || exit_with_error "configuration file $config_file does not exist"
source "$config_file"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
uploader="$script_dir/sync-resource-with-xnat.py"

[ -f "$uploader" ] || exit_with_error "upload helper script $uploader does not exist"
[[ "${SERVER:-}" != "" && "${PROJECT:-}" != "" ]] || exit_with_error "Config file must define SERVER and PROJECT"
(( UPLOAD_ENV == 1 || UPLOAD_USER_SOURCE == 1 )) || exit_with_error "Nothing to do: both UPLOAD_ENV and UPLOAD_USER_SOURCE are disabled"

#1. User environment.
if (( UPLOAD_ENV == 1 )); then
    [[ "${LOCAL_ENV_FOLDER:-}" != "" ]] || exit_with_error "LOCAL_ENV_FOLDER must be defined in config file with UPLOAD_ENV=1"

    remote_env_resource="${REMOTE_ENV_RESOURCE:-ENVS}"
    [ -d "$LOCAL_ENV_FOLDER" ] || exit_with_error "Folder in LOCAL_ENV_FOLDER $LOCAL_ENV_FOLDER does not exist"

    #credentials
    read -r -p "user:" XNAT_USER
    read -r -s -p "password:" XNAT_PASS

    python3 "$uploader" \
        --level project \
        --project "$PROJECT" \
        --local_resource "$LOCAL_ENV_FOLDER" \
        --remote_resource "$remote_env_resource" \
        --upload 1 \
        --xnat_host "$SERVER" \
        --user "$XNAT_USER" \
        --password "$XNAT_PASS" \
        --create_hierarchy 1    \
        --include_dir_under_resource 1
fi

#2 User (re)sources folder
if (( UPLOAD_USER_SOURCE == 1 )); then
    [[ -d "${LOCAL_USER_SOURCE_FOLDER:-}" ]] || exit_with_error "LOCAL_USER_SOURCE_FOLDER must point to existing folder when UPLOAD_USER_SOURCE=1"

    #credentials
    read -r -p "user:" XNAT_USER
    read -r -s -p "password:" XNAT_PASS

    remote_user_source_resource="${REMOTE_USER_SOURCE_RESOURCE:-SRC}"
    python3 "$uploader" \
        --level project \
        --project "$PROJECT" \
        --local_resource "$LOCAL_USER_SOURCE_FOLDER" \
        --remote_resource "$remote_user_source_resource" \
        --upload 1 \
        --xnat_host "$SERVER" \
        --user "$XNAT_USER" \
        --password "$XNAT_PASS" \
        --create_hierarchy 1  \
        --include_dir_under_resource 1        
fi

echo "Complete!"