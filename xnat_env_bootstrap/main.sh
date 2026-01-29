#!/bin/bash

#1 inputs session resource name with runtime configuration

#2 download, import and activate the specified micromamba environment from project resource

#3 read/interpret the runtime configuration from the session resource;

#4 either execute scripts directory from runtime confi(if they are stored verbatim), 
# or generate scripts and run them in this environment.


if [ -z "$1" ]; then 
    echo "usage: main.sh <project> <subject> <experiment> <runtime_resource> [options]"
    echo "options:"
    echo "  -user       <XNAT_USER>     take user from CL [\$XNAT_USER environment variable]"
    echo "  -pass       <XNAT_PASS>     take password from CL [\$XNAT_PASS environment variable]"
    echo "  -host       <XNAT_HOST>     take xnat host from CL [\$XNAT_HOST environment variable]"
    echo "  -job_yaml   <yaml file>     specify job Yaml file directly. <runtime_resource> will be used to store logs only."
    echo "  -job_script <bash script>   specify job batch script directly. <runtime_resource> will be used to store logs only."
    echo "  -microenv   <proj resource> XNAT project resource that contains micromamba environment for execution. Required."
    #option for nnUnet as an alternative to proj_resource, will be added later
    exit -1
fi

#1. Check inputs

#2. Define log file. 

#3. Download runtime environment from -microenv and import it into microconda. 
# is it possible to get its name automatically, and not ask for it? Or, assign a default name?

#4. If /input directory does not exist or is empty, download runtime_resource to temp dir if neither -job_yaml nor -job_script are specified.
# If /input directory exists, copy main.yaml and main.sh (if there) to workdir. 
# if at this point main.yaml is not found, exit with error.

# If 'main.sh' script is not in resource dir, a) generate it from main.yaml and b) set a flag variable that it was generated. Put main.sh in /workdir.
# How main.sh is generated: TODO


#5. Execute main.sh. That script would upload outputs back to xnat

#6. we don't need a separate log for now. stdout and stderr will be uploaded via container service
#7, However, we need to upload main.sh to the runtime resource if it was generated, so do that.
#8. Output successful completion message to stdout.
