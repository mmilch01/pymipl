# This script does the following: 
# 0. import workflow_adapters.py

# 2. input arguments are stored in global_vars dict based on options as follows:

# the following arguments are optional if default is given in [], otherwise mandatory

# --job_yaml:           job yaml file for the workflow_to_batch function
# --output_batch_file:   output batch file for the workflow_to_batch function
# --error_message: ['ERROR running worflow_to_batch on {}']
# --input_mount_path []: maps to 'g_input_mount_path' in global_vars
# --local_workdir_path []: maps to 'g_local_workdir_path'
# --pymipl_dir []: maps to 'g_pymipl_dir'
# --project: maps to 'g_project', XNAT project
# --main_repo_dir: maps to 'g_main_repo_dir', points to the repo dir of the script(s) executed at the run step of the job.

# 3. after global_vars dict is initialized, calls 'workflow_to_batch' which will write to output batch file.
# 4. outputs completion message