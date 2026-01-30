# This script does the following:
# 1. import workflow_adapters.py
# 2. input arguments are stored in global_vars dict based on options as follows:

# the following arguments are optional if default is given in [], otherwise mandatory

# --job_yaml:           job yaml file for the workflow_to_batch function
# --output_batch_file:   output batch file for the workflow_to_batch function
# --error_message: ['ERROR running worflow_to_batch on {}']
# --input_mount_path [/input]: maps to 'g_input_mount_path' in global_vars
# --local_workdir_path [/workdir]: maps to 'g_local_workdir_path'
# --pymipl_dir [/opt/packages/pymipl]: maps to 'g_pymipl_dir'
# --project: maps to 'g_project', XNAT project
# --main_repo_dir: maps to 'g_main_repo_dir', points to the repo dir of the script(s) executed at the run step of the job.

# 3. after global_vars dict is initialized, calls 'workflow_to_batch' which will write to output batch file.
# 4. outputs completion message

import argparse
import sys
from pathlib import Path

from workflow_adapters import workflow_to_batch  # :contentReference[oaicite:0]{index=0}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--job_yaml", required=True, help="Job YAML file for workflow_to_batch")
    p.add_argument("--output_batch_file", required=True, help="Output batch file (appended)")
    p.add_argument("--error_message", default="ERROR running worflow_to_batch on {}", help="Error message template")
    p.add_argument("--input_mount_path", default="/input")
    p.add_argument("--local_workdir_path", default="/workdir")
    p.add_argument("--pymipl_dir", default="/opt/packages/pymipl")
    p.add_argument("--main_env_dir", default="/opt/packages/user/env_repo")
    p.add_argument("--project", required=True)
    p.add_argument("--main_repo_dir", required=False)
    args = p.parse_args()

    global_vars = {
        "g_input_mount_path": args.input_mount_path,
        "g_local_workdir_path": args.local_workdir_path,
        "g_pymipl_dir": args.pymipl_dir,
        "g_project": args.project,
        "g_main_repo_dir": args.main_repo_dir,
    }

    job_yaml_path = Path(args.job_yaml)
    job_yaml_text = job_yaml_path.read_text()

    err = args.error_message.format(str(job_yaml_path))
    workflow_to_batch(job_yaml_text, global_vars, args.output_batch_file, error_message=err,step_command_prefix=f'micromamba run -p {args.main_env_dir}')

    print(f"Appended batch script: {args.output_batch_file}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"run_container_adapter.py failed: {e}", file=sys.stderr)
        raise
