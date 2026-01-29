from pathlib import Path
import yaml
from pyxnat import Interface
import requests
import re
import shutil
import logging
import tempfile

def paths_to_str(x):
    if isinstance(x, Path): return str(x)
    if isinstance(x, dict): return {k: paths_to_str(v) for k, v in x.items()}
    if isinstance(x, list): return [paths_to_str(v) for v in x]
    return x

def workflow_to_batch(job_yaml, global_vars, output_batch_file, error_message=None):
    '''
    cmd must be single command with no |, >, &&, ;, variables, command substitutions, etc.
    '''
        
    d = yaml.safe_load(job_yaml) if isinstance(job_yaml, str) else dict(job_yaml)
    job = dict(d)

    for k, v in global_vars.items(): job[k] = str(v) if isinstance(v, Path) else v

    out = Path(output_batch_file)
    out.parent.mkdir(parents=True, exist_ok=True)

    exists = out.exists()
    with out.open("a") as f:
        if not exists: f.write("#!/bin/bash\n")
        if 'job_title' in d:
            title=str(d['job_title']).format(**d)
            f.write(f'echo "############ JOB: {title}"\n')
        for s in d.get("steps", []):
            s_job = dict(job)
            for k, v in s.items():
                if isinstance(v, Path): s_job[k] = str(v)
                elif isinstance(v, (str, int, float, bool, dict)): s_job[k] = v
            if "step_title" in s:
                title = str(s["step_title"]).format(**s_job)
                f.write(f'echo "********** STEP: {title}"\n')
            #write step command with substituted values to batch file
            if "step_command" in s: 
                cmd = str(s["step_command"]).format(**s_job)
                if error_message is None: error_message=f'command failed: {cmd}'
                f.write(f'cmd=({cmd})\n')
                f.write(f'if ! "${{cmd[@]}}"; then\n    echo "{error_message}"\n    exit 1\nfi\n')
            args = [
                f'--xnat_project "{s_job["g_project"]}"',
                f'--subject "{s_job["job_subject"]}"',
                f'--experiment "{s_job["job_exp_label"]}"',
                f'--scan "{s_job["job_scan_id"]}"',
            ]
            #form upload command syntax
            upl = str(Path(s_job["g_pymipl_dir"]) / "upload_resource_to_xnat.py").format(**s_job)
            
            #write step commands to upload each specified local file to the resource.
            for res, files in (s.get("step_upload_files_to_resource") or {}).items():
                for p in files:
                    p = (str(p) if isinstance(p, Path) else str(p)).format(**s_job)
                    cmd = f'python "{upl}" {" ".join(args)} --source_loc "{p}" --resource_name "{res}"'
                    f.write(f'if ! [[ -f "{p}" ]]; then\n    echo "Missing file: {p}" \n    exit 1\nfi\n')
                    f.write(f'cmd=({cmd})\n')
                    f.write(f'if ! "${{cmd[@]}}"; then\n    printf \'%s\\n\' "failed command: $cmd"\n    exit 1\nfi\n')                    

            #write step commands to upload each specified local directory to the resource.
            for res, dirs in (s.get("step_upload_dirs_to_resource") or {}).items():
                for p in dirs:
                    p = (str(p) if isinstance(p, Path) else str(p)).format(**s_job)
                    cmd = f'python "{upl}" {" ".join(args)} --source_loc "{p}" --resource_name "{res}"'
                    f.write(f'if ! [[ -d "{p}" ]]; then\n   echo "Missing dir: {p}"\n    exit 1\nfi\n')
                    f.write(f'if ! [[ "$(ls -A "{p}" 2>/dev/null)" ]]; then\n    echo "Empty dir: {p}"\n    exit 1\nfi\n')
                    
                    f.write(f'cmd=({cmd})\n')
                    f.write(f'if ! "${{cmd[@]}}"; then\n    printf \'%s\\n\' "failed command: ${{cmd[*]}}"\n    exit 1\nfi\n')
                               

        #f.write(f'echo "############ END JOB: {title}"\n')
    out.chmod(0o755)