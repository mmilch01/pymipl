from pathlib import Path
import yaml
from pyxnat import Interface
import requests
import re
import shutil
import logging
import tempfile
import zipfile
import os 

def sync_resource_xnat(local_resource, resource_name, project, subject=None, experiment=None, scan=None, upload=True, level="scan", XNAT_HOST=None, username=None, password=None, create_hierarchy=False, debug=False):

    if XNAT_HOST is None: XNAT_HOST = os.environ.get("XNAT_HOST")
    if username is None or password is None: username = os.environ.get("XNAT_USER"); password = os.environ.get("XNAT_PASS")
    if not XNAT_HOST or not username or not password: logging.error("Missing XNAT credentials (XNAT_HOST, user, password)"); return 2

    PROJECT_ID = str(project)
    params = { "inbody": "true", "PROJECT_ID": PROJECT_ID, "overwrite": "true" }

    subj = "" if subject is None else str(subject)
    exp = "" if experiment is None else str(experiment)
    scan_id = "" if scan is None else str(scan)

    levels = {'project': 1, 'subject': 2, 'experiment': 3, 'scan': 4}

    try:
        level_ord = levels[level]
        xnat = Interface(server=str(XNAT_HOST), user=username, password=password)
        project = xnat.select.project(PROJECT_ID)

        if level_ord >= levels['subject']:
            subj_obj = project.subject(subj); params['SUBJECT_ID'] = subj
            if not subj_obj.exists():
                if not upload or not create_hierarchy: 
                    logging.error(f"Subject {subj} does not exist")
                    return 1
                subj_obj.create()
                logging.info(f"Created subject '{subj}'")

        if level_ord >= levels['experiment']:
            exp_obj = subj_obj.experiment(exp); params['EXPT_LABEL'] = exp
            if not exp_obj.exists():
                if not upload or not create_hierarchy: 
                    logging.error(f"Experiment '{exp}' does not exist")
                    return 2
                exp_obj.create(); logging.info(f"Created experiment '{exp}'")

        if level_ord >= levels['scan']:
            scan_obj = exp_obj.scan(scan_id)
            if not scan_obj.exists():
                if not upload or not create_hierarchy: 
                    logging.error(f"Scan {scan_id} does not exist under experiment {exp}")
                    return 3
                scan_obj.create(); logging.info(f"Created scan '{scan_id}'")

        if level == "project": res_obj = project.resource(resource_name)
        elif level == "subject": res_obj = subj_obj.resource(resource_name)
        elif level == "experiment": res_obj = exp_obj.resource(resource_name)
        elif level == "scan": res_obj = scan_obj.resource(resource_name)

        if upload and not res_obj.exists(): res_obj.create()
        if not upload and not res_obj.exists(): 
            logging.error(f"Resource '{resource_name}' does not exist at level {level}")
            return 4

    except Exception as e:
        print(e)
        logging.error(e)
        return 5

    if level == "project": base_url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/resources/{resource_name}/files"
    if level == "subject": base_url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/subjects/{subj}/resources/{resource_name}/files"
    if level == "experiment": base_url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/subjects/{subj}/experiments/{exp}/resources/{resource_name}/files"
    if level == "scan": base_url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/subjects/{subj}/experiments/{exp}/scans/{scan_id}/resources/{resource_name}/files"

    if upload:
        tmp_dir, tmp_zip = None, None
        try:
            src = Path(local_resource)
            if not src.exists(): 
                logging.error(f"Source path does not exist: {local_resource}")
                return 6

            upload_items = []
            if src.is_file(): upload_items = [src]
            else:
                tmp_dir = Path(tempfile.mkdtemp(prefix="xnat_upload_"))
                base = tmp_dir / "payload"
                shutil.copytree(src, base)
                tmp_zip = Path(str(base) + ".zip")
                shutil.make_archive(str(base), "zip", root_dir=tmp_dir, base_dir="payload")
                upload_items = [tmp_zip]

            with requests.Session() as s:
                #main upload loop
                s.auth = (username, password)
                for item in upload_items:
                    with item.open("rb") as f:
                        is_zip_upload = (tmp_zip is not None and item == tmp_zip)
                        if is_zip_upload: 
                            params["extract"] = "true"
                            headers = {"Content-Type": "application/zip"}
                        else: 
                            headers = {"Content-Type": "application/octet-stream"}
                        r = s.put(f"{base_url}/{item.name}", params=params, data=f, headers=headers)
                        logging.info(f"Resource upload OK: {r.status_code} {base_url}")
                        if r.status_code >= 400: logging.error(f"Upload failed: {r.status_code} {base_url} {r.text[:1000]}"); return 2
                return 0                            
        except Exception as e:
            print(e)
            logging.error(e)
            return 7
        finally:
            if tmp_zip and tmp_zip.exists():
                try: tmp_zip.unlink()
                except Exception: pass
            if tmp_dir and tmp_dir.exists():
                try: shutil.rmtree(tmp_dir, ignore_errors=True)
                except Exception: pass
        return 0
        
    else: #download
        dst = Path(local_resource)
        wants_dir = str(local_resource).endswith("/") or str(local_resource).endswith("\\") or (dst.exists() and dst.is_dir())
    
        tmp_zip_path = Path(tempfile.mkstemp(prefix="xnat_resource_", suffix=".zip")[1])
        try:
            with requests.Session() as s:
                s.auth = (username, password)
                r = s.get(base_url, params={"format": "zip"}, stream=True)
                if r.status_code >= 400: 
                    logging.error(f"Download failed: {r.status_code} {base_url} {r.text[:1000]}")
                    return 8
                with tmp_zip_path.open("wb") as fo:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk: fo.write(chunk)
    
            with zipfile.ZipFile(tmp_zip_path, "r") as z:
                members = [m for m in z.namelist() if not m.endswith("/")]
                if not members: 
                    logging.error(f"Resource '{resource_name}' is empty")
                    return 9
    
                if wants_dir:
                    dst.mkdir(parents=True, exist_ok=True)
                    z.extractall(dst)
                    logging.info(f"Downloaded resource '{resource_name}' to directory {dst}")
                    return 0
    
                if len(members) != 1: logging.error(f"Resource '{resource_name}' has {len(members)} files; local_resource is a file path"); return 2
                dst.parent.mkdir(parents=True, exist_ok=True)
                with z.open(members[0], "r") as fi, dst.open("wb") as fo: shutil.copyfileobj(fi, fo)
                logging.info(f"Downloaded resource '{resource_name}' single file to {dst}")
                return 0
    
        except Exception as e:
            logging.error(e)
            return 10
        finally:
            if tmp_zip_path.exists():
                try: tmp_zip_path.unlink()
                except Exception: pass
        return 0


def init_global_vars_bootstrap_image(global_vars,xnat_project):
    global_vars['g_input_mount_path']=Path('/input')
    global_vars['g_local_workdir_path']=Path('/workdir')
    global_vars['g_pymipl_dir']=Path('/opt/packages/pymipl')
    global_vars['g_alg_repo_dir']=Path('/opt/packages/user/alg_repo')
    global_vars['g_project']=xnat_project

def workflow_to_batch(job_yaml, global_vars, output_batch_file, error_message=None, step_command_prefix=None):
    '''
    cmd must be single command with no |, >, &&, ;, variables, command substitutions, etc.
    Use step_command_prefix only for calls that need to execute in a custom environment (e.g. duneai);
    otherwise, environment with pymipl is assumed to be active at the time of batch execution.
    So, if your step runs a pymipl script, no step_command_prefix should be specified.
    For example, step_command_prefix can be 'micromamba -p /opt/packages/user/env_repo' for runs inside containers.
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
                #insert command prefix if not pymipl command.
                if step_command_prefix is not None and cmd.find('{g_pymipl_dir}') == -1:
                    f.write(f'cmd=({step_command_prefix} {cmd})\n')
                else: 
                    f.write(f'cmd=({cmd})\n')
                f.write(f'if ! "${{cmd[@]}}"; then\n    echo "{error_message}"\n    exit 1\nfi\n')
            args = [
                f'--xnat_project "{s_job["g_project"]}"',
                f'--subject "{s_job["job_subject"]}"',
                f'--experiment "{s_job["job_exp_label"]}"',
                f'--scan "{s_job["job_scan_id"]}"',
            ]
            #form upload command syntax
            upl = str(Path(s_job["g_pymipl_dir"]) / "sync_resource_with_xnat.py").format(**s_job)
            
            #write step commands to upload each specified local file to the resource.
            for res, files in (s.get("step_upload_files_to_resource") or {}).items():
                for p in files:
                    p = (str(p) if isinstance(p, Path) else str(p)).format(**s_job)
                    cmd = f'python "{upl}" {" ".join(args)} --source_loc "{p}" --resource_name "{res}"'
                    f.write(f'if ! [[ -f "{p}" ]]; then\n    echo "Missing file: {p}" \n    exit 1\nfi\n')
                    #note that there's no command prefix, because environment with pymipl will be expected to be active during execution.
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