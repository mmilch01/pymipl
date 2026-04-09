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
import sys

# XNAT-Jupyter workflow helper functions
import pyxnat
import json
import requests
import datetime
import yaml

def populate_job_fields(env_type, global_vars, workflow_id, data_dict, subject=None, exp_label=None, scan_id=None,
                       subj_key='Subject',exp_key='Experiment',scan_key='Scan'):
    '''
    fills the 'job' dict with the values from the data_list csv. 
    special fields job_subject, job_exp_label, job_scan_id are 
    populated either from direcly specified subject, exp_label, scan_id parameters,
    or extracted using specified keys from data_dict. 
    '''

    job={'steps': []}
    for key,val in data_dict.items(): job[key]=val

    #define job_id, job_title, job_subject, job_exp_label, job_scan_id
    keys=data_dict.keys()
    job_id=workflow_id
    job_title=f'Workflow: {workflow_id}'
    job_workdir=global_vars['g_local_workdir_path']
    
    for var,csv_key,job_key in zip ((subject, exp_label, scan_id),
                                    (subj_key, exp_key, scan_key),
                                    ('job_subject','job_exp_label','job_scan_id')):
        if var is not None:         job[job_key]=var
        elif csv_key in keys:       job[job_key]=data_dict[csv_key]        
        if job_key in job.keys():
            val=job[job_key]
            job_title=f'{job_title}, {csv_key}: {val}'
            job_id=f'{job_id}_{val}'
            #e.g. /workdir/<subject>/<experiment>/... etc.
            job_workdir=job_workdir / val
    job['job_id'],job['job_title'],job['job_workdir']=job_id,job_title,job_workdir
    
    #Do not change the next two lines to correctly preserve the scan context
    job_scan_context=global_vars['g_input_mount_path'] 
    if env_type == 'jupyter': job_scan_context = job_scan_context / job['job_exp_label']
    job_scan_context=job_scan_context / 'SCANS'
    job['job_scan_context']=job_scan_context    
            
    return job

def add_job_step(job,step_command,step_title=None):
    step={ 'step_command': step_command }
    if step_title is not None:
        step['step_title']=step_title


def set_logger():
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)

#helper functions
def paths_to_str(x):
    if isinstance(x, Path): return str(x)
    if isinstance(x, dict): return {k: paths_to_str(v) for k, v in x.items()}
    if isinstance(x, list): return [paths_to_str(v) for v in x]
    return x

def resource_to_xnat(local_resource, xnat_scan_resource, xnat_project, xnat_subject, xnat_experiment, xnat_interface):
    return sync_resource_xnat(local_resource, xnat_scan_resource, xnat_project, xnat_subject, xnat_experiment, \
        level="experiment", create_hierarchy=True,xnat_interface=xnat_interface)

def get_xnat_interface(project):
    host,user,passw,=os.environ.get('XNAT_HOST'),os.environ.get('XNAT_USER'),os.environ.get('XNAT_PASS')
    xnat = pyxnat.Interface(server=host, user=user, password=passw)
    xnat.select.project(project)
    return xnat

def launch_cs_command(xnat_interface,project,subject,experiment,session,workflow_id,xnat_command_id,xnat_wrapper_name,verbose=False):
    host=os.environ.get('XNAT_HOST')
    url = f"{host}/xapi/projects/{project}/commands/{xnat_command_id}/wrappers/{xnat_wrapper_name}/launch"
    try:
        body = {
            "PROJECT": project,
            "SUBJECT": subject,
            "EXPERIMENT": experiment,
            "WORKFLOW_ID": workflow_id,
            "MICROENV": "NONE",
            "REPO_GIT": "NONE", 
            "session": session            
        }
        resp = xnat_interface._http.post(url, json=body)
        if verbose:
            print("Request body: ", body)
            print("Request url: ", url)
            print("Status:", resp.status_code)
            print("Reason:", resp.reason)
            print("URL:", resp.url)
            print("Headers:\n", resp.headers)
            print("Cookies:\n", resp.cookies)
            print("Elapsed:", resp.elapsed)
            print("Text:\n", resp.text)
        return resp.status_code in (200, 201, 202)
    except Exception as e:
        print(e)
        return False


def sync_resource_xnat(local_resource, resource_name, project, subject=None, 
       experiment=None, scan=None, upload=True, level="scan", XNAT_HOST=None, 
       username=None, password=None, xnat_interface=None, create_hierarchy=False, debug=False,include_dir_under_resource=False):
    '''
    Synchronizes a local file or directory with an XNAT resource at project, subject, experiment, or scan level.
    Supports both upload and download modes, with optional automatic creation of missing hierarchy elements. 
    Handles authentication via provided credentials or existing pyxnat Interface. 
    For uploads, directories are archived to ZIP (optionally preserving top-level folder) and extracted server-side. 
    For downloads, resources are retrieved as ZIP and extracted locally, with support for single-file or directory targets.     
    '''

    if XNAT_HOST is None: XNAT_HOST = os.environ.get("XNAT_HOST")
    if not XNAT_HOST: logging.error("Missing XNAT_HOST"); return 2
        
    if xnat_interface is None:
        if username is None or password is None:
            username = os.environ.get("XNAT_USER"); password = os.environ.get("XNAT_PASS")
        if not username or not password:
            logging.error("Missing XNAT credentials (user, password)"); return 2

    PROJECT_ID = str(project)
    params = { "inbody": "true", "PROJECT_ID": PROJECT_ID, "overwrite": "true" }

    subj = "" if subject is None else str(subject)
    exp = "" if experiment is None else str(experiment)
    scan_id = "" if scan is None else str(scan)

    levels = {'project': 1, 'subject': 2, 'experiment': 3, 'scan': 4}

    try:
        level_ord = levels[level]
        new_session = False
        
        if xnat_interface is not None: 
            xnat=xnat_interface
        else:
            xnat = Interface(server=str(XNAT_HOST), user=username, password=password)
            new_session = True

        project = xnat.select.project(PROJECT_ID)

        if level_ord >= levels['subject']:
            subj_obj = project.subject(subj); params['SUBJECT_ID'] = subj
            if not subj_obj.exists():
                if not upload or not create_hierarchy: 
                    logging.error(f"Subject {subj} does not exist")
                    return 1
                subj_obj.create()
                logging.debug(f"Created subject '{subj}'")

        if level_ord >= levels['experiment']:
            exp_obj = subj_obj.experiment(exp); params['EXPT_LABEL'] = exp
            if not exp_obj.exists():
                if not upload or not create_hierarchy: 
                    logging.error(f"Experiment '{exp}' does not exist")
                    return 2
                exp_obj.create(); logging.debug(f"Created experiment '{exp}'")

        if level_ord >= levels['scan']:
            scan_obj = exp_obj.scan(scan_id)
            if not scan_obj.exists():
                if not upload or not create_hierarchy: 
                    logging.error(f"Scan {scan_id} does not exist under experiment {exp}")
                    return 3
                scan_obj.create(); logging.debug(f"Created scan '{scan_id}'")
    
        if level == "project": 
            res_obj = project.resource(resource_name)            
            rlist=project.resources().get()
            logging.debug(f"Level project, resources:{rlist}")
            with xnat._http as s:
                    url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/resources?format=json"
                    r = s.get(url)
                    logging.debug(f"GET {url} -> {r.status_code}")
                    logging.debug(r.text[:2000])

        
        elif level == "subject": res_obj = subj_obj.resource(resource_name)
        elif level == "experiment": res_obj = exp_obj.resource(resource_name)
        elif level == "scan": res_obj = scan_obj.resource(resource_name)

        if upload and not res_obj.exists(): res_obj.create()
        if not upload and not res_obj.exists(): 
            logging.error(f"Resource '{resource_name}' does not exist at level {level}, res_obj:{res_obj}")
            return 4

    except Exception as e:
        print(e)
        logging.error(e)
        return 5
    finally:
        if new_session: xnat.disconnect()

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
                logging.debug(f"archiving directory {src} for upload")
                tmp_dir = Path(tempfile.mkdtemp(prefix="xnat_upload_"))
                tmp_zip_base = tmp_dir / "upload"
                tmp_zip = tmp_zip_base.with_suffix(".zip")
                if not include_dir_under_resource:
                    shutil.make_archive(str(tmp_zip_base), "zip", root_dir=src, base_dir=".")
                else:
                    shutil.make_archive(str(tmp_zip_base), "zip", root_dir=Path(src).parent, base_dir=Path(src).name)
                upload_items = [tmp_zip]

            with xnat._http as s:
                #main upload loop
                for item in upload_items:
                    with item.open("rb") as f:
                        is_zip_upload = (tmp_zip is not None and item == tmp_zip)
                        if is_zip_upload: 
                            params["extract"] = "true"
                            headers = {"Content-Type": "application/zip"}
                        else: 
                            headers = {"Content-Type": "application/octet-stream"}
                        logging.debug(f"uploading to {base_url}/{item.name}")
                        r = s.put(f"{base_url}/{item.name}", params=params, data=f, headers=headers)
                        logging.debug(f"Resource upload OK: {r.status_code} {base_url}")
                        if r.status_code >= 400: logging.error(f"Upload failed: {r.status_code} {base_url} {r.text[:1000]}"); return 2
                logging.info('Upload successful')
                return 0                            
        except Exception as e:
            print(e)
            logging.error(e)
            return 7
        finally:
            logging.debug("Cleaning up")
            if tmp_zip and tmp_zip.exists():
                try: tmp_zip.unlink()
                except Exception: pass
            if tmp_dir and tmp_dir.exists():
                try: shutil.rmtree(tmp_dir, ignore_errors=True)
                except Exception: pass
            if new_session: xnat.disconnect()
        logging.info('Upload successful')
        return 0
        
    else: #download
        dst = Path(local_resource)
        wants_dir = str(local_resource).endswith("/") or str(local_resource).endswith("\\") or (dst.exists() and dst.is_dir())
    
        tmp_zip_path = Path(tempfile.mkstemp(prefix="xnat_resource_", suffix=".zip")[1])
        try:
            with xnat._http as s:
                logging.debug(f"Downloading {base_url}")
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
                    logging.debug(f"Downloaded resource '{resource_name}' to directory {dst}")
                    logging.info('Download successful')
                    return 0
    
                if len(members) != 1: logging.error(f"Resource '{resource_name}' has {len(members)} files; local_resource is a file path"); return 2
                dst.parent.mkdir(parents=True, exist_ok=True)
                with z.open(members[0], "r") as fi, dst.open("wb") as fo: shutil.copyfileobj(fi, fo)
                logging.debug(f"Downloaded resource '{resource_name}' single file to {dst}")
                logging.info('Download successful')
                return 0
    
        except Exception as e:
            logging.error(e)
            return 10
        finally:
            logging.debug("Cleaning up")
            if tmp_zip_path.exists():
                try: tmp_zip_path.unlink()
                except Exception: pass
            if new_session: xnat.disconnect()
        logging.info('Download successful')
        return 0

##########################################################################################
# g_user_env_repo: the folder name stem of the user-supplied micromamba environment repository 
# stored as subfolder in the project 'ENVS' resource. Currently not supported due to XNAT limitations.
# if this is set to 'NONE', in container mode an environment built into Docker image (if any) will be used, 
# from the default location /opt/packages/user/user_env. 

##########################################################################################
# g_env_repo_dir: actual location of the the user-supplied micromamba environment repository. 
# 'jupyter' mode: if g_user_env_repo is set, this is set to 
# /data/projects/<PROJECT>/RESOURCES/ENVS/<g_user_env_repo>
# if g_user_env_repo='NONE', this can be set manually to point to the local environment repo.
# 'container' mode: default to '/opt/packages/user/user_env'. If g_user_env_repo is set, this location 
# is mounted from /data/projects/<PROJECT>/RESOURCES/ENVS/<g_user_env_repo>

##########################################################################################
# g_user_src_repo: folder name stem of the user source code/resources that will be used by the 
# workflow, stored as subfolder in the project 'SRC' resource.
# This can be 'NONE', then no resources directory will be available during runtime.

##########################################################################################
# g_alg_repo_dir: the actual location of the user source code/resources dir. 
# In 'jupyter' mode, default to '/data/projects/<PROJECT>/RESOURCES/SRC/<g_user_src_repo>'.
# Can be overridden if that repo does not exist.
# In 'container' mode, default to '/opt/packages/user/alg_repo'

##########################################################################################
# g_input_mount_path: the input folder with the list of input project experiments, 
# mounted from the XNAT archive.
# for 'jupyter', default to '/data/projects/<project>/experiments'
# for 'container', default to '/input

##########################################################################################
# g_local_workdir_path: local workdir where configuration files, logs, scripts, 
# and outputs will be written. 
# in 'jupyter' mode, this should point to a local writable dir.
# in 'container' mode, default to '/workdir'

#########################################################################################
# g_pymipl_dir: local directory with XNAT workflow generator scripts.
# 'jupyter' mode: any local dir where 'pymipl' library is cloned
# 'container' mode: auto-set to /opt/packages/pymipl

def init_global_vars(env_type, project, workflow_id, 
                     g_user_env_repo='NONE',
                     g_env_repo_dir=None,
                     g_user_src_repo='NONE',
                     g_alg_repo_dir=None,
                     g_input_mount_path=None,
                     g_local_workdir_path=None,
                     g_pymipl_dir=None
    ):
    gv={'g_project': project,'g_workflow_id': workflow_id}
    if env_type.lower()=='jupyter':
        gv['g_user_env_repo']=g_user_env_repo        
        if g_env_repo_dir is not None: 
            gv['g_env_repo_dir']=g_env_repo_dir
        elif g_user_env_repo != 'NONE':
            gv['g_env_repo_dir']=Path('/data/projects') / project /'RESOURCES/ENVS' / g_user_env_repo
        else:
            gv['g_env_repo_dir']='NONE'
        gv['g_user_src_repo']=g_user_src_repo
        if g_alg_repo_dir is not None: 
            gv['g_alg_repo_dir']=g_alg_repo_dir
        elif g_user_src_repo != 'NONE':
            gv['g_alg_repo_dir']=Path('/data/projects') / project /'RESOURCES/SRC' / g_user_src_repo
        else:
            gv['g_alg_repo_dir']='NONE'        
        gv['g_input_mount_path']= (
            g_input_mount_path if g_input_mount_path is not None
            else Path('/data/projects') / project / 'experiments'
        )
        if g_local_workdir_path is not None: 
            gv['g_local_workdir_path']=g_local_workdir_path 
        else:
            raise ValueError('g_local_workdir_path cannot be empty if env_type=jupyter')
        if g_pymipl_dir is not None: 
            gv['g_pymipl_dir']=g_pymipl_dir
        else:
            raise ValueError('g_pymipl_dir_path cannot be empty if env_type=jupyter')                            
    
    elif env_type.lower()=='container':
        gv['g_user_env_repo']=g_user_env_repo        
        gv['g_env_repo_dir']=Path('/opt/packages/user/env_repo')
        gv['g_user_src_repo']=g_user_src_repo
        gv['g_alg_repo_dir']=Path('/opt/packages/user/alg_repo')
        gv['g_input_mount_path']=Path('/input')
        gv['g_local_workdir_path']=Path('/workdir')
        gv['g_pymipl_dir']=Path('/opt/packages/pymipl')
    else:
        raise ValueError(f'Unknown environment type: {env_type}, expected jupyter or container')    
    return gv

def init_global_vars_bootstrap_image(global_vars,xnat_project):
    global_vars['g_input_mount_path']=Path('/input')
    global_vars['g_local_workdir_path']=Path('/workdir')
    global_vars['g_pymipl_dir']=Path('/opt/packages/pymipl')
    global_vars['g_alg_repo_dir']=Path('/opt/packages/user/alg_repo')
    global_vars['g_env_repo_dir']=Path('/opt/packages/user/env_repo')
    global_vars['g_project']=xnat_project

def workflow_to_batch(job_yaml, global_vars, output_batch_file, error_message=None):
    '''
    Generates an executable bash batch script from a workflow YAML definition.
    Merges job-level variables with provided globals and formats all placeholders.
    
    Writes script incrementally, adding shebang only if file does not exist.
    Prints job and step titles as annotated echo statements for traceability.
    
    For each step, resolves parameters and renders a single safe command.
    Wraps execution with error checks and exits on failure with message.
    
    Builds XNAT upload commands for specified files and directories per step.
    Validates existence (and non-empty dirs) before invoking upload helper.
    
    Uses positional arguments (project, subject, experiment, optional scan).
    Assumes target environment (e.g., pymipl) is already active at runtime.
    
    Ensures output directory exists and marks script as executable.    
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
                if error_message is None: error_message=f'Command failed, exiting'
                #insert command prefix if not pymipl command.
                f.write(f'cmd=({cmd})\n')
                f.write(f'echo Running command: {cmd}\n')
                f.write(f'if ! "${{cmd[@]}}"; then\n    echo "{error_message}"\n    exit 1\nfi\n')
            args = [
                f'--project "{s_job["g_project"]}"',
                f'--subject "{s_job["job_subject"]}"',
                f'--experiment "{s_job["job_exp_label"]}"'
            ]
            if 'job_scan_id' in s_job.keys(): args+=[f'--scan "{s_job["job_scan_id"]}"']
            #form upload command syntax
            upl = str(Path(s_job["g_pymipl_dir"]) / "sync-resource-with-xnat.py").format(**s_job)
            
            #write step commands to upload each specified local file to the resource.
            for res, files in (s.get("step_upload_files_to_resource") or {}).items():
                for p in files:
                    p = (str(p) if isinstance(p, Path) else str(p)).format(**s_job)
                    cmd = f'python "{upl}" {" ".join(args)} --source_loc "{p}" --resource_name "{res}"'
                    f.write(f'if ! [[ -f "{p}" ]]; then\n    echo "Missing file: {p}" \n    exit 1\nfi\n')
                    #note that there's no command prefix, because environment with pymipl will be expected to be active during execution.
                    f.write(f'cmd=({cmd})\n')
                    f.write(f'echo Running command: {cmd}\n')
                    f.write(f'if ! "${{cmd[@]}}"; then\n    printf \'%s\\n\' "command failed, exiting"\n    exit 1\nfi\n')     

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

