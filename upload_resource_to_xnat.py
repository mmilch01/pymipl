#!/usr/bin/env python3
import os
import sys
import logging
import shutil
import tempfile
from pathlib import Path
import argparse

import requests
from pyxnat import Interface


def upload_resource(xnat_project, subject, experiment, scan, source_loc, resource_name, level="scan", XNAT_HOST=None, username=None, password=None, create_hierarchy=False):
    
    if XNAT_HOST is None:
        XNAT_HOST = os.environ.get("XNAT_HOST")
    if username is None or password is None:
        username = os.environ.get("XNAT_USER")
        password = os.environ.get("XNAT_PASS")

    if not XNAT_HOST or not username or not password:
        logging.error("Missing XNAT credentials (XNAT_HOST, user, password)")
        return 2

    PROJECT_ID = str(xnat_project)
    params = {
        "inbody": "true", 
        "PROJECT_ID": PROJECT_ID, 
        "overwrite": "true"
    }
    
    subj = "" if subject is None else str(subject)
    exp = "" if experiment is None else str(experiment)
    scan_id = "" if scan is None else str(scan)

    levels={'project': 1,'subject':2,'experiment':3,'scan':4}        
    
    try:
        level_ord=levels[level]
        xnat = Interface(server=str(XNAT_HOST), user=username, password=password)
        project = xnat.select.project(PROJECT_ID)
        #create hierarchy if needed
        if level_ord >= levels['subject']:
            subj_obj=project.subject(subj)
            params['SUBJECT_ID']=subj
            if not subj_obj.exists(): 
                if not create_hierarchy: logging.error(f"Subject {subj} does not exist"); return -1
                else: subj_obj.create(); logging.info(f"Created subject '{subj}'")
        if level_ord >= levels['experiment']:
            exp_obj=subj_obj.experiment(exp)
            params['EXPT_LABEL']=exp
            if not exp_obj.exists():
                if not create_hierarchy: logging.error(f"Experiment '{exp}' does not exist"); return -1
                else: exp_obj.create(); logging.info(f"Created experiment '{exp}'")
        if level_ord >=levels['scan']:
            scan_obj=exp_obj.scan(scan_id)
            if not scan_obj.exists():
                if not create_hierarchy: logging.error(f"Scan {scan_id} does not exist under experiment {exp}"); return -1
                else: scan_obj.create(); logging.info(f"Created scan '{scan_id}'")

        if level == "project": res_obj=project.resource(resource_name)
        elif level == "subject": res_obj=subj_obj.resource(resource_name)
        elif level == "experiment": res_obj=exp_obj.resource(resource_name)
        elif level == "scan": res_obj=scan_obj.resource(resource_name)

        #create the resource.
        if not res_obj.exists(): res_obj.create()
        
        src = Path(source_loc)
        if not src.exists(): logging.error(f"Source path does not exist: {source_loc}"); return 2

        upload_items = []
        tmp_dir, tmp_zip = None,None

        if src.is_file(): upload_items = [src]
        else:
            tmp_dir = Path(tempfile.mkdtemp(prefix="xnat_upload_"))
            base = tmp_dir / "payload"
            shutil.copytree(src, base)
            tmp_zip = Path(str(base) + ".zip")
            shutil.make_archive(str(base), "zip", root_dir=tmp_dir, base_dir="payload")
            upload_items = [tmp_zip]

    except Exception as e:
        logging.error(e)
        return 2

    if level == "project":
        base_url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/resources/{resource_name}/files"
    if level == "subject":
        base_url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/subjects/{subj}/resources/{resource_name}/files"
    if level == "experiment":
        base_url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/subjects/{subj}/experiments/{exp}/resources/{resource_name}/files"
    if level == "scan":
        base_url = f"{XNAT_HOST}/data/archive/projects/{PROJECT_ID}/subjects/{subj}/experiments/{exp}/scans/{scan_id}/resources/{resource_name}/files"

    try:
        with requests.Session() as s:
            s.auth = (username, password)
            for item in upload_items:
                with item.open("rb") as f:
                    is_zip_upload = (tmp_zip is not None and item == tmp_zip)
                    if is_zip_upload:
                        params["extract"]="true"
                        headers={"Content-Type": "application/zip"}
                    else:
                        headers={"Content-Type": "application/octet-stream"}
                    
                    r = s.put(f"{base_url}/{item.name}", params=params, data=f, headers=headers)
                    #print("Resource upload OK:", r.status_code, r.text[:500], "...")
                    logging.info(f"Resource upload OK: {r.status_code} {base_url} {r.text[:500]}")
                    if r.status_code >= 400:
                        logging.error(f"Upload failed: {r.status_code} {base_url} {r.text[:1000]} ")
                        return 2
    finally:
        if tmp_zip and tmp_zip.exists():
            try:
                tmp_zip.unlink()
            except Exception:
                pass
        if tmp_dir and tmp_dir.exists():
            try:
                shutil.rmtree(tmp_dir, ignore_errors=True)
            except Exception:
                pass

    return 0


def main():
    ap = argparse.ArgumentParser(description="Upload a file or directory to an XNAT resource at project/subject/experiment/scan level")
    ap.add_argument("--level", default="scan", choices=["project", "subject", "experiment", "scan"])
    ap.add_argument("--xnat_project", required=True)
    ap.add_argument("--subject")
    ap.add_argument("--experiment")
    ap.add_argument("--scan")
    ap.add_argument("--source_loc", required=True)
    ap.add_argument("--resource_name", required=True)
    ap.add_argument("--xnat_host")
    ap.add_argument("--user")
    ap.add_argument("--password")
    ap.add_argument("--create_hierarchy", type=int, choices=[0, 1], default=0)
    ap.add_argument("--logfile")
    args = ap.parse_args()

    if args.level == "subject" and not args.subject:
        ap.error("--subject is required when --level subject")
    if args.level == "experiment" and (not args.subject or not args.experiment):
        ap.error("--subject and --experiment are required when --level experiment")
    if args.level == "scan" and (not args.subject or not args.experiment or not args.scan):
        ap.error("--subject, --experiment, and --scan are required when --level scan")

    if args.logfile:
        logging.basicConfig(filename=args.logfile + ".log", encoding="utf-8", filemode="a", format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M")
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M")

    return upload_resource(args.xnat_project, args.subject, args.experiment, args.scan, args.source_loc, args.resource_name, args.level, args.xnat_host, args.user, args.password, bool(args.create_hierarchy))


if __name__ == "__main__":
    sys.exit(main())
