#!/usr/bin/env python3
import sys
import argparse
import logging
from workflow_adapters import sync_resource_xnat

def main():
    ap = argparse.ArgumentParser(description="Upload or download a file/directory to/from an XNAT resource at project/subject/experiment/scan level")
    ap.add_argument("--level", default="scan", choices=["project", "subject", "experiment", "scan"])
    ap.add_argument("--project", required=True)
    ap.add_argument("--subject")
    ap.add_argument("--experiment")
    ap.add_argument("--scan")
    ap.add_argument("--local_resource", required=True, help="Local file or directory to upload or download")
    ap.add_argument("--remote_resource", required=True)
    ap.add_argument("--upload", type=int, choices=[0, 1], default=1, help="1=upload local_resource to XNAT; 0=download resource to local_resource")
    ap.add_argument("--xnat_host")
    ap.add_argument("--user")
    ap.add_argument("--password")
    ap.add_argument("--create_hierarchy", type=int, choices=[0, 1], default=0)
    ap.add_argument("--logfile")
    args = ap.parse_args()

    if args.level == "subject" and not args.subject: ap.error("--subject is required when --level subject")
    if args.level == "experiment" and (not args.subject or not args.experiment): ap.error("--subject and --experiment are required when --level experiment")
    if args.level == "scan" and (not args.subject or not args.experiment or not args.scan): ap.error("--subject, --experiment, and --scan are required when --level scan")

    if args.logfile: logging.basicConfig(filename=args.logfile + ".log", encoding="utf-8", filemode="a", format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M")
    else: logging.basicConfig(stream=sys.stdout, level=logging.INFO, format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M")

    return sync_resource_xnat(args.local_resource, args.remote_resource, args.project, args.subject, args.experiment, args.scan, bool(args.upload), args.level, args.xnat_host, args.user, args.password, bool(args.create_hierarchy))

if __name__ == "__main__":
    sys.exit(main())