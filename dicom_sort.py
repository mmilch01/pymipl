import os, argparse, json, csv
from pathlib import Path
import pydicom

def get_dicom_features(file,d):
    try:
        ds=pydicom.dcmread(file,stop_before_pixels=True)
    except Exception as e:
        return False
    d['PatName']=str(ds.PatientName)
    try:
        d['SeriesDescription']=str(ds.SeriesDescription)
    except Exception as e:
        d['SeriesDescription']=None
        
    if ds.SOPClassUID=='1.2.840.10008.5.1.4.1.1.66.4':
        d['SOPClass']='Seg'
        try: 
            d['ReferencedSeriesInstanceUID']=ds.ReferencedSeriesSequence[0].SeriesInstanceUID
        except Exception as e:
            d['ReferencedSeriesInstanceUID']=None

    elif ds.SOPClassUID=='1.2.840.10008.5.1.4.1.1.481.3':
        d['SOPClass']='RTStruct'
        try:
            d['FrameOfReferenceUID']=ds.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID
        except Exception as e:
            d['FrameOfReferenceUID']=None
        try:
            d['ReferencedSeriesInstanceUID']=ds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID
        except Exception as e:
            d['ReferencedSeriesInstanceUID']=None

    elif ds.SOPClassUID=='1.2.840.10008.5.1.4.1.1.2':
        d['SOPClass']='CTImageStorage'
        try:
            d['SeriesInstanceUID']=ds.SeriesInstanceUID
        except Exception as e:
            d['SeriesInstanceUID']=None
    elif ds.SOPClassUID=='1.2.840.10008.5.1.4.1.1.4':
        d['SOPClass']='MRImageStorage'
        try:
            d['SeriesInstanceUID']=ds.SeriesInstanceUID
        except Exception as e:
            d['SeriesInstanceUID']=None        
    elif ds.SOPClassUID=='1.2.840.10008.5.1.4.1.1.130':
        d['SOPClass']='PETImageStorage'
    else:
        d['SOPClass']=ds.SOPClassUID
    return True
    

def process_subdir(d:dict,root:Path):
    is_first_file=True
    if 'children' not in d.keys(): d['children']=[]
    for entry in Path(root / Path(d['path'])).iterdir():
        if entry.is_dir():
            entry_dict={}
            entry_dict['path']=entry.relative_to(root).as_posix()
            #print(entry.as_posix())
            d['children']+=[entry_dict]
            #entry_dict['parent']=d
            level=d['level']+1
            entry_dict['level']=level
            if level==1: print(entry.as_posix())
            process_subdir(entry_dict,root)
            
        elif entry.is_file():
            file_dict={}
            if not get_dicom_features(entry, file_dict) or not is_first_file: continue
            is_first_file=False            
            file_dict['path']=entry.relative_to(root).as_posix()
            #file_dict['parent']=d
            file_dict['level']=d['level']+1
            d['children']+=[file_dict]

def analyze_dir(dir, save_to_file=None):
    root=Path(dir).absolute()
    d={'path':root.as_posix(),'parent':None, 'level': 0}
    process_subdir(d,root)
    if save_to_file:
        with open(save_to_file,'w') as txt:
            json.dump(d,txt,indent=2)
    return d


def get_subject(file_entry):
    '''
        Subject field location may be project-specific. For now, return 'PatName' field (equals subject for RIDER project).
    '''
    return file_entry['PatName']
    
def save_scanlist(subjects,csv_file):
    """
    Save structural scans and their derived segmentations to CSV.

    One row per structural scan.

    Segmentation columns are assigned globally by exact
    (SeriesDescription, SOPClass) match:
      SegScan1, SegScan1_SeriesDescription, SegScan1_SOPClass, ...
    """
    
    # Global segmentation slot registry:
    #   (SeriesDescription, SOPClass) -> slot number
    rows = []
    seg_scans = {}
    n_seg = 0
    
    for subject_id, subject_data in subjects.items():
        xnat_subject = subject_data.get("id", subject_id)

        for struct_uid, struct_entry in subject_data.get("struct_uids", {}).items():
            structural = struct_entry.get("structural", {}) or {}
            segs = struct_entry.get("segmentations", []) or []

            struct_path = structural.get("path", "") or ""
            struct_series_desc = structural.get("SeriesDescription", "")
            struct_series_desc = "" if struct_series_desc is None else struct_series_desc

            # Experiment = first path component
            path_parts = struct_path.split("/")
            xnat_experiment = path_parts[0] if path_parts else ""

            # Structural scan number = token after /SCANS/
            struct_scan = ""
            if "/SCANS/" in struct_path:
                after_scans = struct_path.split("/SCANS/", 1)[1]
                struct_scan = after_scans.split("/", 1)[0]

            row = {
                "Subject": xnat_subject,
                "Experiment": xnat_experiment,
                "StructScan": struct_scan,
                "StructScanSerDesc": struct_series_desc,
            }

            for seg in segs:
                sop_class = seg.get("SOPClass", "")
                #if sop_class not in ("RTStruct", "Seg"): continue

                seg_path = seg.get("path", "") or ""
                seg_series_desc = seg.get("SeriesDescription", "")
                seg_series_desc = "" if seg_series_desc is None else seg_series_desc

                seg_scan = ""
                if "/SCANS/" in seg_path:
                    after_scans = seg_path.split("/SCANS/", 1)[1]
                    seg_scan = after_scans.split("/", 1)[0]

                seg_key = (seg_series_desc, sop_class)

                if seg_key in seg_scans:
                    seg_num = seg_scans[seg_key]
                else:
                    n_seg += 1
                    seg_num = n_seg
                    seg_scans[seg_key] = seg_num

                row[f"SegScan{seg_num}"] = seg_scan
                row[f"SegScan{seg_num}_SerDesc"] = seg_series_desc
                row[f"SegScan{seg_num}_SOPClass"] = sop_class

            rows.append(row)

    fieldnames = [
        "Subject",
        "Experiment",
        "StructScan",
        "StructScanSerDesc",
    ]

    for i in range(1, n_seg + 1):
        fieldnames.extend([
            f"SegScan{i}",
            f"SegScan{i}_SerDesc",
            f"SegScan{i}_SOPClass",
        ])

    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            
    return rows

def reindex_to_structurals_and_segs(d, save_to_file=None,csv_file=None):
    '''
    reindex analyzed dir to structurals and corresponding segmentations.
    '''
    #first pass: build array of subjects, to each subject add array of structurals.
    #list of subjects. each subject is dict with serinstUIDs as key.
    subjects={}
    struct_scans={}
    
    #pass 1, populate subjects and struct_scans    
    for exp in d['children']:
        for subd1 in exp['children']:
            if Path(subd1['path']).stem == 'SCANS': 
                for scan in subd1['children']:
                    for scan_resource in scan['children']:
                        try:                     
                            #if Path(scan_resource['path']).name.lower() != 'dicom': continue
                            if len(scan_resource['children']) < 1: continue
                            file_entry=scan_resource['children'][0]
                            sopclass = file_entry['SOPClass']
                            uid=None
                            
                            if sopclass in ['CTImageStorage','MRImageStorage']:
                                uid=file_entry['SeriesInstanceUID']
                                if uid not in struct_scans:
                                    struct_scans[uid]={'structural':file_entry, 'segmentations': []}
                                else:
                                    struct_scans[uid]['structural']=file_entry                                
                                    
                            elif sopclass in ['RTStruct','Seg']:
                                uid=file_entry['ReferencedSeriesInstanceUID']
                                if uid not in struct_scans:
                                    struct_scans[uid]={'structural':None,'segmentations':[]}
                                struct_scans[uid]['segmentations']+=[file_entry]
                                
                            if uid is not None:
                                subject=get_subject(file_entry)                      
                                if subject not in subjects:
                                    subjects[subject]={'id':subject,'struct_uids':{}}
                                if uid not in subjects[subject]['struct_uids']:
                                    subjects[subject]['struct_uids'][uid]={}
                        except Exception as e:
                            print (f'WARNING: skipping scan {file_entry}: {e}')
                            continue
                        
    #pass 2, link subjects and structural scans.
    for uid, entry in struct_scans.items():
        for subject,subj_info in subjects.items():
            subj_uid_dict=subj_info['struct_uids']
            if uid in subj_uid_dict:
                subj_uid_dict[uid]=struct_scans[uid]
    rows=None        
    if save_to_file:
        with open(save_to_file,'w') as txt:
            json.dump(subjects,txt,indent=2)
    if csv_file:
        rows=save_scanlist(subjects,csv_file)
    return subjects, rows
'''
parses a directory with project data containing structural DICOM's, DICOMRTs and DICOMSeg's. 
'''
if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("root_dir",type=str,help="Directory to search")
    args=parser.parse_args()
    d=analyze_dir(args.root_dir,'tree.json')
    
