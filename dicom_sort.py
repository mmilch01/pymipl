import os, argparse, json
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
            print(entry.as_posix())
            d['children']+=[entry_dict]
            #entry_dict['parent']=d
            entry_dict['level']=d['level']+1
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

'''
parses a directory with project data containing structural DICOM's, DICOMRTs and DICOMSeg's. 
'''
if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("root_dir",type=str,help="Directory to search")
    args=parser.parse_args()
    d=analyze_dir(args.root_dir,'tree.json')
    
