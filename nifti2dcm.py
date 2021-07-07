'''
Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2021, Computational Imaging Lab, Washington University School of Medicine

Redistribution and use in source and binary forms, for any purpose, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import os, pydicom, argparse
import nibabel as nib
import numpy as np
import ipywidgets as ipw
from utils import write_rec_file

def sort_dcms_by_slice_pos(input_dicom_path,dcm_files,stop_before_pixels=True):
    '''
    Sort DICOMs from an input directory (assumed to contain a single study) according to slice position
    Output: a sorted list of DICOM datasets
    '''
    dcmss=[]
    for idx,dcm in enumerate(dcm_files):
        ds = pydicom.dcmread(os.path.join(input_dicom_path,dcm), stop_before_pixels)
        if idx==0:           
            if 'ImagePositionPatient' in ds: sortTag='ImagePositionPatient'
            elif 'SliceLocation' in ds: sortTag='SliceLocation'
            else: return None
        if not sortTag in ds: return None
        if sortTag=='ImagePositionPatient': z=ds.ImagePositionPatient[2]
        else: z=ds.SliceLocation
        dcmss+=[dict(file=dcm,dataset=ds,z=z)]
    return sorted(dcmss, key=lambda dcms: dcms['z'])

def voxel_array_from_sorted_dicoms(dicomsSorted):
    '''
    extract the 3D voxel array from a list of sorted DICOM objects.
    '''
    ind=0
    
    if len(dicomsSorted) < 1: return None
    ds0=dicomsSorted[0]['dataset']
    
    imwidth,imheight,imdepth=ds0.Rows,ds0.Columns,len(dicomsSorted)
    pixeldata_type=ds0.pixel_array.dtype
    voxels=np.zeros([imwidth,imheight,imdepth],dtype=pixeldata_type)
    
    for i in range(len(dicomsSorted)):
        voxels[:,:,i]=np.transpose(dicomsSorted[i]['dataset'].pixel_array)
        
    return voxels

def convert_nifti_to_dcm(input_dcm:str, input_nifti:str, output_dcm:str, newSeriesDescription:str,\
                         newSeriesInstanceUID:str,newSeriesNumber:int,flipX:bool,flipY:bool,flipZ:bool):
    '''
    Main routine to read NIFTI and DICOM images, replace voxels and specified meta tags in DICOM dataset, 
    and write back synthetic DICOM.
    '''
    #load NIFTI image
    nii=nib.load(input_nifti)
    dcm_in_files=next(os.walk(input_dcm))[2]
    numberOfDicomImages = len(dcm_in_files)
    
    #get a list of DICOM datasets, one dataset per slice
    dcm_in_sorted=sort_dcms_by_slice_pos(input_dcm,dcm_in_files,stop_before_pixels=False)

    ds0=dcm_in_sorted[0]['dataset']
    dcm_pixeldata_type=ds0.pixel_array.dtype

    #read voxel arrays
    dcm_in_voxels=voxel_array_from_sorted_dicoms(dcm_in_sorted)
    nii_in_voxels=nii.get_fdata().astype(ds0.pixel_array.dtype)
    if dcm_in_voxels.shape != nii_in_voxels.shape:
        print ('NIFTI and DICOM image shapes don\'t match!')    
    
    #pre-flip NIFTI voxels
    if flipX: nii_in_voxels=np.flip(nii_in_voxels,0)
    if flipY: nii_in_voxels=np.flip(nii_in_voxels,1)
    if flipZ: nii_in_voxels=np.flip(nii_in_voxels,2)

    #initialize optional metadata
    siUID=pydicom.uid.generate_uid() if newSeriesInstanceUID is None else newSeriesInstanceUID
    sDescr=ds.SeriesDescription if newSeriesDescription is None else newSeriesDescription
    sNumber=ds.SeriesNumber if newSeriesNumber is None else newSeriesNumber        

    #make output dcm dir
    try:
        os.mkdir(output_dcm)
    except OSError as error:
        print(error)      
    
    #cycle through input DICOM datasets, replace voxels and metadata, and save in output DICOM dir
    for i in range(len(dcm_in_sorted)):
        ds=dcm_in_sorted[i]['dataset']
        ds.PixelData=np.transpose(nii_in_voxels[:,:,i]).tobytes()
        ds[0x0020, 0x000e].value=siUID
        ds[0x0008,0x103e].value=sDescr
        ds[0x0020, 0x0011].value=sNumber
        ds.save_as(output_dcm+'/'+str(i)+'.dcm')    


def get_parser():
    """
    Parse input arguments.
    """
    parser = argparse.ArgumentParser(description='Replace voxels in a DICOM image with those from a NIFTI image')

    # Positional arguments.
    parser.add_argument("input_dicom", help="path to input DICOM dir")
    parser.add_argument("input_nifti", help="path to input NIFTI image")    
    parser.add_argument("output_dicom", help="path to output DICOM dir")
    
    parser.add_argument("--series_description",metavar="<string>",type=str,default=None,help='new series description [keep original]')
    parser.add_argument("--series_uid",metavar="<string>",type=str,default=None,help='new series instance uid [auto-generated]')
    
    parser.add_argument("--series_number",metavar="<string>",type=str,default=None,help='new series number [keep original]')
    parser.add_argument("--flip_x",action="store_true",default=False,help='flip X axis')
    parser.add_argument("--flip_y",action="store_true",default=False,help='flip Y axis')    
    parser.add_argument("--flip_z",action="store_true",default=False,help='flip Z axis')
    
    return parser.parse_args() 

if __name__ == "__main__":
    
    p = get_parser()
    print(p)
    
    convert_nifti_to_dcm(p.input_dicom,p.input_nifti,p.output_dicom,p.series_description, \
                        p.series_uid,p.series_number,p.flip_x,p.flip_y,p.flip_z)
    write_rec_file(p.output_dicom,infiles=[p.input_dicom,p.input_nifti])
