'''
Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2026, Computational Imaging Lab, Washington University School of Medicine

Redistribution and use in source and binary forms, for any purpose, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import sys, os, pydicom, argparse
import nibabel as nib
from nibabel.orientations import io_orientation, ornt_transform
import numpy as np
import ipywidgets as ipw
#from utils import write_rec_file

def sort_dcms_by_slice_pos(input_dicom_path,dcm_files,stop_before_pixels=True):
    '''
    Sort DICOMs from an input directory (assumed to contain a single study) according to slice position
    Output: a sorted list of DICOM datasets
    '''
    dcmss=[]
    for idx,dcm in enumerate(dcm_files):
        ds = pydicom.dcmread(os.path.join(input_dicom_path,dcm), stop_before_pixels)
        if idx==0:           
            #if 'ImagePositionPatient' in ds: sortTag='ImagePositionPatient'
            if 'SliceLocation' in ds: sortTag='SliceLocation'
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
                         newSeriesInstanceUID:str,newSeriesNumber:int,flipX:bool,flipY:bool,flipZ:bool,\
                         intensityScaling:bool=False):
    '''
    Main routine to read NIFTI and DICOM images, replace voxels and specified meta tags in DICOM dataset, 
    and write back synthetic DICOM.
    '''
    #load NIFTI image
    nii0=nib.load(input_nifti)
    
    #this line works around apparent bug in nib's Nifti1Header.set_dim_info function
    nii0.header.set_dim_info(None,None,None)
    
    dcm_in_files=next(os.walk(input_dcm))[2]
    numberOfDicomImages = len(dcm_in_files)
    
    
    #get a list of DICOM datasets, one dataset per slice
    dcm_in_sorted=sort_dcms_by_slice_pos(input_dcm,dcm_in_files,stop_before_pixels=False)

    ds0=dcm_in_sorted[0]['dataset']
    dcm_pixeldata_type=ds0.pixel_array.dtype

    #reorient NIfTI to match DICOM axes derived from ImageOrientationPatient
    if 'ImageOrientationPatient' not in ds0:
        print('ERROR: DICOM missing ImageOrientationPatient; cannot determine orientation.')
        return -1
    iop=np.array([float(x) for x in ds0.ImageOrientationPatient])
    row_dir=iop[:3]
    col_dir=iop[3:]
    slice_dir=np.cross(row_dir, col_dir)
    lps_to_ras=np.diag([-1.,-1.,1.])
    dcm_axes_ras=lps_to_ras @ np.column_stack([row_dir, col_dir, slice_dir])
    dcm_affine = np.eye(4)
    dcm_affine[:3, :3] = dcm_axes_ras
    target_ornt=io_orientation(dcm_affine)
    nii_ornt=io_orientation(nii0.affine)

    print(dcm_axes_ras)
    print(target_ornt)
    print(nii_ornt)

    transform=ornt_transform(nii_ornt, target_ornt)
    try:
        nii=nii0.as_reoriented(transform)
    except Exception as e:
        print(e)
        return -1
    print('NIfTI reoriented to match DICOM axes')

    #read voxel arrays
    dcm_in_voxels=voxel_array_from_sorted_dicoms(dcm_in_sorted)
    dcm_dtype=ds0.pixel_array.dtype
    dtype_info=np.iinfo(dcm_dtype)
    if intensityScaling:
        # scale NIfTI range linearly onto DICOM integer range;
        # encode the mapping in RescaleSlope/Intercept so readers recover real values.
        raw=nii.get_fdata()
        vmin,vmax=raw.min(),raw.max()
        slope=(vmax-vmin)/(dtype_info.max-dtype_info.min) if vmax!=vmin else 1.0
        intercept=vmin-dtype_info.min*slope
        nii_in_voxels=np.round((raw-intercept)/slope).astype(dcm_dtype)
    else:
        # cast NIfTI to DICOM dtype with explicit clipping; avoid float64
        # round-trip when NIfTI is already integer.
        raw=nii.get_data() if np.issubdtype(nii.header.get_data_dtype(),np.integer) else nii.get_fdata()
        clipped=np.clip(raw,dtype_info.min,dtype_info.max)
        if not np.array_equal(raw,clipped):
            print(f'WARNING: NIfTI values clipped to [{dtype_info.min},{dtype_info.max}] to fit DICOM dtype {dcm_dtype}')
        nii_in_voxels=clipped.astype(dcm_dtype)
        slope,intercept=None,None
    nii_in_voxels=nii_in_voxels.squeeze()

    if dcm_in_voxels.shape != nii_in_voxels.shape:
        print ('NIFTI and DICOM image shapes don\'t match!')
        print ('NIFTI shape:',nii_in_voxels.shape)
        print ('DICOM shape:',dcm_in_voxels.shape)
        return -1
    
    #pre-flip NIFTI voxels
    if flipX: nii_in_voxels=np.flip(nii_in_voxels,0)
    if flipY: nii_in_voxels=np.flip(nii_in_voxels,1)
    if flipZ: nii_in_voxels=np.flip(nii_in_voxels,2)

    #initialize optional metadata
    siUID=pydicom.uid.generate_uid() if newSeriesInstanceUID is None else newSeriesInstanceUID
    sDescr=ds0.SeriesDescription if newSeriesDescription is None else newSeriesDescription
    sNumber=ds0.SeriesNumber if newSeriesNumber is None else newSeriesNumber        

    #make output dcm dir
    try:
        os.mkdir(output_dcm)
    except OSError as error:
        print(error)      
    
    #cycle through input DICOM datasets, replace voxels and metadata, and save in output DICOM dir
    for i in range(len(dcm_in_sorted)):
        ds=dcm_in_sorted[i]['dataset']
        print(f'ImagePositionPatient: {ds.ImagePositionPatient}, SliceLocation: {ds.SliceLocation}')
        slice_data=np.transpose(nii_in_voxels[:,:,i])
        ds.PixelData=slice_data.tobytes()
        ds[0x0020, 0x000e].value=siUID
        ds[0x0008,0x103e].value=sDescr
        ds[0x0020, 0x0011].value=sNumber
        if intensityScaling:
            ds.RescaleSlope=slope
            ds.RescaleIntercept=intercept
            ds.RescaleType='US'
        ds.SmallestImagePixelValue=int(slice_data.min())
        ds.LargestImagePixelValue=int(slice_data.max())
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
    parser.add_argument("--intensity_scaling",metavar="<bool>",type=str,default='False',choices=['True','False'],
                        help='scale NIfTI intensities to DICOM integer range and encode mapping in RescaleSlope/Intercept [False]')

    return parser.parse_args()

if __name__ == "__main__":

    p = get_parser()
    print(p)
    #write_rec_file(p.output_dicom, infiles=[p.input_dicom,p.input_nifti])

    sys.exit (convert_nifti_to_dcm(p.input_dicom,p.input_nifti,p.output_dicom,p.series_description, \
                        p.series_uid,p.series_number,p.flip_x,p.flip_y,p.flip_z,\
                        p.intensity_scaling=='True'))
    
