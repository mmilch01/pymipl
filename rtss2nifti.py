'''
Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2021, Computational Imaging Lab, School of Medicine, Washington University in Saint Louis

Redistribution and use in source and binary forms, for any purpose, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import re,random, nibabel as nib, argparse, numpy as np, os, json
from datetime import datetime
from skimage import measure
from PIL import Image,ImageDraw

import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from nibabel.nifti1 import Nifti1Image,Nifti1Header
import nibabel.nifti1

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

def get_rasterized_poly_slice(poly2d, imwid, imht):
    '''
    Create a binary mask from a closed 2D polygon.
    Input: 2d polygon
    output: 2d image
    '''
    img=Image.new('L',(imwid,imht),0)
    ImageDraw.Draw(img).polygon(poly2d,outline=1,fill=1)
    return np.transpose(np.array(img))

def pts2poly(pts, origin, mm2vox_size):
    '''
    Reformat an array of polygon points in RTSS format into a numeric 2D points array    
    Input:     
    pts: array of 3D points in RTSS space, string encoded e.g. ['-11.2', '-75.2', '4.3','3','5','8']
    origin: numpy array of size 3, origin in RTSS space
    voxel_size: numpy array of size 3, number of voxels in 1mm

    Output: a tuple of array of 2D points in image raster space and an integer z coordinate.
    ''' 
    npts=len(pts) //  3
    if npts < 1: return None,None
    poly=np.zeros([2*npts])
    ind=0
    origin2d=origin[:2]
    mm2vox_size2d=mm2vox_size[:2]
    
    z=int(round((pts[2] - origin[2])*mm2vox_size[2]))
    
    for i in range(0,len(pts),3):
        pt=np.rint( np.multiply( np.array([float(pts[i]),float(pts[i+1])]) - origin2d , mm2vox_size2d ) )
        poly[ind:ind+2]=pt
        ind+=2
    
    return list(poly),z

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

def rtss_to_nifti(input_rtstruct_dicom:str, input_structural_dicom:str,output_rtss_nii:str,
                  output_struct_nii:str, exclude_labels:list, write_one_roi_per_file:bool):
    
    '''
    Convert RTSTRUCT and structural DICOM to a NIFTI mask.    
    '''
    
    #1. read the structural image.
    dicomFiles = next(os.walk(input_structural_dicom))[2]
    numberOfDicomImages = len(dicomFiles)
    dicomsSorted=sort_dcms_by_slice_pos(input_structural_dicom,dicomFiles,stop_before_pixels=False)
    ds_struct=dicomsSorted[0]['dataset']
    struct_voxels=voxel_array_from_sorted_dicoms(dicomsSorted)

    xPixelSize,yPixelSize=ds_struct.PixelSpacing[0],ds_struct.PixelSpacing[1]
    xyPixelSize=0.5*(xPixelSize+yPixelSize)
    #for now; remember to update with DistanceBetweenSlices
    zPixelSize=ds_struct.SliceThickness

    imwidth,imheight,imdepth=ds_struct.Rows,ds_struct.Columns,len(dicomsSorted)
    voxel_vol_mm3=xPixelSize*yPixelSize*zPixelSize
    

    if not write_one_roi_per_file:
        rtss_voxels=[np.zeros([imwidth,imheight,imdepth],dtype=np.uint16)]
    else: 
        rtss_voxels=[]

    print('Voxel size: {} by {} by {} mm'.format(xPixelSize,yPixelSize,zPixelSize))

    # Find position of first slice
    patientPosition = ds_struct.ImagePositionPatient    
    patientStartingZ = dicomsSorted[0]['z']

    mm_vox_size=np.array([1./xPixelSize,1./yPixelSize,1./zPixelSize])
    origin=np.array([patientPosition[0],patientPosition[1],patientStartingZ])

    print('Patient position is ', patientPosition[:2])
    print('First slice at ', patientStartingZ)

    #2. read the RTSTRUCT. 
    ds_rtss=pydicom.dcmread(input_rtstruct_dicom, stop_before_pixels=False)
    if [0x3006,0x0020] not in ds_rtss: 
        raise ValueError('Cannot find (0x3006,0020) StructureSetROISequence tag in RTSS file')
    
    structure_set_roi_sequence=ssrs=ds_rtss[(0x3006,0x0020)]._value
    roi_contour_sequence=rcs=ds_rtss[(0x3006,0x0039)]

    n=len(ssrs)
    print ('Found {} structures'.format(n))    

    roi_list=[]
    current_region_code=1
    roi_file_index=0
    
    for ind in range(n):
        ss=ssrs[ind]
        roi_number=ss[(0x3006,0x0022)].value
        roi_name=re.sub(r'\W+', '', ss[(0x3006,0x0026)].value)
        print('ROI number {}, name {}'.format(roi_number,roi_name))
        
        if roi_name.lower() in exclude_labels: 
            print('skipping structure',roi_name)
            continue

        roi_contour=None
        
        #now find the matching ROI contour.
        for rc in roi_contour_sequence:
            #print ('comparing',rc[0x3006, 0x0084]._value,'and', roi_number)
            if rc[0x3006, 0x0084]._value==roi_number: 
                roi_contour=rc
                break

        if not roi_contour: 
            print('WARNING: no matching roi contour sequence for this ROI')
            continue

        display_color = roi_contour[0x3006,0x002a] if (0x3006,0x002a) in roi_contour else None    
        print('display color:',display_color)

        if (0x3006,0x0040) in roi_contour:
            contour_sequence=roi_contour[0x3006,0x0040]

        else:
            print('WARNING: no matching contour sequence for this ROI')
            continue
        nContours=len(contour_sequence._value)

        #write rasterized data to image array.
        nptsAcc=0
        
        
        if write_one_roi_per_file:
            new_roi_voxels=np.zeros([imwidth,imheight,imdepth],dtype=np.uint16)
            rtss_voxels.append(new_roi_voxels)
            current_voxels=new_roi_voxels
        else:
            current_voxels=rtss_voxels[0]
            
        for k in range(nContours):
            contour=contour_sequence[k]
            npts=contour.NumberOfContourPoints
            nptsAcc+=npts
            print('Contour {}, number of points: {}'.format(k+1,npts))
            pts=contour.ContourData
            poly2d,z=pts2poly(pts,origin,mm_vox_size)

            imslice=current_region_code*get_rasterized_poly_slice(poly2d,imwidth,imheight)
            
            rtss_voxels[roi_file_index][:,:,z]=imslice
             
        #get region properties.
        try:
            p=measure.regionprops((current_voxels==current_region_code).astype(int))[0]
            vol_mm3=p.area*voxel_vol_mm3
        except:
            vol_mm3=0
            
        print('volume:',vol_mm3,'mm3')

        v=display_color.value
        color='0x{:02X}{:02X}{:02X}'.format(int(v[0]),int(v[1]),int(v[2]))
        
        if not write_one_roi_per_file:
            out_file=output_rtss_nii
        else: 
            out_file=output_rtss_nii+'_'+roi_name
            
        roi_descriptor=dict(roi_number=roi_number,
                            roi_name=roi_name,
                            display_color=color,
                            num_contours=nContours,
                            points_in_all_contours=nptsAcc,
                            intensity_value=current_region_code,
                            out_file_root=out_file,
                            volume_mm3=vol_mm3
                           )
        
        roi_list.append(roi_descriptor)
        
        #update appropriate counters
        if not write_one_roi_per_file:
            current_region_code+=1
        else:
            roi_file_index+=1
            
    
    #create and save nifti images.
    #flip axes to orient from DICOM (LPS) to RAS space

    flips=np.array([-1.,-1.,1.,1.])
    nifti_affine=np.diag(np.array([xPixelSize,yPixelSize,zPixelSize,1])*flips)
    nifti_affine[:3,3]=origin*flips[:-1]

    nifti_image_struct=Nifti1Image(struct_voxels,nifti_affine)    
        
    if not write_one_roi_per_file:
        nifti_image_roi=Nifti1Image(rtss_voxels[0],nifti_affine)
        print ('writing',output_rtss_nii)
        nibabel.nifti1.save(nifti_image_roi,output_rtss_nii)
    else:
        for i in range(len(rtss_voxels)):
            nifti_image_roi=Nifti1Image(rtss_voxels[i],nifti_affine)
            outfile=roi_list[i]['out_file_root']
            print('writing',outfile)
            nibabel.nifti1.save(nifti_image_roi,outfile)

    print('writing',output_rtss_nii+'.json')
    with open(output_rtss_nii+'.json', 'w') as fout:
        json.dump(roi_list,fout)

    print('writing',output_struct_nii)
    nibabel.nifti1.save(nifti_image_struct,output_struct_nii)
    print('done')    
    
def get_parser():
    """
    Parse input arguments.
    """
    parser = argparse.ArgumentParser(description='Convert DICOM RT structure images to NIFTI')

    # Positional arguments.
    parser.add_argument("in_rtss", help="Input DICOM RTSTRUCT file")
    parser.add_argument("in_struct_dir", help="Input structural DICOM directory")
    parser.add_argument("out_roi_mask", help="Output ROI mask file root")
    parser.add_argument("--out_struct", metavar="<string>",type=str,default=None, 
                        help="Output structural image root [output_nifti_rtss+struct.nii]")
    parser.add_argument("--exclude_labels", metavar="<string>",type=str,default=None,
                        help="Comma separated list of ROI labels to exclude, case insensitive [None]")
    parser.add_argument("--separate_masks", action="store_true", default=False, help="write each ROI mask in a separate file [False]")

    return parser.parse_args()
    
if __name__ == "__main__":
    p = get_parser()
    structural=p.out_roi_mask+'_struct.nii' if p.out_struct is None else p.out_struct
    
    exc_labels=[] if p.exclude_labels is None else p.exclude_labels.split(',')
    
    for i in range(len(exc_labels)):
        exc_labels[i]=exc_labels[i].lower()
        
    rtss_to_nifti(p.in_rtss, p.in_struct_dir,p.out_roi_mask,
                  structural,exc_labels,p.separate_masks)
    
    write_rec_file(p.out_roi_mask,main_extension='nii',infiles=[p.in_rtss,p.in_struct_dir])
    write_rec_file(structural,main_extension='nii',infiles=[p.in_rtss,p.in_struct_dir])
    
    print('done')
           