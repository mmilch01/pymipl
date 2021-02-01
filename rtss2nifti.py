import re,random, nibabel as nib, argparse, numpy as np, os, json
from datetime import datetime
from skimage import measure
from PIL import Image,ImageDraw

import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from nibabel.nifti1 import Nifti1Image,Nifti1Header
import nibabel.nifti1

def sort_dcms_by_slice_pos(input_dicom_path,dcm_files,stop_before_pixels=True):
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
    Input: 2d polygon
    output: 2d image
    '''
    img=Image.new('L',(imwid,imht),0)
    ImageDraw.Draw(img).polygon(poly2d,outline=1,fill=1)
    return np.transpose(np.array(img))

def pts2poly(pts, origin, mm2vox_size):
    '''
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
                  output_struct_nii:str,exclude_labels:list)
    
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

    rtss_voxels=np.zeros([imwidth,imheight,imdepth],dtype=np.uint16)

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

for ind in range(n):
    ss=ssrs[ind]
    roi_number=ss[(0x3006,0x0022)].value
    roi_name=re.sub(r'\W+', '', ss[(0x3006,0x0026)].value)
    print('ROI number {}, name {}'.format(roi_number,roi_name))
    if roi_name in exclude_labels: 
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
    for k in range(nContours):
        contour=contour_sequence[k]
        npts=contour.NumberOfContourPoints
        nptsAcc+=npts
        print('Contour {}, number of points: {}'.format(k+1,npts))
        pts=contour.ContourData
        poly2d,z=pts2poly(pts,origin,mm_vox_size)
        #print(z)
        #print (poly2d)
        #continue
        imslice=current_region_code*get_rasterized_poly_slice(poly2d,imwidth,imheight)
        #print(poly2d)
        #print(imslice)
        rtss_voxels[:,:,z]=imslice
        #break
    
    v=display_color.value
    color='0x{:02X}{:02X}{:02X}'.format(int(v[0]),int(v[1]),int(v[2]))
    roi_descriptor=dict(roi_number=roi_number,
                        roi_name=roi_name,
                        display_color=color,
                        num_contours=nContours,
                        points_in_all_contours=nptsAcc,
                        intensity_value=current_region_code)
                        
    roi_list.append(roi_descriptor)        
    current_region_code+=1    
    
#create and save nifti image.
#flip axes to orient to RAS space

flips=np.array([-1.,-1.,1.,1.])
nifti_affine=np.diag(np.array([xPixelSize,yPixelSize,zPixelSize,1])*flips)
nifti_affine[:3,3]=origin*flips[:-1]

nifti_image_rtss=Nifti1Image(rtss_voxels,nifti_affine)
nifti_image_struct=Nifti1Image(struct_voxels,nifti_affine)

print ('writing',output_rtss_nii)
nibabel.nifti1.save(nifti_image_rtss,output_rtss_nii)

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
    parser.add_argument("input_rtstruct_dicom_file", help="Path to input DICOM RTSTRUCT file")
    parser.add_argument("input_structural_dicom_dir", help="Path to input structural DICOM directory")
    parser.add_argument("output_nifti_rtstruct", help="Path to output NIFTI image")
    parser.add_argument("--out_structural", metavar="<string>",type=str,default=None, 
                        help="Path to output NIFTI image [output_nifti_rtss+struct.nii]")
    parser.add_argument("--exclude_labels", metavar="<string>",type=str,default=None,
                        help="Exclude labels, comma separated [None]")
    

    return parser.parse_args()
    
if __name__ == "__main__":
    p = get_parser()
    print(p)
    out_structural=p.output_nifti_rtstruct+'_struct.nii' if out_structural is None else p.out_structural
    
    exc_labels=[] if p.exclude_labels is None else p.exclude_labels.split(',')
    rtss_to_nifti(input_rtstruct_dicom:str, input_structural_dicom:str,output_rtss_nii:str,
                  output_struct_nii:str,exclude_labels:list)
    
    