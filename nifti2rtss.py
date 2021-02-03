'''
This script was derived from: https://github.com/wanderine/nnunetdocker/convert_to_RTSTRUCT.py v. 07.31.2020

Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2021, Computational Imaging Lab, Washington University School of Medicine

Redistribution and use in source and binary forms, for any purpose, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import random, nibabel as nib, argparse, numpy as np, os
from datetime import datetime

from skimage import measure

import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

def concatenate_coordinates(coordinates_x, coordinates_y, coordinates_z):

    vector = np.zeros((len(coordinates_x)*3,1))

    for i in range(len(coordinates_x)):
        vector[i*3+0] = coordinates_x[i] 
        vector[i*3+1] = coordinates_y[i]
        vector[i*3+2] = coordinates_z[i]

    return vector



def sort_dcms_by_slice_pos(input_dicom_path,dcm_files):
    dcmss=[]
    for idx,dcm in enumerate(dcm_files):
        ds = pydicom.dcmread(os.path.join(input_dicom_path,dcm), stop_before_pixels=True)
        if idx==0:           
           if 'ImagePositionPatient' in ds: sortTag='ImagePositionPatient'
           elif 'SliceLocation' in ds: sortTag='SliceLocation'
           else: return None
        if not sortTag in ds: return None
        if sortTag=='ImagePositionPatient': z=ds.ImagePositionPatient[2]
        else: z=ds.SliceLocation
        dcmss+=[dict(file=dcm,dataset=ds,z=z)]
    return sorted(dcmss, key=lambda dcms: dcms['z'])

def create_rtss_dataset(dicoms_sorted,structure_label):
    rf=dicoms_sorted[0]['dataset']

    SOP_class_UID='1.2.840.10008.5.1.4.1.1.481.3'
    SOP_inst_UID,ser_inst_UID=generate_uid(),generate_uid()
    dt0=datetime.min
    date0,time0=dt0.strftime("%Y%m%d"),dt0.strftime("%H%M%S")
    dt=datetime.now()
    date,time=dt.strftime("%Y%m%d"),dt.strftime("%H%M%S")

    meta=Dataset()
    meta.FileMetaInformationGroupLength = 198
    meta.FileMetaInformationVersion = bytes('01', 'utf-8') # '\x00\x01'
    meta.MediaStorageSOPClassUID = SOP_class_UID
    meta.MediaStorageSOPInstanceUID = SOP_inst_UID
    meta.TransferSyntaxUID = '1.2.840.10008.1.2'
    meta.ImplementationClassUID = '1.2.40.0.13.1.1.1'
    meta.ImplementationVersionName = u'1.0'

    r=Dataset()
    r.Manufacturer,r.StructureSetLabel,r.file_meta=u'NRG',structure_label,meta
    r.OperatorsName=u'nifti2rtss'
    r.is_implicit_VR,r.is_little_endian=True,True
    r.SpecificCharacterSet = 'ISO_IR 100'
    r.InstanceCreationDate = date
    r.InstanceCreationTime = time
    r.SOPClassUID=SOP_class_UID
    r.SOPInstanceUID=SOP_inst_UID
    r.InstanceNumber='1'
    r.SeriesNumber=None

    r.StudyDate=rf.StudyDate if 'StudyDate' in rf else date0
    r.StudyTime=rf.StudyTime if 'StudyTime' in rf else time0
    
    r.AccessionNumber=rf.AccessionNumber if 'AccessionNumber' in rf else None
    r.StudyDescription,r.StudyInstanceUID,r.StudyID=rf.StudyDescription,rf.StudyInstanceUID,rf.StudyID

    r.PatientName,r.PatientID=rf.PatientName,rf.PatientID
    r.PatientBirthDate=''
    r.PatientSex,r.ReferringPhysicianName=rf.PatientSex,rf.ReferringPhysicianName

    r.Modality='RTSTRUCT'
    r.SeriesInstanceUID=ser_inst_UID
    r.SeriesDescription=u'RTSS generated by nifti2rtss'
    r.SeriesDate,r.SeriesTime=date,time
    r.StructureSetDate,r.StructureSetTime=date,time
    
    #1. referenced frame of reference sequence
    referenced_frame_of_ref_seq=Sequence()
    r.ReferencedFrameOfReferenceSequence=referenced_frame_of_ref_seq

    #2. Referenced frame of reference #1
    referenced_frame_of_ref1=Dataset()
    referenced_frame_of_ref1.FrameOfReferenceUID=rf.FrameOfReferenceUID
    
    #3. RT referenced study sequence
    rt_referenced_study_seq=Sequence()
    referenced_frame_of_ref1.RTReferencedStudySequence=rt_referenced_study_seq
    

    #4. RT referenced study sequence, study #1
    rt_referenced_study1=Dataset()
    rt_referenced_study1.ReferencedSOPClassUID=rf.SOPClassUID
    rt_referenced_study1.ReferencedSOPInstanceUID=rf.StudyInstanceUID
    

    #5. RT referenced series sequence
    rt_referenced_series_seq=Sequence()
    rt_referenced_study1.RTReferencedSeriesSequence=rt_referenced_series_seq

    #6. RT referenced series 1
    rt_referenced_series1=Dataset()
    rt_referenced_series1.SeriesInstanceUID=rf.SeriesInstanceUID

    #7. Contour image sequence
    contour_image_sequence=Sequence()
    rt_referenced_series1.ContourImageSequence=contour_image_sequence
 
    #Loop over all DICOM images
    for dcms in dicoms_sorted:  #range(1,numberOfDicomImages+1):
        dstemp = dcms['dataset']
        # Contour Image Sequence: Contour Image
        contour_image = Dataset()
        contour_image.ReferencedSOPClassUID = dstemp.SOPClassUID
        contour_image.ReferencedSOPInstanceUID = dstemp.SOPInstanceUID
        contour_image_sequence.append(contour_image)
   
    #append all sequences
    rt_referenced_series_seq.append(rt_referenced_series1)
    #print('rt_referenced_series_seq', rt_referenced_series_seq)
    rt_referenced_study_seq.append(rt_referenced_study1)
    referenced_frame_of_ref_seq.append(referenced_frame_of_ref1)
    
    #8. Structure set ROI sequence
    structure_set_roi_sequence=Sequence()
    r.StructureSetROISequence=structure_set_roi_sequence

    # Structure set ROI #1
    #structure_set_roi1=Dataset(); ssr1=structure_set_roi1
    #ssr1.ROINumber,ROIName,ROIDescription="1","na","na"
    #ssr1.ROIGenerationAlgorithm='AUTOMATIC'
    #structure_set_roi_sequence.append(ssr1)

    return r
    

def convert(input_nifti_path: str, input_dicom_path: str, output_dicom_path: str, structure_label,poly_approx_tol,min_poly_pts):

    tol=poly_approx_tol
    
    #---------------
    # First DICOM part
    #---------------

    # Get number of DICOM files in DICOM path
    dicomFiles = next(os.walk(input_dicom_path))[2]

    numberOfDicomImages = len(dicomFiles)
    numberOfROIs = 1   # The whole volume is 1 ROI, assuming 1 tumour per patient
    
    # Load template DICOM file header (first file)
    dicomsSorted=sort_dcms_by_slice_pos(input_dicom_path,dicomFiles)

    ds = dicomsSorted[0]['dataset']
    #ds.dir()
    #return
    #print(ds)
    
    xPixelSize = ds.PixelSpacing[0]
    yPixelSize = ds.PixelSpacing[1]

    xyPixelSize=0.5*(xPixelSize+yPixelSize)
    poly_approx_tol/=xyPixelSize

    zPixelSize = ds.SliceThickness
    
    print("Each voxel is ",xPixelSize," x ",yPixelSize," x ",zPixelSize,'tolerance:',poly_approx_tol,'voxels')

    # Find position of first slice
    patientPosition = ds.ImagePositionPatient
    
    patientStartingZ = dicomsSorted[0]['z']
    
    print('Patient position is ', patientPosition[:2])
    print('First slice at ', patientStartingZ)

    #---------------
    # NIFTI part
    #---------------

    # Load nifti volume
    
    #nii = nib.load(input_nifti_path)
    
    nii0 = nib.load(input_nifti_path)
    flips=np.sign(nii0.affine)
    nii=nii0.as_reoriented([[0,-1*flips[0,0]],[1,-1*flips[1,1]],[2,flips[2,2]]])
    print("axes flips:", [[0,-1*flips[0,0]],[1,-1*flips[1,1]],[2,flips[2,2]]])
    
    volume = nii.get_fdata()
    volume = volume.astype(float)

    print('Nifti file dimensions:',volume.shape)

    AllCoordinates = []

    if len(volume.shape)==4: 
        volume = volume[...,0]
        print('Assuming the first channel of the input nifti is the seg mask.')
    elif len(volume.shape)==3:
        print('Segmentation mask has the same number of dimensions as the input volume.')
    else:
        print('Dimension not supported.')
            
    # Loop over slices in volume, get contours for each slice
    #DEBUG: inited=0
    for slice in range(volume.shape[2]):
        AllCoordinatesThisSlice = []
	
        image = volume[:,:,slice]
        
        # Get contours in this slice using scikit-image
        contours = measure.find_contours(image, 0.5)
	    
	    
        #if len(contours)>0: print("slice ",slice,", contours:", len(contours))

        # Save contours for later use
        for n, contour in enumerate(contours):
        
            #DEBUG: if inited==2: continue
            #print("n is ",n,"for slice ",slice)
            nc=len(contour[:,0])
            cont1=measure.approximate_polygon(contour,2)
            nCoordinates = len(cont1[:,0])
            #if nc>1000: print('Large contour before approximation: ',nc,'pts, after:',nCoordinates)
            if nCoordinates<min_poly_pts: continue
            #DEBUG: inited+=1
            #print("number of coordinates is ",len(contour[:,0])*3," for contour ",n," for slice ",slice)
            zcoordinates = slice * np.ones((nCoordinates,1)) 
            
            # Add patient position offset
            reg_contour = np.append(cont1, zcoordinates, -1)
            # Assume no other orientations for simplicity
            reg_contour[:,0] = reg_contour[:,0] * xPixelSize + patientPosition[0]
            reg_contour[:,1] = reg_contour[:,1] * yPixelSize + patientPosition[1]
            #z coordinate will be fixed later.
            reg_contour[:,2] = reg_contour[:,2] * zPixelSize + patientStartingZ

            # Storing coordinates as mm instead of as voxels
            #coordinates = concatenate_coordinates(contour[:,0] * xPixelSize, contour[:,1] * yPixelSize, zcoordinates * zPixelSize)
            coordinates = concatenate_coordinates(*reg_contour.T)
            coordinates = np.squeeze(coordinates)
            
            AllCoordinatesThisSlice.append(coordinates)

        AllCoordinates.append(AllCoordinatesThisSlice)

    #---------------
    # Second DICOM part (RTstruct)
    #---------------
    rtds=create_rtss_dataset(dicomsSorted,structure_label)

    # Structure Set ROI Sequence
    structure_set_roi_sequence = rtds.StructureSetROISequence
    rtds.StructureSetLabel = structure_label

    print('Number of ROIs:',numberOfROIs)
    # Loop over ROIs
    for ROI in range(1,numberOfROIs+1):
        # Structure Set ROI Sequence: Structure Set ROI
        structure_set_roi = Dataset()
        structure_set_roi.ROINumber = str(ROI)
        structure_set_roi.ReferencedFrameOfReferenceUID = ds.FrameOfReferenceUID 
        structure_set_roi.ROIName = 'ROI_' + str(ROI)
        structure_set_roi.ROIGenerationAlgorithm = 'AUTOMATIC'
        structure_set_roi_sequence.append(structure_set_roi)
	
    # ROI Contour Sequence
    roi_contour_sequence = Sequence()
    rtds.ROIContourSequence = roi_contour_sequence

    # Loop over ROI contour sequences
    for ROI in range(1,numberOfROIs+1):

        # ROI Contour Sequence: ROI Contour 1
        roi_contour = Dataset()
        roi_contour.ROIDisplayColor = [0, 230, 0]

        # Contour Sequence
        contour_sequence = Sequence()
        roi_contour.ContourSequence = contour_sequence
        cnumber=0

        # Loop over slices in volume (ROI)
        for slice in range(volume.shape[2]):

            # Should Contour Sequence be inside this loop?
            #roi_contour.ContourSequence = contour_sequence

            # Loop over contour sequences in this slice
            numberOfContoursInThisSlice = len(AllCoordinates[slice])
            if numberOfContoursInThisSlice < 1: continue

            # Contour Image Sequence
            contour_image_sequence = Sequence()
            contour_image1 = Dataset()
            contour_image1.ReferencedSOPClassUID = ds.SOPClassUID
            contour_image1.ReferencedSOPInstanceUID=dicomsSorted[slice]['dataset'].SOPInstanceUID
            contour_image_sequence.append(contour_image1) #one image per contour    

            for c in range(numberOfContoursInThisSlice):

                currentCoordinates = AllCoordinates[slice][c]
                
                # Contour Sequence: Contour 1
                contour = Dataset()
                contour.ContourImageSequence = contour_image_sequence
                
                currentCoordinates[2::3]=dicomsSorted[slice]['z']
                cnumber+=1
                contour.ContourGeometricType = 'CLOSED_PLANAR'                
                contour.NumberOfContourPoints = len(currentCoordinates)/3
                contour.ContourData = currentCoordinates.tolist()
                contour_sequence.append(contour)
                #DEBUG print('contourData',currentCoordinates.tolist())

        roi_contour.ReferencedROINumber = ROI
        roi_contour_sequence.append(roi_contour)


    # RT ROI Observations Sequence
    rtroi_observations_sequence = Sequence()
    rtds.RTROIObservationsSequence = rtroi_observations_sequence

    # Loop over ROI observations
    for ROI in range(1,numberOfROIs+1):
        # RT ROI Observations Sequence: RT ROI Observations 1
        rtroi_observations = Dataset()
        rtroi_observations.ObservationNumber = str(ROI)
        rtroi_observations.ReferencedROINumber = str(ROI)
        rtroi_observations.RTROIInterpretedType = 'ORGAN'
        
        rtroi_observations.ROIInterpreter = ''
        rtroi_observations_sequence.append(rtroi_observations)

    rtds.ApprovalStatus='UNAPPROVED'
    RTDCM_name = output_dicom_path
    #ds.is_implicit_VR,ds.is_little_endian=True,True
    print('is_implicit_VR=',ds.is_implicit_VR,'is_little_endian=',ds.is_little_endian)
    pydicom.filewriter.dcmwrite(RTDCM_name,rtds,write_like_original=False)
    #rtds.save_as(RTDCM_name)
    print('RTSTRUCT saved as %s'%RTDCM_name)
    
def get_parser():
    """
    Parse input arguments.
    """
    parser = argparse.ArgumentParser(description='Convert nifti images to RTSTRUCT file')

    # Positional arguments.
    parser.add_argument("input_nifti", help="Path to input NIFTI image")
    parser.add_argument("input_dicom", help="Path to input DICOM images")
    parser.add_argument("output_dicom", help="Path to output DICOM image")
    parser.add_argument("--structure_label",metavar="<string>",type=str,default="ROI1",help='structure set label')
    parser.add_argument("--tolerance",metavar="<float>", type=float, default=1,help="polygon approximation tolerance (mm)")
    parser.add_argument("--min_poly_pts", metavar="<int>",type=int,default=3,help="minimum number of points in polygon")

    return parser.parse_args()

if __name__ == "__main__":
    p = get_parser()
    print(p)
    convert(p.input_nifti, p.input_dicom, p.output_dicom, p.structure_label,p.tolerance,p.min_poly_pts)
