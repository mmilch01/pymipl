from rt_utils import RTStructBuilder
import matplotlib.pyplot as plt
import argparse
import nibabel as nib
import numpy as np
import SimpleITK as sitk
import re

def buildMaskArray(seriesPath, labelPath):
    """
    Helper for the following function: taken from rt_utils
    """
    rtstruct = RTStructBuilder.create_from(
        dicom_series_path=seriesPath, rt_struct_path=labelPath)
    
    rois = rtstruct.get_roi_names()
    print(f'ROIS defined (all summed up): {rois}')
    good_rois = []
    
    for roi in rois:
#        if not(bad_word_is_in_roi(roi)):
         good_rois.append(roi)
    
    masks = []

    for roi in good_rois:
        mask_3d = np.moveaxis(rtstruct.get_roi_mask_by_name(roi).astype(int),[0, 1, 2], [1, 2, 0])
        masks.append(mask_3d)

    return good_rois, masks

    #final_mask = sum(masks)  # sums element-wise
    #final_mask = np.where(final_mask>=1, 1, 0)
    # Reorient the mask to line up with the reference image
    #final_mask = np.moveaxis(final_mask, [0, 1, 2], [1, 2, 0])

    return final_mask

def buildMasks(structPath, rtPath, outroot):
    """
    To convert the gt to the correct .nii.gz file
    """
    roi_names,masks = buildMaskArray(structPath, rtPath)
    # SUV_* BP Liver
    # Load original DICOM image for reference
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(structPath)
    reader.SetFileNames(dicom_names)
    ref_img = reader.Execute()
    sitk.WriteImage(ref_img,outroot+'_struct.nii')

    for roi_name,mask in zip(roi_names,masks):
        mask_img=sitk.GetImageFromArray(mask)
        mask_img.CopyInformation(ref_img)
        sitk.WriteImage(mask_img,outroot+'_roi_'+re.sub(r'[^a-zA-Z0-9-.]', '_', roi_name)+'.nii',imageIO="NiftiImageIO")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RT Struct converter")

    parser.add_argument("dicom_series_path", type=str, help="Path to the DICOM series directory")
    parser.add_argument("rt_struct_path", type=str, help="Path to the RT Struct DICOM file")
    parser.add_argument("out_label", type=str, help="Label prefix for output files")
    args = parser.parse_args()
    buildMasks(args.dicom_series_path,args.rt_struct_path,'./'+args.out_label)

    #dicom_series_path = "/data/ADAPT/RIDER-1129164940/09-20-2006-1-NA-96508/scans/4-unknown/resources/DICOM"
    #rt_struct_path="/data/ADAPT/RIDER-1129164940/09-20-2006-1-NA-96508/scans/9-TEST/resources/secondary/1-1.dcm"
    #outfile="/data/ADAPT/RIDER-1129164940/09-20-2006-1-NA-96508/outfile.nii"
    #buildMaskArray(dicom_series_path,rt_struct_path)
    #buildMasks(dicom_series_path,rt_struct_path,outfile)
'''
    rtstruct = RTStructBuilder.create_from(
        dicom_series_path=dicom_series_path,
        rt_struct_path=rt_struct_path
    )
    print("Available ROI Names:")
    for roi_name in rtstruct.get_roi_names():
        print(f"saving {roi_name}.nii")
        mask_3d = rtstruct.get_roi_mask_by_name(roi_name)
        nifti_img = nib.Nifti1Image(mask_3d.astype('uint8'), affine=np.eye(4))
        nib.save(nifti_img, f"{roi_name}.nii")
'''
