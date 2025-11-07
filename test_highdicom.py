"""
Converts a DICOM SEG file to a NIfTI file using a reference nifti file 
(the image to which the segmentation applies)
Requirements
------------
- highdicom
- nibabel
- SimpleITK
"""

import os
from typing import Generator

import numpy as np
import nibabel as nib
import SimpleITK as sitk
import highdicom
import argparse
import re


def format_output_path(output_dir: str, nii_name: str, seg_name: str) -> str:
    """
    Formats the output path for a NIfTI file.
    Args:
        output_dir (str): Output directory.
        nii_name (str): Name of the source NIfTI file.
        seg_name (str): Name of the segmentation.
    Returns:
        str: Formatted output path.
    """
    seg_str=re.sub(r'[^a-zA-Z0-9-.]', '_', seg_name)
    output_path = os.path.join(output_dir, f"{nii_name}_{seg_str}.nii.gz")
    return output_path


def read_nii(nii_path: str) -> sitk.Image:
    """
    Reads a NIfTI file and returns its data and affine matrix.
    Args:
        nii_path (str): Path to the NIfTI file.
    Returns:
        tuple[np.ndarray, np.ndarray]: NIfTI data array and affine transformation matrix.
    """
    reader = sitk.ImageFileReader()
    reader.SetFileName(nii_path)
    return reader.Execute()


def read_seg(dicom_seg_path: str) -> highdicom.seg.Segmentation:
    """
    Reads a DICOM SEG file and returns its data and affine matrix.
    Args:
        dicom_seg_path (str): Path to the DICOM SEG file.
    Returns:
        highdicom.seg.Segmentation: DICOM SEG object.
    """
    dicom_seg = highdicom.seg.segread(dicom_seg_path)
    return dicom_seg


def get_affine_from_sitk(image) -> np.ndarray:
    """Extract affine matrix from SimpleITK image."""

    # Get components
    direction = image.GetDirection()
    origin = image.GetOrigin()
    spacing = image.GetSpacing()

    # Reshape direction to 3x3 matrix
    direction_matrix = np.array(direction).reshape(3, 3)

    # Create affine matrix
    affine = np.eye(4)

    # Fill rotation and spacing
    for i in range(3):
        for j in range(3):
            affine[i, j] = direction_matrix[i, j] * spacing[j]

    # Fill translation
    affine[:3, 3] = origin

    return affine


def flip_based_on_affine(seg_data: sitk.Image, ref_image: sitk.Image) -> sitk.Image:
    """
    Checks if the data array needs to be flipped based on the affine matrix.
    """
    src_x, src_y, _ = nib.aff2axcodes(  # pylint: disable=unbalanced-tuple-unpacking
        get_affine_from_sitk(seg_data)
    )
    dest_x, dest_y, _ = nib.aff2axcodes(  # pylint: disable=unbalanced-tuple-unpacking
        get_affine_from_sitk(ref_image)
    )
    flips = [src_x != dest_x, src_y != dest_y, False]
    return sitk.Flip(seg_data, flips)


def format_nifti(seg_data: sitk.Image, ref_nii: sitk.Image) -> sitk.Image:
    """
    Formats a DICOM SEG image into a NIFTI image.
    Args:
        seg_data (sitk.Image): DICOM SEG image.
        ref_nii (sitk.Image): Reference NIFTI image.
    Returns:
        sitk.Image: Formatted NIFTI image.
    """
    # flip the data array based on the affine matrix
    seg_data = flip_based_on_affine(seg_data, ref_nii)
    # account for how NIFTI uses opposite z indexing
    seg_data = sitk.Flip(seg_data, [False, False, True])
    return seg_data


def write_nifti(nii: sitk.Image, nii_path: str):
    """
    Writes a NIfTI object to a file.
    Args:
        nii (sitk.Image): NIfTI object.
        nii_path (str): Path to the output NIfTI file.
    """
    writer = sitk.ImageFileWriter()
    writer.SetFileName(nii_path)
    writer.Execute(nii)


def split_seg_channels(
    seg_data: highdicom.seg.Segmentation,
) -> Generator[np.ndarray, None, None]:
    """
    Splits a DICOM SEG data array into separate channels for each segment.
    Args:
        seg_data (np.ndarray): DICOM SEG data array.
    Returns:
        np.ndarray: Array of separate channels for each segment.
    """
    print(f"Number of segments: {seg_data.number_of_segments}")
    uids = [x[2] for x in seg_data.get_source_image_uids()]
    for i in range(seg_data.number_of_segments):
        yield seg_data.get_pixels_by_source_instance(
            uids,
            segment_numbers=[i + 1],
            ignore_spatial_locations=True,
            assert_missing_frames_are_empty=True,
        )[..., 0]


def copy_sitk_image_info(src: sitk.Image, dst: sitk.Image) -> sitk.Image:
    """
    Copy image information from one SimpleITK image to another.
    Args:
        src (sitk.Image): Source image.
        dst (sitk.Image): Destination image.
    Returns:
        sitk.Image: Destination image with copied information.
    """
    # won't always have the same dims so can't use CopyInformation
    dst.SetSpacing(src.GetSpacing())
    dst.SetOrigin(src.GetOrigin())
    dst.SetDirection(src.GetDirection())
    for key in src.GetMetaDataKeys():
        dst.SetMetaData(key, src.GetMetaData(key))
    return dst


def dicomseg2nii(dicom_seg_path, nii_path, output_path):
    dicom_seg = read_seg(dicom_seg_path)  # (z, x, y)
    sitk_dcm_seg = sitk.ReadImage(dicom_seg_path)
    nii = read_nii(nii_path)  # (x, y, z)
    out_paths = []
    for i, seg_channel in enumerate(split_seg_channels(dicom_seg)):
        nii_out = sitk.GetImageFromArray(seg_channel)
        if nii_out.GetSize() != nii.GetSize():
            raise ValueError(
                f"Segmentation size {nii_out.GetSize()} does not match NIfTI size {nii.GetSize()}"
            )
        #nii_out = copy_sitk_image_info(sitk_dcm_seg, nii_out)
        nii_out=copy_sitk_image_info(nii, nii_out)
        nii_out = format_nifti(nii_out, nii)
        nii_code = (
            dicom_seg.SegmentSequence[i]
            .SegmentedPropertyTypeCodeSequence[0]
            .CodeMeaning
        )
        nii_basename = os.path.basename(nii_path).split(".")[0]
        channel_output_path = format_output_path(
            output_path,
            nii_basename,
            nii_code,
        )
        write_nifti(nii_out, channel_output_path)
        out_paths.append(channel_output_path)

    return out_paths

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="RT Struct converter")

    parser.add_argument("nifti_ref", type=str, help="Path to the DICOM series directory")
    parser.add_argument("dcm_seg_dir", type=str, help="Path to the RT Struct DICOM file")
    parser.add_argument("out_dir", type=str, help="Label prefix for output files")
    args = parser.parse_args()
    dicomseg2nii(args.dcm_seg_dir,args.nifti_ref,args.out_dir)
    
