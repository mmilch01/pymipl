# pymipl

A collection of tools for working with regions of interest (ROI's) in DICOM RT and NIFTI formats. Includes extraction of subvolumes and format interconversion utilities. All tools require Python 3.5+ to run and are tested on 64-bit Linux systems.

## dcmrt_to_subvol

Convert a DICOM RT structure set to a set of NIFTI images, each containing individual structure set.

Input: DICOM dirs of a structural MRI and RTSTRUCT <br> 
Output: individual subvolumes of all structures (except for Skull) in NIFTI format<br>

usage: dcmrt_to_subvol \<structural DICOM dir\> \<RTSTRUCT DICOM dir\> \<output_file_name_prefix\> 

## convert_subvol

Convert a NIFTI ROI subimage to a reference image space, optionally convert to DICOM RTSTRUCT

Input: subimage, reference image <br>
Output: subimage in reference image space, (optionally) RTSS DICOM file with synthesized contour set.<br>
usage: convert_subvol \<NIFTI subimage\> \<original NIFTI image\> [options] <br>

## nifti2rtss.py

Create RTSTRUCT from a NIFTI binary volume and structural MRI. 

Input: reference DICOM directory, NIFTI mask volume <br> 
Output: RTSTRUCT with conturs created from this NIFTI mask referencing the reference DICOM series.

usage: nifti2rtss.py [-h] [--structure_label <string>] [--tolerance <float>] [--min_poly_pts <int>] input_nifti input_dicom output_dicom<br>

## rtss2nifti.py
Convert DICOM RT structure images to NIFTI
usage: <br>
Input: DICOM RTSTRUCT file, referenced structural DICOM scan
Output: NIFTI files for structural and mask files and metadata in JSON format.

rtss2nifti.py [-h] [--out_struct <string>] [--exclude_labels <string>] [--separate_masks] in_rtss in_struct_dir out_roi_mask

## nifti2mesh.py
Convert a NIFTI binary mask to a mesh file.
usage: python nifti2mesh.py [--min_mask_value <int>] [--no_mesh_smoothing] in_nifti_file out_mesh_file

input: 3D binary mask
output: 3D mesh file. Output formats are those supported by <a href="https://pypi.org/project/meshio">meshio</a>
