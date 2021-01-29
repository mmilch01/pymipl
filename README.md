# pymipl

A collection of tools for working with regions of interest (ROI's) in DICOM RT and NIFTI formats. Includes extraction of subvolumes and format interconversion utilities.

## dcmrt_to_subvol

Convert a DICOM RT structure set to a set of NIFTI images, each containing individual structure set.

Input: DICOM dirs of a structural MRI and RTSTRUCT <br> 
Output: individual subvolumes of all structures (except for Skull) in NIFTI format<br>

usage: dcmrt_to_subvol \<structural DICOM dir\> \<RTSTRUCT DICOM dir\> \<output_file_name_prefix\> 

Requrements: <br>
64-bit Linux <br>
Python 3.5+ with the following packages:<br>

argparse,numpy,skimage,nibabel

## convert_subvol

Convert a NIFTI ROI subimage to a reference image space, optionally convert to DICOM RTSTRUCT

Input: subimage, reference image <br>
Output: subimage in reference image space, (optionally) RTSS DICOM file with synthesized contour set.<br>
usage: convert_subvol \<NIFTI subimage\> \<original NIFTI image\> [options] <br>

options:<br>
    -m                      subimage is a binary ROI mask (required to convert to RTSS)<br>
    -dcm \<dcm_dir\>        referenced structural DICOM series dir (-rtss required)<br>
    -rtss \<file\>          output RT Structure Set file name (-dcm required)<br>
    --min_poly_pts \<int\>  minimum number of points in polygon<br>
    --tolerance \<float\>   polygon approximation tolerance (mm)<br>

## nifti2rtss.py

Create RTSTRUCT from a NIFTI binary volume and structural MRI. 

Input: reference DICOM directory, NIFTI mask volume <br> 
Output: RTSTRUCT with conturs created from this NIFTI mask referencing the reference DICOM series.

usage: nifti2rtss.py [-h] [--structure_label <string>] [--tolerance <float>] [--min_poly_pts <int>] input_nifti input_dicom output_dicom<br>
parameters:<br>
--tolerance:     maximum error (mm) in approximating the original contour [1]
--min_poly_pts:  minimum number of points in included contour [3]

Requrements: <br>
64-bit Linux <br>
Python 3.5+ with the following packages:<br>

nibabel, pydicom
