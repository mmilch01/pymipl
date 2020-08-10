# pymipl

## dcmrt_to_subvol

Input: DICOM dirs of a structural MRI and RTSTRUCT <br> 
Output: individual subvolumes of all structures (except for Skull) in NIFTI format<br>

usage: dcmrt_to_subvol \<structural DICOM dir\> \<RTSTRUCT DICOM dir\> \<output_file_name_prefix\> 

Requrements: <br>
64-bit Linux <br>
Python 3.5+ with the following packages:<br>

argparse,numpy,skimage,nibabel

## nifti2rtss.py

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

