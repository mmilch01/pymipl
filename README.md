# pymipl

## dcmrt_to_subvol

Input: DICOM dirs of a structural MRI and RTSTRUCT <br> 
Output: individual subvolumes of all structures (except for Skull) in NIFTI format<br>

usage: dcmrt_to_subvol \<structural DICOM dir\> \<RTSTRUCT DICOM dir\> \<output_prefix\> 

Requrements: <br>
64-bit Linux <br>
Python 3.5+ with the following packages:<br>

argparse,numpy,skimage,nibabel

