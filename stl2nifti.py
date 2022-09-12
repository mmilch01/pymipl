'''
Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2022, Computational Imaging Lab, School of Medicine, Washington University in Saint Louis

Redistribution and use in source and binary forms, for any purpose, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import argparse,stltovoxel,numpy as np
from nibabel.nifti1 import Nifti1Image,Nifti1Header
import nibabel.nifti1
from utils import write_rec_file
from stl import mesh

def stl2nifti(infile:str,outfile:str, resolution:int, padding_fraction:float):
    #read stl mesh from file
    print("Creating mesh")
    mesh_obj=mesh.Mesh.from_file(infile)    
    
    #make 'organized mesh', format to input to stltovoxel.
    org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
    meshdim=(np.amax(org_mesh,0)-np.amin(org_mesh,0))[0]  

    print('meshdim',meshdim)
    #convert mesh to raster
    print("Converting mesh to raster")
    vol,_,_=stltovoxel.convert_mesh(org_mesh,resolution=resolution)

    #create output nifti file
    nifti_affine=np.diag(np.array([meshdim[2]/vol.shape[0],meshdim[1]/vol.shape[1],meshdim[0]/vol.shape[2],1]))
    
    #pad raster with zeroes, pad size calculated as fraction of each linear dimension
    sh=np.array(vol.shape);
    pad=(sh*padding_fraction).astype(int)        
    vol1=np.pad(vol, [(pad[0],pad[0]),(pad[1],pad[1]),(pad[2],pad[2])] )
    
    #write out
    print("Writing",outfile)
    nifti_image=Nifti1Image(vol1,nifti_affine)
    nibabel.nifti1.save(nifti_image,outfile)

def get_parser():
    """
    Parse input arguments.
    """
    parser = argparse.ArgumentParser(description='Convert stereolithography file to NIFTI binary mask')

    # Positional arguments.
    parser.add_argument("in_stl", help="Input stl file")
    parser.add_argument("out_nii", help="Output NIFTI file")
    parser.add_argument("--resolution", metavar="<int>",type=int,default=512,
                        help="size in voxels of the biggest dimension of the output raster [512]")
    parser.add_argument("--padding_fraction", metavar="<float>",type=float,default=0.25,
                        help="fraction of the output raster size for bilateral zero padding [0.25]")
    
    return parser.parse_args()

if __name__ == "__main__":
    p = get_parser()
    stl2nifti(p.in_stl,p.out_nii,p.resolution,p.padding_fraction)
    write_rec_file(p.out_nii,main_extension='nii',infiles=[p.in_stl])