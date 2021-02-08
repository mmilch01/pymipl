'''
Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2021, Computational Imaging Lab, School of Medicine, Washington University in Saint Louis

Redistribution and use in source and binary forms, for any purpose, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import re, nibabel as nib, argparse, numpy as np, os, json
from skimage import measure
from PIL import Image,ImageDraw
from nibabel.nifti1 import Nifti1Image,Nifti1Header
import meshio, vtk
from vtk.numpy_interface import dataset_adapter as dsa
from utils import write_rec_file

def get_triangular_mesh_from_vtkPolyData(polydata):
    '''
    Form a triangular mesh understood by meshio, converted from vtkPolyData object.
    Returns arrays of points and polygon indices.
    '''
    
    polys=dsa.WrapDataObject(polydata)
    polys1=np.array(polys.Polygons)
    indices=np.array(range(polys1.shape[0]))
    
    #check that all polygons are actually triangles.
    poly_codes=polys1[ indices[indices % 4 == 0] ]

    if poly_codes.shape[0] != poly_codes[poly_codes == 3].shape[0]:
        raise(ValueError('polydata contains non-triangular shapes!'))    
    
    #extract triangle point indices
    triangles=polys1[ indices[indices % 4 !=0] ]
        
    return np.array(polys.Points),np.reshape(triangles,[triangles.shape[0]//3,3])

def read_NIFTI_into_vtkImage(file:str):
    '''
    Read a NIFTI file into a vtkImage object
    '''
    reader=vtk.vtkNIFTIImageReader()
    reader.SetFileName(file)
    reader.TimeAsVectorOn()
    reader.Update()
    return reader.GetOutput()
                                            
def extract_mesh_from_vtkImage(vtkImage, isovalue:float):  
    '''
    run marching cubes to extract triangulation from image. 
    '''
    vtkCF=vtk.vtkContourFilter()
    vtkCF.SetInputData(vtkImage)
    vtkCF.SetValue(0,isovalue)
    vtkCF.Update()
    return vtkCF.GetOutput()

def smoothen_mesh_vtkPolys(vtkPolys):
    '''
    apply mesh smoothing to a vtk mesh
    '''
    vtkSF=vtk.vtkWindowedSincPolyDataFilter()
    vtkSF.SetNumberOfIterations(10)
    vtkSF.SetInputData(vtkPolys)
    vtkSF.Update()
    return vtkSF.GetOutput() 

def vtk_write_stl(vtkPolys, filename:str):
    '''
    test to write vtkPolys into stl file
    '''
    stlWriter=vtk.vtkSTLWriter()
    stlWriter.SetInputData(vtkPolys)
    stlWriter.SetFileTypeToBinary()
    stlWriter.SetFileName(filename)
    stlWriter.Write()
    
def get_parser():
    """
    Parse input arguments.
    """
    parser = argparse.ArgumentParser(description='Convert NIFTI binary mask to a mesh file.')

    # Positional arguments.
    parser.add_argument("in_nii", help="Input NIFTI file")
    parser.add_argument("out_file", help="Output mesh file. Supported formats listed here: https://pypi.org/project/meshio")
    parser.add_argument("--mask_value", metavar="<int>",type=int,default=1,
                        help="input mask isovalue to create surface [1]")
    return parser.parse_args()
    
if __name__ == "__main__":
    p = get_parser()
    
    print ('reading',in_nii)
    vtkImage=read_NIFTI_into_vtkImage(input_nii)
    vtkPolys=extract_mesh_from_vtkImage(vtkImage,1)
    vtkPolysSmoothed=smoothen_mesh_vtkPolys(vtkPolys)
    points, nodes=get_triangular_mesh_from_vtkPolyData(vtkPolysSmoothed)
    cells=[("triangle",nodes)]
    mesh=meshio.mesh(points,cells)
    print('writing',p.out_file)
    meshio.write(p.out_file,mesh)   
    write_rec_file(p.out_file,infiles=[p.in_rtss,p.in_struct_dir])    
    print('done')
                
    