'''
Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2021, Computational Imaging Lab, School of Medicine, Washington University in Saint Louis

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import nibabel as nib, numpy as np, sys, argparse
from skimage.morphology import label
from skimage.measure import regionprops
import nibabel.nifti1
from utils import write_rec_file

def label_max_prob(mask_weighted,label,nlabel): 
    weights=[]
    for val in range(1,nlabel+1):
        weights.append(np.sum((label==val)*mask_weighted))
    return np.argmax(np.array(weights))+1

class DefParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
        
if __name__=="__main__":
    
    p=DefParser(description='Extract single connected component with highest probability mask')
    p.add_argument('input_binary_mask',type=str,help='input binary mask with multiple segments')
    p.add_argument('input_weights_image',type=str,help='probability map for the given mask')
    p.add_argument('output_binary_mask',type=str,help='output mask with most probable connected component')

    a=p.parse_args()
    mask_file,atlas,out=a.input_binary_mask,a.input_weights_image,a.output_binary_mask

    mask=nib.load(mask_file)
    atl=nib.load(atlas)
    mask_raw=mask.get_fdata()
    atl_raw=atl.get_fdata()

    lb,nlabel=label(np.squeeze(mask_raw),return_num=True,connectivity=1)
    
    mask_weighted=np.squeeze(mask_raw*atl_raw)
    lmax=label_max_prob(mask_weighted,lb,nlabel)
    out_label_image=nibabel.nifti1.Nifti1Image(lb==lmax,None,header=mask.header)

    print ('saving ',out)
    nibabel.nifti1.save(out_label_image,out)
    write_rec_file(out,'nii',[mask_file,atlas])
                