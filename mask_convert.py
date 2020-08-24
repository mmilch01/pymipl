'''
Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2020, Computational Imaging Lab, School of Medicine, Washington University in Saint Louis

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''


import sys, nibabel as nib, numpy as np, argparse, sys

def split_masks(files, targs, outfile):
    try:
        img1=nib.load(files[0])
        img2=nib.load(files[1])
        t1,t2=int(targs[0]),int(targs[1])        
    except:
        print('ERROR: cannot read input file(s)')
        return False
    if img1.shape != img2.shape: 
        print('Dimensions of input images don\'t match.')
        return False
    msk1,msk2=np.where(img1.get_fdata()>0,1,0)*t1,np.where(img2.get_fdata()>0,1,0)*t2
    msk3=np.multiply(msk1,np.where(msk2==0,1,0))+msk2
    
    img3=nib.Nifti1Image(msk3,img1.affine,header=img1.header)
    try:
        nib.save(img3,outfile)
    except:
        print('ERROR: cannot write output file')
        print(sys.last_traceback)
        return False
    return True    
    
class DefParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
        
if __name__=="__main__":
    p=DefParser(description='Convert two masks to a multi-mask')
    p.add_argument('mask1',type=str,help='first mask')
    p.add_argument('mask1_target',type=str,help='target value for mask1')    
    p.add_argument('mask2',type=str,help='second mask')
    p.add_argument('mask2_target',type=str,help='target value for mask2')
    p.add_argument('outfile',type=str,help='output mask')    
    a=p.parse_args()
    print ('split_masks ',a.mask1,a.mask2,a.mask1_target,a.mask2_target,a.outfile)
    sys.exit(split_masks([a.mask1,a.mask2],[a.mask1_target,a.mask2_target],a.outfile))
    
