'''
Author: Mikhail Milchenko, mmilchenko@wustl.edu
Copyright (c) 2020, Computational Imaging Lab, School of Medicine, Washington University in Saint Louis

Redistribution and use in source and binary forms, for any purpose, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import json,os,os.path,sys,argparse
import numpy as np,skimage,nibabel as nib, nibabel.processing, nibabel.funcs
from skimage import measure, filters, morphology
from skimage.transform import rescale, resize

def get_cube_type(max_size):
    tum_size_map=[
        {'range':[0,5],'cube_dim':15,'label':'small'},
        {'range':[5,10],'cube_dim':20,'label':'medium'},
        {'range':[10,20],'cube_dim':30,'label':'large'},
        {'range':[20,60],'cube_dim':70,'label':'ex_large'},
        {'range':[60,120],'cube_dim':120,'label':'max_size'}
    ]
    for e in tum_size_map:
        rng=e['range']
        if max_size>=rng[0] and max_size<rng[1]: return e
    print('ROI size {}mm exceeded maximum ROI size ({}mm)!'.format(max_size,tum_size_map[-1]['range'][1]))
    return None
    
def resample_image_111(img_in,is_mask):
    order=0 if is_mask else 3
    return nib.processing.resample_to_output(img_in,[1,1,1],order=order)

def conform_image_111(img_in,img_ref,is_mask):
    order=0 if is_mask else 3
    #print(img_in.affine)
    
    #Output orientation for now is RAS.
    out=nibabel.processing.conform(img_in,img_ref.shape,voxel_size=(1.0,1.0,1.0),order=order)
    #out.affine=np.diag(np.sign(np.diagonal(img_in.affine)))
    #out.set_qform(np.diag(np.sign(np.diagonal(img_in.affine))))
    #print(out.affine)
    return out

def split_image(file_im,files_masks):
    try:
        im0=nib.load(file_im)
        if len(im0.shape)<3 or len(im0.shape)>4:
            print('input image shape is <3 or >4!'); return None      
        if len(im0.shape)==4:
            im=resample_image_111(nibabel.funcs.four_to_three(im0)[0],False)
        else:
            im=resample_image_111(im0,False)   
        masks=[]
        T=np.array(im.shape)
        
        for mf in files_masks:
            msk=nib.load(mf)
            if msk.shape[:3] != im0.shape[:3]:
                print('ERROR: Dimensions of input file {} and mask {} do not match!'.format(im.shape,msk.shape))
                return None
            masks+=[conform_image_111(msk,im,True)]
        return im,masks,T
    except:
        print('ERROR: cannot read input file(s)')
        print(sys.exc_info()[1])
        return None
        
def get_subimages(img,mask,T):
    props=skimage.measure.regionprops(mask.get_fdata().astype('int'))
    bb=props[0].bbox
    bb_range=[(bb[i],bb[i+3]) for i in range(0,3) ]
    bb_center=[ bb_range[i][0]+(bb_range[i][1]-bb_range[i][0])/2 for i in range (0,3) ]
    
    #print('bbox',bb_range)
    dims=np.array([bb[3]-bb[0],bb[4]-bb[1],bb[5]-bb[2]]).astype('int')
    #print('bbox size:',dims)
    #print('bbox center:',bb_center)
    
    tp=get_cube_type(dims.max())
    if tp is None:
        print('get_subimages: invalid input object size')
        return None,None
    
    cube_dim=tp['cube_dim']
    #print('cube_dim',cube_dim)
 
    rng=[]
    for i in range (3):
        d0=int(bb_center[i]-(cube_dim/2))
        d1=int(d0+cube_dim)
        #print('axis',i,'bounds',d0,d1)
        rng+=[[d0,d1]]
 
    print('output region:', rng)
    sub,hdr=nifti_subimage_111(img,rng,T,False)
    msk,_=nifti_subimage_111(mask,rng,T,True)
    return sub,msk,hdr

#bounds are a 2d array [[st,en],[st,en],[st,en]]
def nifti_subimage_111(img,bounds,T,is_mask):
    voxels=img.get_fdata()
    ish=np.array(voxels.shape)
    dims=np.array([bounds[i][1]-bounds[i][0] for i in range(3)])
    
    subim=np.zeros(dims)
    ist,ien=[0,0,0],[0,0,0]
    sst,sen=[0,0,0],[0,0,0]
    
    for i in range(3):
        i0,i1=0,ish[i]
        s0,s1=bounds[i][0],bounds[i][1]
        if s1<i0 or s0>i1: print('Intersection of mask and image is empty!'); return None #the sub-range is outside the original image
        ist[i]=s0 if s0>=i0 else 0
        ien[i]=s1 if s1<=i1 else i1
        sst[i]=0 if s0>=i0 else i0-s0
        sen[i]=int(dims[i]) if s1<=i1 else int(dims[i]-(s1-i1))
    
    subim[sst[0]:sen[0],sst[1]:sen[1],sst[2]:sen[2]]=voxels[ist[0]:ien[0],ist[1]:ien[1],ist[2]:ien[2]]
    img=nib.nifti1.Nifti1Image( subim, np.diag(np.sign( np.diagonal(img.affine) ) ) )
    
    #print ("image shape", ish, "subrange", bounds, "sst", sst, "sen", sen,"ist",ist,"ien",ien)
    return img,dict(scaled_src_dims=T.tolist(),targ_range=[ [sst[0],sen[0]],[sst[1],sen[1]],[sst[2],sen[2]] ],scaled_src_range=[[ist[0],ien[0]],[ist[1],ien[1]],[ist[2],ien[2]]])

def get_slice(r):
    '''
    Create slice object form 2-d array representing range in a 3D image.    
    '''
    if isinstance(r[0],list):
        return slice(r[0][0],r[0][1]),slice(r[1][0],r[1][1]),slice(r[2][0],r[2][1])
    else:
        return slice(0,r[0]),slice(0,r[1]),slice(0,r[2])

def subimage2image(sub_image,ref_image,json_header,is_mask=True, overlay=False):
    sub,ref,js=sub_image,ref_image,json_header
    #debug
    #sub.to_filename('/data/temp/sub.nii')
    #ref.to_filename('/data/temp/ref')
    print(js)
    
    #1. Create a resampled copy of the reference image.   
    d=np.zeros(js['scaled_src_dims'])
    #2. Init source and target ranges
    ssrc=get_slice(js['scaled_src_range'])
    starg=get_slice(js['targ_range'])
    
    d[ssrc]=sub.get_fdata()[starg]
    ref_sampled_im=nib.Nifti1Image(d,np.diag(np.sign(np.diagonal(ref_image.affine))))
    #debug
    #ref_sampled_im.to_filename('/data/temp/ref_sampled_im')    
    pd=ref._header['pixdim']
    
    #interpolation order
    order=0 if is_mask else 3
    res=nibabel.processing.conform(ref_sampled_im,ref.shape,voxel_size=(pd[1],pd[2],pd[3]),order=order)
    res.set_sform(ref.affine)
    #print(res.affine)
    return res



class DefParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
        
if __name__=="__main__":
    p=DefParser(description='Convert between mask based ROI and original image',formatter_class=argparse.RawTextHelpFormatter,epilog='')
                
    p.add_argument('command',type=str,help='\nroi2subim\tconvert ROI mask to subimage\n\
\t\treqiured: --roi, optional: --img. \n\
\t\te.g. subimage_convert.py subim2roi --roi my_roi\
        \nsubim2roi\tconvert subimage to ROI mask\n\
\t\trequired args: --sub_img or --sub_roi; --img\n\
        \te.g. subimage_convert.py roi2subim my_image none none my_sub_roi ')
    
    p.add_argument('--img',metavar='<nifti file>', type=str,help='root of the original image')
    p.add_argument('--roi',metavar='<nifti file>', type=str,help='root of ROI mask in the space of the original image')
    p.add_argument('--sub_img',metavar='<nifti file>',type=str,help='subimage root')
    p.add_argument('--sub_roi',metavar='<nifti file>',type=str,help='subimage based ROI root')
    p.add_argument('--suffix',metavar='<string>', type=str,help='output file suffix [_roi2subim, _subim2roi]')
                
    a=p.parse_args()
    suff=a.suffix if a.suffix is not None else '_'+a.command.replace('--','')
                
    if a.command=='roi2subim':  
        if a.roi=='none': print('unsupported input!'); sys.exit(-1)
        if a.img=='none': a.img=a.roi
              
        print('reading',a.img,a.roi)
        im,masks,T=split_image(a.img,[a.roi])
        subim,submask,header_dict=get_subimages(im,masks[0],T)
        
        out_sub_roi=a.roi.replace('.nii','')+suff
        out_sub_img=a.img.replace('.nii','')+suff
        
        print('writing',out_sub_roi+'.nii'); submask.to_filename(out_sub_roi+'.nii')
        jout=out_sub_roi+'.json'; print('writing',jout)
        with open(jout,'w') as f: json.dump(header_dict,f)
                
        print('writing',out_sub_img+'.nii'); subim.to_filename(out_sub_img+'.nii')
        jout=out_sub_img+'.json'; print('writing',jout)
        with open(jout,'w') as f: json.dump(header_dict,f)
            
    elif a.command=='subim2roi':
        img,img_roi,sub_img=None,None,None
        
        if a.img is not None: print('reading',a.img); img=nib.load(a.img)        
        jin=None
        if a.sub_img is not None: 
                print('reading',a.sub_img); sub_img=nib.load(a.sub_img)
                jin_name=a.sub_img.replace('.nii','')+'.json'
                if os.path.isfile(jin_name): 
                    with open(jin_name,'r') as f: jin=json.loads(f.read())
                
        if a.sub_roi is not None:
                print('reading',a.sub_roi); sub_roi=nib.load(a.sub_roi)
                jin_name=a.sub_roi.replace('.nii','')+'.json'
                if os.path.isfile(jin_name):
                    with open(jin_name,'r') as f: jin=json.loads(f.read())
                
        if jin is None: 
                print('a json file with the same root as sub_img or sub_roi is required')
                sys.exit(-1)
        
        if a.sub_img is not None:
            res=subimage2image(sub_img,img,jin,False)
            fout=a.sub_img.replace('.nii','')+suff; print('writing',fout)
            res.to_filename(fout)
        
        if a.sub_roi is not None:
            res=subimage2image(sub_roi,img,jin,True)
            fout=a.sub_roi.replace('.nii','')+suff; print('writing',fout)
            res.to_filename(fout)
    else:
        print('command not understood'); exit (-1)
                
 
