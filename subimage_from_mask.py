import os,sys,argparse,numpy as np,skimage,nibabel as nib, nibabel.processing, nibabel.funcs
from skimage import measure, filters, morphology
from skimage.transform import rescale, resize

def get_cube_type(max_size):
    tum_size_map=[
        {'range':[0,5],'cube_dim':15,'label':'small'},
        {'range':[5,10],'cube_dim':20,'label':'medium'},
        {'range':[10,20],'cube_dim':30,'label':'large'},
        {'range':[20,60],'cube_dim':70,'label':'ex_large'}
    ]
    for e in tum_size_map:
        rng=e['range']
        if max_size>=rng[0] and max_size<rng[1]: return e
    return None
    
def resample_image_111(img_in,is_mask):
    order=0 if is_mask else 3
    return nib.processing.resample_to_output(img_in,[1,1,1],order=order)

def conform_image_111(img_in,img_ref,is_mask):
    order=0 if is_mask else 3
    return nibabel.processing.conform(img_in,img_ref.shape,voxel_size=(1.0,1.0,1.0),order=order)

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
        for mf in files_masks:
            msk=nib.load(mf)
            if msk.shape[:3] != im0.shape[:3]:
                print('ERROR: Dimensions of input file {} and mask {} do not match!'.format(im.shape,msk.shape))
                return None
            masks+=[conform_image_111(msk,im,True)]
        return im,masks
    except:
        print('ERROR: cannot read input file(s)')
        print(sys.exc_info()[1])
        return None
       
def get_subimages(img,mask):
    props=skimage.measure.regionprops(mask.get_fdata().astype('int'))
    bb=props[0].bbox
    bb_range=[(bb[i],bb[i+3]) for i in range(0,3) ]
    bb_center=[ bb_range[i][0]+(bb_range[i][1]-bb_range[i][0])/2 for i in range (0,3) ]
    
    #print('bbox',bb_range)
    dims=np.array([bb[3]-bb[0],bb[4]-bb[1],bb[5]-bb[2]]).astype('int')
    #print('bbox size:',dims)
    #print('bbox center:',bb_center)
    
    tp=get_cube_type(dims.max())
    if tp is None: return None,None
    
    cube_dim=tp['cube_dim']
    #print('cube_dim',cube_dim)
 
    rng=[]
    for i in range (3):
        d0=int(bb_center[i]-(cube_dim/2))
        d1=int(d0+cube_dim)
        #print('axis',i,'bounds',d0,d1)
        rng+=[[d0,d1]]
 
    print('output region:', rng)
    return nifti_subimage_111(img,rng,False),nifti_subimage_111(mask,rng,True)

#bounds are a 2d array [[st,en],[st,en],[st,en]]
def nifti_subimage_111(img,bounds,is_mask):
    voxels=img.get_fdata()
    ish=np.array(voxels.shape)
    dims=np.array([bounds[i][1]-bounds[i][0] for i in range(3)])
    
    subim=np.zeros(dims)
    ist,ien=[0,0,0],[0,0,0]
    sst,sen=[0,0,0],[0,0,0]
    
    for i in range(3):
        i0,i1=0,ish[i]
        s0,s1=bounds[i][0],bounds[i][1]
        if s1<i0 or s0>i1: print('exiting'); return None #the sub-range is outside original image
        ist[i]=s0 if s0>=i0 else 0
        ien[i]=s1 if s1<=i1 else i1
        sst[i]=0 if s0>=i0 else i0-s0
        sen[i]=dims[i] if s1<=i1 else dims[i]-(s1-i1)
    
    subim[sst[0]:sen[0],sst[1]:sen[1],sst[2]:sen[2]]=voxels[ist[0]:ien[0],ist[1]:ien[1],ist[2]:ien[2]]
    img=nib.nifti1.Nifti1Image(subim,np.eye(4))
    #print ("image shape", ish, "subrange", bounds, "sst", sst, "sen", sen,"ist",ist,"ien",ien)
    return img

class DefParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
        
if __name__=="__main__":
    p=DefParser(description='Extract subimage based on mask')
    p.add_argument('image',type=str,help='input NIFTI image')
    p.add_argument('mask',type=str,help='input NIFTI mask')    
    p.add_argument('out_image',type=str,help='output subimage file')
    p.add_argument('out_mask',type=str,help='output submask file')
    a=p.parse_args()
    print('reading',a.image,a.mask)
    im,masks=split_image(a.image,[a.mask])
    if im is None or masks is None: sys.exit (-1)
    subim,submask=get_subimages(im,masks[0])
    if subim is None or submask is None: sys.exit(-1)
    print('writing',a.out_image)
    subim.to_filename(a.out_image)
    print('writing',a.out_mask)
    submask.to_filename(a.out_mask)
    
