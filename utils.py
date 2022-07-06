import argparse, getpass, cmdline_provenance as cmdprov, os.path, nibabel as nib

class IFH:
    def __init__(self):
        '''
        interfile=""
        version_of_keys="3.3"
        conversion_program="tools_4dfp.py" #program that generated the file
        name_of_data_file=""               #image data filename"
        number_format="float"
        imagedata_byte_order="littleendian"#bigendian
        number_of_bytes_per_pixel=4
        number_of_dimensions=4
        matrix_size=np.array([0,0,0,0])    #image array dims
        orientation=2                      #2: transverse; 3: coronal; 4: sagittal
        scaling_factor=np.array([1.0,1.0,1.0]) #mm/pixel
        mmpix=np.array([1.0,1.0,1.0])      #signed scaling factor
        center=np.array([0.0,0.0,0.0])     #image center
        '''
        self._vals=dict()
        self._int_scalars=['number of bytes per pixel','number of dimensions',\
                              'orientation', 'matrix size [1]','matrix size [2]', \
                             'matrix size [3]','matrix size [4]']
        self._float_scalars=['scaling factor (mm/pixel) [1]','scaling factor (mm/pixel) [2]',\
                             'scaling factor (mm/pixel) [3]']
        self._float_arrays=['mmppix','center']
        
    def readIFH(self,root):
        ifhfile=root+'.4dfp.ifh'
        v,ins,fs,fa=self._vals,self._int_scalars,self._float_scalars,self._float_arrays
        with open(ifhfile,'r') as ifh:
            for line in ifh:
                keyval=line.split(':=')
                key,val=' '.join(keyval[0].strip().split(' ')),keyval[1].strip()
                #print(key,val)
                if key in ins:
                    v[key]=int(val)
                elif key in fs:
                    v[key]=float(val)
                elif key in fa:
                    v[key]=[float(s) for s in val.split()]
                else:
                    v[key]=val              
        #print (v)

    def __str__(self):
        v=self._vals
        with io.StringIO() as out:
            for key in v.keys(): print(key,v[key],file=out)
            return out.getvalue()
        
    def writeIFH(self,root):
        v,ins,fs,fa=self._vals,self._int_scalars,self._float_scalars,self._float_arrays        
        ifhfile=root+'.4dfp.ifh'
        print('writing',ifhfile)
        with open(ifhfile,'w') as ifh:
            for key in v.keys():
                val=v[key]
                if key in ins:
                    ifh.write("{} := {}\n".format(key,val))
                elif key in fs:
                    ifh.write("{} := {:f}\n".format(key,val))
                elif key in fa:
                    ifh.write("{} := {}\n".format(key, ' '.join([str(f) for f in val]) ) )
                else:
                    ifh.write("{} := {}\n".format(key,val))
            
class _4DFP:
    def __init__(self,root=None):
        self.analyze_image=None
        self.ifh=IFH()
        self.infile_list=[]
        if root is not None:
            self.read(root)
            self.infile_list=[root]
            
    def get_copy(self):
        sai=self.analyze_image
        newim=_4DFP()
        newim.analyze_image=nib.AnalyzeImage(get_voxels().copy(),sai.affine,sai.header)
        newim.ifh=self.ifh
        newim.infile_list=self.infile_list
        return newim
        
    def get_voxels(self):
        return self.analyze_image.get_fdata() if self.analyze_image is not None else None
    
    def set_voxels(self,new_voxels):
        ai=self.analyze_image
        self.analyze_image=nib.AnalyzeImage(new_voxels,ai.affine,ai.header)
        
    def read(self, root):
        try:
            print('reading '+root)
            self.analyze_image=nib.AnalyzeImage.load(root+".4dfp.img")
            self.ifh.readIFH(root)
        except Exception as e:
            print (e.strerror)
    
    def write(self, root, infile_list=[]):
        try:
            print('writing',root+'.4dfp.img')
            self.analyze_image.to_filename(root+'.4dfp.img')
            self.ifh.writeIFH(root)
            write_rec_file_4dfp(root,self.infile_list+infile_list)
            
        except Exception as e:
            print(e)

def get_rec_file_root(root,main_extension):
    '''
    Construct .rec file from input file
    root: input file with or without extension
    main_extension: input file extension to be included in rec file name, e.g. 'nii'
    '''
    print('get_rec_file_root',root)
    if main_extension is None:
        return root+'.rec'
    
    out,ext=root,main_extension.replace('.','')
    if not root.endswith('.'+ext):
        out+='.'+ext
    return out+'.rec'


def get_inlog(root, inlogs):
    '''
    Read input file log, if exists
    root: input file
    inlogs: dict containing all input logs
    '''
    recfile=root+'.rec'
    if os.path.isfile(recfile):
        inlogs[root]=cmdprov.read_log(recfile)    
    else:
        print('log file {} does not exist'.format(recfile))
        
def get_inlog_4dfp(root, inlogs):
    '''
    Read input file log, if exists
    root: input file
    inlogs: dict containing all input logs
    '''    
    if os.path.isfile(root+'.rec'):
        inlogs[root]=cmdprov.read_log(root+'.rec')
    elif os.path.isfile(root+'.4dfp.rec'):
        inlogs[root]=cmdprov.read_log(root+'.4dfp.rec')
    elif os.path.isfile(root+'.4dfp.img.rec'):
        inlogs[root]=cmdprov.read_log(root+'.4dfp.img.rec')    
    else:
        print('log file with root {} and rec extension does not exist'.format(root))
        
def write_rec_file_4dfp(root:str,infiles=[]):
    '''
    Write log file for input extension.
    root: input file root
    main_extension: input file extension
    infiles: list of all input files that may have logs attached
    '''    
    inlogs={}
    for file in infiles:
        get_inlog_4dfp(file,inlogs) 
    outlog=cmdprov.new_log(infile_logs=inlogs,extra_notes=["user: "+getpass.getuser(),"node: "+os.uname()[1]])    
    outlog_file=root+'.4dfp.img.rec'    
    print('writing',outlog_file)
    cmdprov.write_log(outlog_file,outlog)
    
def write_rec_file(root:str,main_extension=None,infiles=[]):
    '''
    Write log file for input extension.
    root: input file root
    main_extension: input file extension
    infiles: list of all input files that may have logs attached
    '''
    
    inlogs={}
    for file in infiles:
        get_inlog(file,inlogs)
    outlog=cmdprov.new_log(infile_logs=inlogs,extra_notes=["user: "+getpass.getuser(),"node: "+os.uname()[1]])
    outlog_file=get_rec_file_root(root,main_extension)
    print('writing',outlog_file)
    cmdprov.write_log(outlog_file,outlog)