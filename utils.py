import argparse, getpass, cmdline_provenance as cmdprov, os.path

def get_rec_file_root(root,main_extension):
    '''
    Construct .rec file from input file
    root: input file with or without extension
    main_extension: input file extension to be included in rec file name, e.g. 'nii'
    '''
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