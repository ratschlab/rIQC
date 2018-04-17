import pysam
import pdb

import os
import scipy as sp
import fnmatch
import gzip
import pandas


def readExpData(fn, mode = 'raw'):

    assert mode in ['rpkm','raw'], 'Function is called with unknown mode'
    if fn.endswith('gz'):
        IN     = gzip.open(fn ,'r')
    else:
        IN     = open(fn, 'r')
    header = IN.readline()
    IN.close()

    header  = header.strip('\n').split('\t')
    if mode == 'rpkm':
        myCols = range(3,len(header),3) 
    elif mode == 'raw':
        print "CAREFUL: You chose mode 'raw', so all results will be based on countst only"
        myCols = range(1,len(header), 3)
    exonpos = sp.array(pandas.read_csv(fn, delim_whitespace = True, usecols = [0], skiprows =2, engine = 'c'))#sp.loadtxt(fn, delimiter = '\t', dtype = 'string', usecols = [0], skiprows = 2))
    data    = sp.array(pandas.read_csv(fn, delim_whitespace = True, dtype = 'float',  usecols = myCols, skiprows = 2, engine = 'c'))#sp.loadtxt(fn, delimiter = '\t', dtype = 'float',  usecols = myCols, skiprows = 2))
    header  = sp.array(header[3::3])
    sidx    = sp.argsort(header)
    header  = header[sidx]
    data    = data[:,sidx]

    return exonpos.astype('string').ravel(), header, data

def readExpDataBam(base_dir):
    # base_dir = "/cbio/grlab/projects/TCGA/PanCancer/icgc_qc"
    allfiles = os.listdir(base_dir)
    allfiles = fnmatch.filter(allfiles, '*.tsv')
    for i,f in enumerate(allfiles):
        if i == 0:
            header  = sp.array([f.split('.')[0]])
            data    = sp.array(pandas.read_csv(os.path.join(base_dir, f), delim_whitespace = True))#sp.loadtxt(os.path.join(base_dir, f), delimiter = '\t', dtype = 'string')
            exonpos = data[:,0]#data[:,[0,1]].ravel('C')
            data    = data[:,1].astype('float')#data[:,[6,7]].ravel('C').astype('float')
        else:
            header = sp.hstack((header, sp.array([f.split('.')[0]])))
            tmp    = sp.array(pandas.read_csv(os.path.join(base_dir, f), delim_whitespace = True))#sp.loadtxt(os.path.join(base_dir, f), delimiter = '\t', dtype = 'string')
            tmp    = tmp[:,1].astype('float')#tmp[:,[6,7]].ravel('C').astype('float')
            data   = sp.vstack((data, tmp))
    if len(data.shape) == 1:
        data = data[:, sp.newaxis]
    else:
        data    = data.T
    sidx    = sp.argsort(header)
    header  = header[sidx]
    data    = data[:, sidx]

    ### remove non chromosomal contigs
    iOK     = sp.array([x.startswith('chr') for x in exonpos])
    exonpos = exonpos[iOK]
    data    = data[iOK,:]

    chrm  = sp.array([x.split(':')[0].strip('chr') for x in exonpos])
    start = sp.array([x.split(':')[1].split('-')[0] for x in exonpos]).astype('int')
    sidx  = sp.lexsort((start, chrm))
    
    exonpos = exonpos[sidx]
    data    = data[sidx,:]

    return exonpos, header, data
