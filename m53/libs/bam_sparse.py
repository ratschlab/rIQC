import pysam
import time
import scipy as sp
import scipy.sparse as spst
import sys
import os
import warnings
import pdb
import h5py
from sets import Set

def get_counts_from_single_bam_sparse(fn_bam, regions):

    """This function extracts read counts from a given sparse bam file spanning 
       a set of given intervals."""

    if not os.stat(fn_bam).st_size > 0:
        warnings.warn('WARNING: alignment file %s seems to be empty and will be skipped! \n' % fn_bam)
        dummy = sp.zeros(regions.shape[0] * 2)
        dummy[:] = sp.nan
        return dummy
        
    if fn_bam.lower().endswith('npz'):
        IN = sp.load(fn_bam)
    else:
        IN = h5py.File(fn_bam, 'r')
    refseqs = sp.unique([x.split('_')[0] for x in IN.keys()])
    cnts = sp.zeros((regions.shape[0], 2), dtype='float')
    t0 = time.time()
    
    if len(regions.shape) > 1:
        sidx = sp.argsort(regions[:, 0])
    else:
        sidx = sp.argsort(regions)

    last_chrm = ''
    
    for i, ii in enumerate(sidx):
        rec = regions[ii]
        if i > 0 and i % 100 == 0:
            print '%i rounds to go. ETA %.0f seconds' % (regions.shape[0] - i, (time.time() - t0) / i * (regions.shape[0] - i))
        if len(regions.shape) == 1:
            chrm = rec.split(':')[0]
            if not chrm in refseqs:
                chrm = chrm.strip('chr')
            start1 = int(rec.split(':')[1].split('-')[0])
            end1   = int(rec.split(':')[1].split('-')[1])
            start2 = None
            end2   = None
        else:
            chrm = rec[0].split(':')[0] 
            if not chrm in refseqs:
                chrm = chrm.strip('chr')
            start1 = int(rec[0].split(':')[1].split('-')[0])
            end1   = int(rec[0].split(':')[1].split('-')[1])
            start2 = int(rec[1].split(':')[1].split('-')[0])
            end2   = int(rec[1].split(':')[1].split('-')[1])
        try:
            if last_chrm == '' or chrm != last_chrm:
                cache = spst.coo_matrix((IN[chrm + '_reads_dat'][:], (IN[chrm + '_reads_row'][:], IN[chrm + '_reads_col'][:])), shape=IN[chrm + '_reads_shp'][:], dtype='uint32').tocsc()
                last_chrm = chrm
            cnt1 = int(sp.ceil(sp.sum(cache[0, start1:end1].todense()) / 50.0))
            if start2 is None:
                cnt2 = cnt1
            else:
                cnt2 = int(sp.ceil(sp.sum(cache[0, start2:end2].todense()) / 50.0))
            #print '%s\t%s\tcnt1: %i\tcnt2: %i' % (rec[0], rec[1], cnt1, cnt2)
        except:
            print >> sys.stderr, 'Ignored %s' % chrm
            cnt1 = 1
            cnt2 = 1
        finally:
            #cnts.append([cnt1, cnt2])
            cnts[ii, :] = [cnt1, cnt2]
    IN.close()

    return cnts.ravel('C')
    #return sp.array(cnts, dtype='float').ravel('C')

def get_counts_from_multiple_bam_sparse(fn_bams, regions):
    """ This is a wrapper to concatenate counts for a given list of bam
        files"""
    
    if len(fn_bams) == 1:
        return get_counts_from_single_bam_sparse(fn_bams[0], regions)[:, sp.newaxis]
    else:
        return sp.hstack([get_counts_from_single_bam_sparse(fn_bams[i], regions)[:,sp.newaxis] for i in range(len(fn_bams))])
