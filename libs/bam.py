import pysam
import time
import scipy as sp
import sys
import os
import warnings
import pdb
from sets import Set

def get_counts_from_single_bam(fn_bam, regions):

    """This function extracts read counts from a given bam file spanning 
       a set of given intervals."""

    if not os.path.exists(fn_bam + '.bai'):
        #raise Exception('\nERROR: alignment file %s seems not to be indexed\n' % fn_bam)
        warnings.warn('WARNING: alignment file %s seems not to be indexed and will be skipped! \n' % fn_bam)
        dummy = sp.zeros(regions.shape[0] * 2)
        dummy[:] = sp.nan
        return dummy
    if not os.stat(fn_bam).st_size > 0:
        warnings.warn('WARNING: alignment file %s seems to be empty and will be skipped! \n' % fn_bam)
        dummy = sp.zeros(regions.shape[0] * 2)
        dummy[:] = sp.nan
        return dummy
        
    samfile = pysam.Samfile(fn_bam, 'rb')
    refseqs = samfile.references
    cnts = []
    t0 = time.time()
    
    
    for i, rec in enumerate(regions):
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
            cnt1    = len([1 for read in samfile.fetch(chrm, start1, end1) if not read.is_secondary]) #Otherwise does not match firebrowse            
            if start2 is None:
                cnt2 = cnt1
            else:
                cnt2    = len([1 for read in samfile.fetch(chrm, start2, end2) if not read.is_secondary]) #Otherwise does not match firebrowse
        except ValueError:
            print >> sys.stderr, 'Ignored %s' % chrm
            cnt1 = 1
            cnt2 = 1
        finally:
            cnts.append([cnt1, cnt2])
    samfile.close()
    return sp.array(cnts, dtype='float').ravel('C')

def get_counts_from_multiple_bam(fn_bams, regions):
    """ This is a wrapper to concatenate counts for a given list of bam
        files"""
    
    if len(fn_bams) == 1:
        return get_counts_from_single_bam(fn_bams[0], regions)[:, sp.newaxis]
    else:
        return sp.hstack([get_counts_from_single_bam(fn_bams[i], regions)[:,sp.newaxis] for i in range(len(fn_bams))])
