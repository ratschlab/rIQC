import pysam
import time
import scipy as sp
import numpy as np
import sys
import os
import warnings
import pdb
from sets import Set


def get_counts_from_single_bam(fn_bam, regions):
    """This function extracts read counts from a given bam file spanning
       a set of given intervals."""

    if not os.path.exists(fn_bam + '.bai'):
        # raise Exception('\nERROR: alignment file %s seems not to be indexed\n' % fn_bam)
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
    cnts = dict()
    t0 = time.time()

    # Sort regions by chr
    if len(regions.shape) <= 1:
        exit("regions.shape was not > 1")
    sidx = sp.argsort(regions[:, 1])  # seq_name

    for i, ii in enumerate(sidx):
        if i > 0 and i % 100 == 0:
            print '%i rounds to go. ETA %.0f seconds' % (
            regions.shape[0] - i, (time.time() - t0) / i * (regions.shape[0] - i))

        rec = regions[ii, 0:3]
        rec_exons = regions[ii, 4].split(",")
        exon_data = np.zeros((len(rec_exons), 4), dtype=float)

        chrm = rec[1]

        if chrm not in refseqs:
            chrm = chrm.strip('chr')
            if chrm not in refseqs:
                exit("%s is not in bam-references" % chrm)

        if len(regions.shape) == 1:
            exit("regions.shape == 1")

        for e in range(len(rec_exons)):
            start = int(rec_exons[e].split('-')[0])
            end = int(rec_exons[e].split("-")[1])
            exon_data[e, 0] = start
            exon_data[e, 1] = end
            exon_data[e, 2] = end - start
            cnt = 1
            try:
                cnt = int(sp.ceil(sp.sum(
                    [sp.sum((sp.array(read.positions) >= start) & (sp.array(read.positions) < end)) for read in
                    samfile.fetch(str(chrm), start, end) if not read.is_secondary]) / 50.0))
            except ValueError:
                print >> sys.stderr, 'Ignored %s' % chrm
            finally:
                exon_data[e, 3] = cnt / exon_data[e, 2]
        cnts[rec[0]] = exon_data
    samfile.close()

    return cnts
    # return sp.array(cnts, dtype='float').ravel('C')


def get_counts_from_multiple_bam(fn_bams, regions):
    """ This is a wrapper to concatenate counts for a given list of bam
        files"""

    if len(fn_bams) == 1:
        return [get_counts_from_single_bam(fn_bams[0], regions)]
    else:
        li = []
        for i in range(len(fn_bams)):
            li.append(get_counts_from_single_bam(fn_bams[i], regions))
        return li
