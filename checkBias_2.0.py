import scipy as sp
import scipy.stats as spst
import numpy as np


import sys
import os
import pdb
import glob
import fnmatch
from optparse import OptionParser, OptionGroup

from libs.annotation import *
from libs.counts import *
from libs.viz import *
from libs.bam import *
from libs.bam_sparse import *
from libs.kmer import *

import logging

### may not be necessary

def parse_options(argv):

    parser = OptionParser()
    
    sampleinput = OptionGroup(parser, 'Input')
    sampleinput.add_option('-f', '--file', dest='fn_exonq', metavar= 'FILE', help='Exon quantification file from firehose', default = '-')
    sampleinput.add_option('-b', '--bam_dir', dest = 'dir_bam', metavar = 'FILE', help = 'Directory of bam files', default = '-')
    sampleinput.add_option('-t', '--tab_cnt', dest = 'dir_cnt', metavar = 'FILE', help = 'Directory of tab delimited count files' , default = '-')
    sampleinput.add_option('-F', '--fastq_dir', dest = 'fastq_dir', metavar = 'FILE', help = 'Directory of fastq files', default = '-')
    sampleinput.add_option('-a', '--fn_anno', dest = 'fn_anno', metavar = 'FILE', help = 'Annotation', default = '-')
    sampleinput.add_option('-n', '--fn_bam', dest = 'fn_bam', metavar = 'FIlE', help = 'Specifies singl bam file', default = '-')
    sampleinput.add_option('-w', '--whitelist', dest = 'fn_white', metavar = 'FILE',   help = 'White list with ids',      default = '-')
    sampleinput.add_option('-G', '--genome', dest = 'fn_genome', metavar = 'FILE', help = 'Path to genome file in fasta', default = '-')
    sampleinput.add_option('-i', '--genelist', dest = 'fn_genes', metavar = 'FILE', help = 'file with genenames to use', default = '-')

    sampleoutput = OptionGroup(parser, 'Output')
    sampleoutput.add_option('-o', '--fnout',    dest = 'fn_out',   metavar = 'FILE',   help = 'prefix for output',        default = 'out')
    sampleoutput.add_option('-m', '--fn_anno_tmp', dest = 'fn_anno_tmp', metavar = 'FILE', help = 'Temp file for storing anno info', default = os.path.join(os.path.realpath(__file__).rsplit('/',1)[:-1][0] ,'anno.tmp'))

    opt_gen = OptionGroup(parser, 'General Options')    
    opt_gen.add_option('-q', '--quant', dest = 'qmode', metavar = 'STRING', help = 'What type of quantification to use [rpkm,raw]', default = 'raw')
    opt_gen.add_option('-c', '--pseudocount', dest = 'doPseudo', action = "store_true", help = 'Add Pseudocounts to ratio', default=False)
    opt_gen.add_option('-C', '--count_only', dest = 'count_only', action = "store_true", help = 'Only do counting on given input [off]', default=False)
    opt_gen.add_option('-l', '--length',    dest = 'length',   metavar = 'STRING', help = 'Length filter [uq,mq,lq]', default = 'uq')
    opt_gen.add_option('-g', '--log',       dest = 'fn_log',   metavar = 'FILE',   help = 'Log file',                 default = 'out.log')
    opt_gen.add_option('-v', '--verbose', dest = 'isVerbose', action = "store_true", help = 'Set Logger To Verbose', default=False)
    opt_gen.add_option('', '--sparse_bam', dest = 'sparse_bam', action = "store_true", help = 'Input BAM files are in sparse hdf5 format [off]', default=False)
    opt_gen.add_option('-p', '--plot', dest = 'doPlot', action = "store_true", help = 'Plot figures', default=False)
    opt_gen.add_option('-s', '--fn_sample_ratio', dest = 'fn_sample_ratio', metavar = 'FILE', help = 'Sample Ratios in relation to yours', default = os.path.join(os.path.realpath(__file__).rsplit('/',1)[:-1][0] ,'data','sampleRatios/TCGA_sample_a_ratio_uq.tsv'))
    opt_gen.add_option('-d', '--mask-filter', dest = 'filt', help = 'Mask all readcounts below this integer', default = '0')

    opt_kmer = OptionGroup(parser, 'Options for k-mer counting')
    opt_kmer.add_option('-k', '', dest = 'k', type = 'int', help = 'Length of k-mer for alignmentfree counting [27]', default = 27)
    opt_kmer.add_option('-R', '--reads_kmer', dest = 'kmer_thresh', type = 'float', help = 'Required active reads per sample [50000] / if btw 0 and 1 fraction of input reads considered', default = 50000)
    opt_kmer.add_option('-S', '--step_k', dest = 'step_k', type = 'int', help = 'Step-size for k-mer counting [4]', default = 4)
    
    parser.add_option_group(sampleinput)
    parser.add_option_group(sampleoutput)
    parser.add_option_group(opt_gen)
    parser.add_option_group(opt_kmer)
    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)
    if sp.sum(int(options.fn_exonq != '-') + int(options.dir_bam != '-') + int(options.dir_cnt != '-') + int(options.fn_bam != '-') + int(options.fastq_dir != '-')) != 1:
        print "Please specify exactly one type of input file(s) (e.g.: Exon quantification, Bam Files or Tab delimited count files)"
        parser.print_help()
        sys.exit(2)
    if options.fastq_dir != '-' and options.fn_genome == '-':
        print >> sys.stderr, 'For usage on fastq files a genome file in fasta needs to be provided via -G/--genome'
        sys.exit(2)
    return options
    

def whitelisting(options, header, data):
    whitelist = sp.loadtxt(options.fn_white, delimiter = '\t', dtype = 'string')
    midx_m    = sp.in1d(header, whitelist)
    tags      = sp.array([x.split('-')[3] for x in header])
    midx_n    = np.core.defchararray.startswith(tags, '1')        
    header    = header[midx_m | midx_n]
    data      = data[:, midx_m | midx_n]
    return header, data


def calculateBias(exonTgene, data, exonpos):
    mycounts = sp.zeros((exonTgene.shape[0], data.shape[1], 2))

    for i,rec in enumerate(exonTgene):

        istart = exonpos == rec[0]
        iend   = exonpos == rec[1]
        if sp.sum(istart) == 0:
            continue
        if sp.sum(iend) == 0:
            continue
        if exonpos[istart][0].split(':')[-1] == '-' and int(exonpos[istart][0].split(':')[1].split('-')[0]) < int(exonpos[iend][0].split(':')[1].split('-')[0]):
            istart,iend = iend,istart
      
        mycounts[i,:,0] = data[istart,:]
        mycounts[i,:,1] = data[iend,:]  
    return mycounts


def main():
    ### Parse options
    options = parse_options(sys.argv)
    filt = int(options.filt)


    #### set up logger
    logging.basicConfig(filename = options.fn_log, level = 0, format = '%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger()
    if options.isVerbose:
        consoleHandler = logging.StreamHandler()
        log.addHandler(consoleHandler)
    ### Read annotation from file
    logging.info("Reading Annotation from file")
    exonTgene = getAnnotationTable(options)

    ### Read in or generate quantifications
    logging.info("Reading in quantifications")
    if options.fn_exonq != '-': ### TODO: change this into raw by principle and then generate rpkm which match firebrowse
        exonpos,  header, data  = readExpData(options.fn_exonq, options.qmode)        

        iOK = ~sp.array([x.startswith('chrM') for x in exonpos])
        exonpos = exonpos[iOK]
        data    = data[iOK,:]
        sidx    = sp.argsort(exonpos)
        exonpos = exonpos[sidx]
        data    = data[sidx,:]


    elif options.dir_cnt != '-':
        exonpos, header, data   = readExpDataBam(options.dir_cnt) ### move this over to file lists rather than dirs
        data[data < filt]       = sp.nan
        if options.qmode == 'rpkm':
            data = (data * 1E9) / sp.sum(data , axis = 0)
            exonl = sp.array([int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) + 1 for x in exonpos])
            data /= sp.tile(exonl[:, sp.newaxis], data.shape[1])

    elif options.fastq_dir != '-':
        kmers1, kmers2 = prepare_kmers(options, exonTgene) 
        kmers1, kmers2 = clean_kmers(options, kmers1, kmers2)
        fastq_list = glob.glob(os.path.join(options.fastq_dir))
        header = fastq_list
        data = get_counts_from_multiple_fastq(fastq_list, kmers1, kmers2, options)
        exonpos = exonTgene[:, :2].ravel('C')
    elif options.dir_bam != '-':
        if options.sparse_bam:
            bam_list = glob.glob(os.path.join(options.dir_bam, '*.hdf5'))
            header = bam_list ### change this TODO 
            data = get_counts_from_multiple_bam_sparse(bam_list, exonTgene) ### REMOVE
        else:
            bam_list = glob.glob(os.path.join(options.dir_bam, '*.bam'))
            header = bam_list ### change this TODO 
            data = get_counts_from_multiple_bam(bam_list, exonTgene) ### REMOVE
        exonpos = exonTgene[:, :2].ravel('C')
    elif options.fn_bam != '-':
        if options.count_only:
            print "WARNING: Running only gene counts"
            exonTable = sp.sort(exonTgene[:,[0,1]].ravel())        
            data = get_counts_from_single_bam(options.fn_bam,exonTable)
            sp.savetxt(options.fn_out+'counts.tsv', sp.vstack((exonTable,data[::2])).T, delimiter = '\t', fmt = '%s')
            sys.exit(0)
        else:
            header = [options.fn_bam] ### change this TODO 
            if options.sparse_bam:
                data = get_counts_from_multiple_bam_sparse([options.fn_bam], exonTgene) ### REMOVE
            else:
                data = get_counts_from_multiple_bam([options.fn_bam], exonTgene) ### REMOVE
            exonpos = exonTgene[:, :2].ravel('C')


    ### normalize counts by exon length
    logging.info("Normalize counts by exon length")
    if (options.fn_exonq == '-') | (options.qmode == 'raw'):
        exonl = sp.array([int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) + 1 for x in exonpos], dtype='float') / 1000.
        data /= sp.tile(exonl[:, sp.newaxis], data.shape[1])
        #data /= sp.hstack([exonl[:, sp.newaxis] for x in xrange(data.shape[0] / exonl.shape[0])]).ravel(1)  


    ### Subset to whitelist
    if options.fn_white != '-': 
        logging.info("Subsetting to whitelist")
        header, data = whitelisting(options, header, data)


    ### Calculate 3'/5' Bias
    logging.info("Calculate Bias")
    mycounts = calculateBias(exonTgene, data, exonpos)

    ### subset to high expression ### TODO: need to change this for clarity here....
    logging.info("Make sure I got only reasonably expressed genes")
    if (options.fn_genes == '-') & (options.fn_exonq != '-'): ### assuming that i do not have rpkm and not pre-selected genes anyways
        primeCov = sp.mean(mycounts[:,:,0], axis = 1)   + sp.mean(mycounts[:,:,1], axis = 1)     

        ### ensure average expression of 1 rpkm across samples
        iOK  = (sp.mean(mycounts[:,:,0], axis = 1) > 1) & (sp.mean(mycounts[:,:,1], axis = 1) > 1) 
        mycounts = mycounts[iOK,:,:]
        myLength = myLength[iOK]
        sp.savetxt(options.fn_out+'.geneSet', exonTgene[iOK,:], fmt = '%s', delimiter = '\t')

    logging.info("Find Median")
    vals = []
    for i in xrange(mycounts.shape[1]):        
        if options.doPseudo:
            ### AK: I had to filter for counts with some signal, otherwise the median is always 1.0 ...
            iOK = ((mycounts[:, i, 1] > 0) | (mycounts[:, i, 0] > 0))
            ratio = sp.percentile((mycounts[iOK, i, 1] + 1)  / (mycounts[iOK, i, 0] + 1), 50)
        else:
            ### AK: This version allows for the ocurrence of infs in some cases where there should be nan's 
            #iOK = ~(sp.isnan(mycounts[:,i,0])) & ~(sp.isnan(mycounts[:,i,1]))
            #tmp = ((mycounts[:,i,1] )[iOK]) / ((mycounts[:,i,0] )[iOK] )
            #vals.append(sp.percentile(tmp[~sp.isnan(tmp)],50))

            ### AK: This is my new version
            iOK = ((mycounts[:, i, 1] > 0) & (mycounts[:, i, 0] > 0))
            #iOK = ~(sp.isnan(mycounts[:, i, 0])) & ~(sp.isnan(mycounts[:, i, 1])) & (mycounts[:, i, 0] > 0)
            ratio = sp.percentile(mycounts[iOK, i, 1] / mycounts[iOK, i, 0], 50)
        assert sp.sum(sp.isnan(ratio)) + sp.sum(sp.isinf(ratio)) == 0
        vals.append(ratio)
    vals = sp.array(vals)

    sidx   = sp.argsort(vals)
    iqr    = ((sp.percentile(vals,75) - sp.percentile(vals,25)) * 1.5)

    logging.info("Tukey Filter is estimated to be %f" % (iqr+sp.percentile(vals, 75)))
    if len(vals) > 1:
        print "Tukey Filter is estimated to be %f" % (iqr+sp.percentile(vals, 75))
        print "Tukey Filter is estimated to be %f" % (sp.percentile(vals, 25)-iqr)

    sp.savetxt('%s_sample_a_ratio_%s.tsv' % (options.fn_out,options.length), sp.vstack((header, vals.astype('string'))).T, delimiter = '\t', fmt = '%s')
    
    if options.doPlot:
        logging.info("Plot all samples")

        baselinedata = sp.loadtxt(options.fn_sample_ratio, delimiter = '\t', dtype = 'string')
        baselinedata = baselinedata[:,1].astype('float')

        basePval = sp.hstack((baselinedata, vals))
        midx     = sp.hstack((sp.ones(baselinedata.shape[0]), sp.zeros(vals.shape[0]))).astype('bool')
        plotBias(basePval, '%s_bias_sorted_vline_%s.png' % (options.fn_out,options.length), midx)
        midx     = sp.hstack((sp.ones(baselinedata.shape[0]), sp.zeros(vals.shape[0]))).astype('bool')
        plotBias(basePval, '%s_bias_sorted_vline_log_%s.png' % (options.fn_out,options.length), midx, logScale = True)


if __name__ == "__main__":
    main()
