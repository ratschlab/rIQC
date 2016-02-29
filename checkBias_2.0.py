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

import logging

### may not be necessary

def parse_options(argv):

    parser = OptionParser()
    
    sampleinput = OptionGroup(parser, 'Input')
    sampleinput.add_option('-f', '--file', dest='fn_exonq', metavar= 'FILE', help='Exon quantification file from firehose', default = '-')
    sampleinput.add_option('-b', '--bam_dir', dest = 'dir_bam', metavar = 'FILE', help = 'Directory of bam files', default = '-')
    sampleinput.add_option('-t', '--tab_cnt', dest = 'dir_cnt', metavar = 'FILE', help = 'Directory of tab delimited count files' , default = '-')
    sampleinput.add_option('-a', '--fn_anno', dest = 'fn_anno', metavar = 'FILE', help = 'Annotation', default = '-')
    sampleinput.add_option('-n', '--fn_bam', dest = 'fn_bam', metavar = 'FIlE', help = 'Specifies bam file for counting only', default = '-')

    #### optional arguments
    optional = OptionGroup(parser, 'Options')    
    optional.add_option('-q','--quant', dest = 'qmode', metavar = 'STRING', help = 'What type of quantification to use [rpkm,raw]', default = 'raw')
    optional.add_option('-l', '--length',    dest = 'length',   metavar = 'STRING', help = 'Length filter [uq,mq,lq]', default = 'uq')
    optional.add_option('-w', '--whitelist', dest = 'fn_white', metavar = 'FILE',   help = 'White list with ids',      default = '-')
    optional.add_option('-o', '--fnout',    dest = 'fn_out',   metavar = 'FILE',   help = 'prefix for output',        default = 'out')
    optional.add_option('-g', '--log',       dest = 'fn_log',   metavar = 'FILE',   help = 'Log file',                 default = 'out.log')
    optional.add_option('-p', '--plot', dest = 'doPlot', action = "store_true", help = 'Plot figures', default=False)
    optional.add_option('-v', '--verbose', dest = 'isVerbose', action = "store_true", help = 'Set Logger To Verbose', default=False)
    optional.add_option('-c', '--pseudocount', dest = 'doPseudo', action = "store_true", help = 'Add Pseudocounts to ratio', default=False)
    optional.add_option('-m', '--fn_anno_tmp', dest = 'fn_anno_tmp', metavar = 'FILE', help = 'Temp file for storing anno info', default = os.path.join(os.path.realpath(__file__).rsplit('/',1)[:-1][0] ,'anno.tmp'))
    optional.add_option('-i', '--genelist', dest = 'fn_genes', metavar = 'FILE', help = 'file with genenames to use', default = '-')
    optional.add_option('-s', '--fn_sample_ratio', dest = 'fn_sample_ratio', metavar = 'FILE', help = 'Sample Ratios in relation to yours', default = os.path.join(os.path.realpath(__file__).rsplit('/',1)[:-1][0] ,'data','sampleRatios/TCGA_sample_a_ratio_uq.tsv'))
    optional.add_option('-d', '--mask-filter', dest = 'filt', help = 'Mask all readcounts below this integer', default = '0')
    
    parser.add_option_group(sampleinput)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)
    if sp.sum(int(options.fn_exonq != '-') + int(options.dir_bam != '-') + int(options.dir_cnt != '-') + int(options.fn_bam != '-')) != 1:
        print "Please specify exactly one type of input file(s) (e.g.: Exon quantification, Bam Files or Tab delimited count files)"
        parser.print_help()
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
    myLength = sp.zeros(exonTgene.shape[0])

    for i,rec in enumerate(exonTgene):

        istart = exonpos == rec[0]
        iend   = exonpos == rec[1]
        if sp.sum(istart) == 0:
            continue
        if sp.sum(iend) == 0:
            continue
        if exonpos[istart][0].split(':')[-1] == '-':
            tmp    = istart
            istart = iend
            iend   = tmp

        startcc  = sp.array(exonpos[istart][0].split(':')[1].split('-')).astype('int')
        endcc    = sp.array(exonpos[iend][0].split(':')[1].split('-')).astype('int')

        startlen = startcc[1] - startcc[0]
        endlen   = endcc[1]   - endcc[0]
       
        myLength[i]     = float(rec[4])
        mycounts[i,:,0] = data[istart,:]
        mycounts[i,:,1] = data[iend,:]  
    return mycounts, myLength


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

    elif options.dir_bam != '-':
        bam_list = glob.glob(os.path.join(options.dir_bam, '*.bam'))
        header = bam_list ### change this TODO 
        data = get_counts_from_multiple_bam(bam_list, exonTgene) ### REMOVE
        exonpos = exonTgene[:, :2].ravel('C')
    elif options.fn_bam != '-':
        print "WARNING: Running only gene counts"
        exonTable = getFullAnnotationTable()
        data = get_counts_from_single_bam(options.fn_bam,exonTable)
        sp.savetxt(options.fn_out+'counts.tsv', sp.vstack((exonTable,data[::2])).T, delimiter = '\t', fmt = '%s')
        sys.exit(0)

    ### normalize counts by exon length
    logging.info("Normalize counts by exon length")
    if (options.fn_exonq == '-') | (options.qmode == 'raw'):
        exonl = sp.array([int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) + 1 for x in exonpos])
        data /= sp.tile(exonl[:, sp.newaxis], data.shape[1])
        #data /= sp.hstack([exonl[:, sp.newaxis] for x in xrange(data.shape[0] / exonl.shape[0])]).ravel(1)  



    ### Subset to whitelist
    if options.fn_white != '-': 
        logging.info("Subsetting to whitelist")
        header, data = whitelisting(options, header, data)


    ### Calculate 3'/5' Bias
    logging.info("Calculate Bias")
    mycounts, myLength = calculateBias(exonTgene, data, exonpos)


    ### subset to high expression ### TODO: need to change this for clarity here....
    logging.info("Make sure I got only reasonably expressed genes")
    if (options.fn_genes == '-') & (options.fn_exonq != '-'): ### assuming that i do not have rpkm and not pre-selected genes anyways
        primeCov = sp.mean(mycounts[:,:,0], axis = 1)   + sp.mean(mycounts[:,:,1], axis = 1)     
        ### ensure average expression of 1 rpkm across samples
        if options.length == 'uq':
            iOK  = (sp.mean(mycounts[:,:,0], axis = 1) > 1) & (sp.mean(mycounts[:,:,1], axis = 1) > 1) 
        elif options.length == 'mq':
            iOK  = (sp.mean(mycounts[:,:,0], axis = 1) > 1) & (sp.mean(mycounts[:,:,1], axis = 1) > 1)
        elif options.length == 'lq':
            iOK  = (sp.mean(mycounts[:,:,0], axis = 1) > 1) & (sp.mean(mycounts[:,:,1], axis = 1) > 1)
        mycounts = mycounts[iOK,:,:]
        myLength = myLength[iOK]
        sp.savetxt(options.fn_out+'.geneSet', exonTgene[iOK,:], fmt = '%s', delimiter = '\t')


    if options.doPseudo:
        logging.info("Add Pseudocount and estimate ratios")
        ratio    = ((mycounts[:,:,1] + 1) / (mycounts[:,:,0] + 1))
    else:
        logging.info("Estimate ratios")
        ratio    = ((mycounts[:,:,1])     / (mycounts[:,:,0]))

    logging.info("Find Median")
    vals = []
    for i in xrange(mycounts.shape[1]):        
        iOK = ~(sp.isnan(mycounts[:,i,0])) & ~(sp.isnan(mycounts[:,i,1]))
        tmp = ((mycounts[:,i,1] )[iOK]) / ((mycounts[:,i,0] )[iOK] )
        vals.append(sp.percentile(tmp[~sp.isnan(tmp)],50))
    vals = sp.array(vals)



    sidx   = sp.argsort(vals)
    iqr    = ( (sp.percentile(vals,75) - sp.percentile(vals,25) ) * 1.5)

    logging.info("Tukey Filter is estimated to be %f" % (iqr+sp.percentile(vals, 75)))
    print "Tukey Filter is estimated to be %f" % (iqr+sp.percentile(vals, 75))
    print "Tukey Filter is estimated to be %f" % (sp.percentile(vals, 25)-iqr)

    sp.savetxt('%s_sample_a_ratio_%s.tsv' % (options.fn_out,options.length), sp.vstack((header, vals.astype('string'))).T, delimiter = '\t', fmt = '%s')
    
    ratio = ratio[:,sidx]  
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
