import scipy as sp
import scipy.stats as spst
import numpy as np

import sys
import os
import pdb
import glob
import fnmatch
import cPickle
from optparse import OptionParser, OptionGroup
import logging

from libs.annotation import *
from libs.viz import *
from libs.bam import *
from libs.bam_sparse import *
from libs.kmer import *


def parse_options(argv):
    parser = OptionParser()

    sampleinput = OptionGroup(parser, 'Input')

    sampleinput.add_option('', '--bam_dir', dest='dir_bam', metavar='FILE',
                           help='Directory of bam files', default='-')
    sampleinput.add_option('', '--fastq_dir', dest='dir_fastq', metavar='FILE',
                           help='Directory of fastq files', default='-')
    sampleinput.add_option('', '--bam_fn', dest='fn_bam', metavar='FIlE',
                           help='Specifies single bam file', default='-')
    sampleinput.add_option('', '--cnt_dir', dest='dir_cnt', metavar='FILE',
                           help='Directory of pre-produced tab delimited count files', default='-')

    sampleinput.add_option('', '--anno_fn', dest='fn_anno', metavar='FILE',
                           help='Annotation', default='-')

    sampleinput.add_option('', '--genome', dest='fn_genome', metavar='FILE',
                           help='Path to genome file in fasta', default='-')
    sampleinput.add_option('', '--genelist', dest='fn_genes', metavar='FILE',
                           help='file with genenames to use', default='-')
    sampleinput.add_option('', '--separate_files', dest='separate_files', action="store_true",
                           help='Consider all input files individually [off]', default=False)

    sampleoutput = OptionGroup(parser, 'Output')

    sampleoutput.add_option('', '--out_fn', dest='fn_out', metavar='FILE',
                            help='prefix for output', default='out')
    sampleoutput.add_option('', '--anno_tmp_fn', dest='fn_anno_tmp', metavar='FILE',
                            help='Temp file for storing anno info', default='anno.tmp')
    sampleoutput.add_option('', '--pickle_all', dest='fn_pickle_all', metavar='FILE',
                            help='Pickle file for storing all kmers', default=None)
    sampleoutput.add_option('', '--pickle_filt', dest='fn_pickle_filt', metavar='FILE',
                            help='Pickle file for storing filtered/cleaned kmers', default=None)

    opt_gen = OptionGroup(parser, 'General Options')

    opt_gen.add_option('', '--quant', dest='qmode', metavar='STRING',
                       help='What type of quantification to use [rpkm,raw]', default='raw')
    opt_gen.add_option('', '--pseudocount', dest='doPseudo', action="store_true",
                       help='Add Pseudocounts to ratio', default=False)
    opt_gen.add_option('', '--count_only', dest='count_only', action="store_true",
                       help='Only do counting on given input [off]', default=False)
    opt_gen.add_option('', '--length', dest='length', metavar='STRING',
                       help='Length filter [uq,mq,lq]', default='uq')
    opt_gen.add_option('', '--log', dest='fn_log', metavar='FILE',
                       help='Log file', default='out.log')
    opt_gen.add_option('', '--verbose', dest='isVerbose', action="store_true",
                       help='Set Logger To Verbose', default=False)
    opt_gen.add_option('', '--sparse_bam', dest='sparse_bam', action="store_true",
                       help='Input BAM files are in sparse hdf5 format [off]', default=False)
    opt_gen.add_option('', '--plot', dest='doPlot', action="store_true",
                       help='Plot figures', default=False)
    opt_gen.add_option('', '--fn_sample_ratio', dest='fn_sample_ratio', metavar='FILE',
                       help='Sample Ratios in relation to yours',
                       default=os.path.join(os.path.realpath(__file__).rsplit('/', 1)[:-1][0],
                                            'data', 'sampleRatios/TCGA_sample_a_ratio_uq.tsv'))
    opt_gen.add_option('', '--mask-filter', dest='filt',
                       help='Mask all readcounts below this integer', default='0')
    opt_gen.add_option('', '--protein-coding-filter_OFF', dest="protein_coding_filter", action="store_false",
                       help="Consider only genes that are protein-coding", default=True)

    opt_kmer = OptionGroup(parser, 'Options for k-mer counting')

    opt_kmer.add_option('', '--kmer_length', dest='k', type='int',
                        help='Length of k-mer for alignmentfree counting', default=27)
    opt_kmer.add_option('', '--reads_kmer', dest='kmer_thresh', type='float',
                        help='Required active reads per sample or if in [0, 1] then fraction of input reads considered',
                        default=50000)
    opt_kmer.add_option('', '--step_k', dest='step_k', type='int',
                        help='Step-size for k-mer counting', default=4)

    parser.add_option_group(sampleinput)
    parser.add_option_group(sampleoutput)
    parser.add_option_group(opt_gen)
    parser.add_option_group(opt_kmer)
    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)
    if sp.sum(int(options.dir_bam != '-') + int(options.fn_bam != '-')
              + int(options.dir_cnt != '-') + int(options.dir_fastq != '-')) != 1:
        print "Please specify exactly one type of input file(s) (e.g.: Exon quantification, Bam Files)"
        parser.print_help()
        sys.exit(2)
    if options.dir_fastq != '-' and options.fn_genome == '-':
        print >> sys.stderr, 'For usage on fastq files a genome file in fasta needs to be provided via -G/--genome'
        sys.exit(2)
    return options


def __get_counts_from_marginal_exons(exon_t_gene, data):
    mycounts = sp.zeros((exon_t_gene.shape[0], data.shape[1], 2))

    for i, rec in enumerate(exon_t_gene):

        istart = i*2
        iend = i*2 + 1

        if rec[0].split(':')[-1] == '-' and \
                int(rec[0].split(':')[1].split('-')[0]) \
                < int(rec[1].split(':')[1].split('-')[0]):
            istart, iend = iend, istart

        mycounts[i, :, 0] = data[istart, :]
        mycounts[i, :, 1] = data[iend, :]

    return mycounts


def main():
    ### Parse options
    options = parse_options(sys.argv)
    filt = int(options.filt)

    #### set up logger
    logging.basicConfig(filename=options.fn_log, level=0, format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger()
    if options.isVerbose:
        consoleHandler = logging.StreamHandler()
        log.addHandler(consoleHandler)
    ### Read annotation from file
    logging.info("Reading Annotation from file")


    if(options.dir_cnt == '-'):
        exon_t_gene = get_annotation_table(options)

        if options.dir_fastq != '-':
            if(options.fn_pickle_filt != None and os.path.exists(options.fn_pickle_filt)):
                (kmers1, kmers2) = cPickle.load(open(options.fn_pickle_filt, 'r'))
            elif(os.path.exists('filt_kmers_k%i.pickle' % options.k)):
                (kmers1, kmers2) = cPickle.load(open(('filt_kmers_k%i.pickle' % options.k), 'r'))
            else:
                kmers1, kmers2 = prepare_kmers(options, exon_t_gene)
                kmers1, kmers2 = clean_kmers(options, kmers1, kmers2)

            fastq_list = glob.glob(os.path.join(options.dir_fastq, '*.fastq')) \
                         + glob.glob(os.path.join(options.dir_fastq, '*.fastq.gz')) \
                         + glob.glob(os.path.join(options.dir_fastq, '*.fq')) \
                         + glob.glob(os.path.join(options.dir_fastq, '*.fq.gz'))
            if options.separate_files:
                header = fastq_list
            else:
                header = ','.join(fastq_list)

            data = get_counts_from_multiple_fastq(fastq_list, kmers1, kmers2, options)

        elif options.dir_bam != '-':
            if options.sparse_bam:
                bam_list = glob.glob(os.path.join(options.dir_bam, '*.hdf5'))
                header = bam_list  ### change this TODO
                data = get_counts_from_multiple_bam_sparse(bam_list, exon_t_gene)  ### REMOVE
            else:
                bam_list = glob.glob(os.path.join(options.dir_bam, '*.bam'))
                header = bam_list  ### change this TODO
                data = get_counts_from_multiple_bam(bam_list, exon_t_gene)  ### REMOVE
        elif options.fn_bam != '-':
            header = [options.fn_bam]  ### change this TODO
            if options.sparse_bam:
                data = get_counts_from_multiple_bam_sparse([options.fn_bam], exon_t_gene)  ### REMOVE
            else:
                data = get_counts_from_multiple_bam([options.fn_bam], exon_t_gene)  ### REMOVE


        ### normalize counts by exon length
        logging.info("Normalize counts by exon length")
        if options.qmode == 'raw':
            exonl = sp.array([int(x.split(':')[1].split('-')[1])
                            - int(x.split(':')[1].split('-')[0]) + 1 for x in exon_t_gene[:, :2].ravel('C')],
                            dtype='float') / 1000.
            data /= sp.tile(exonl[:, sp.newaxis], data.shape[1])

        ### Calculate 3'/5' Bias
        logging.info("Get counts from marginal exons")
        mycounts = __get_counts_from_marginal_exons(exon_t_gene, data)

    else: #MM options.dir_cnt != '-'
        count_files = 0
	cnt_file = None
        for cnt_file in glob.glob(options.fn_out + "counts*.tsv"):
            count_files = count_files + 1
        
	if(cnt_file != None and count_files > 0):
            cnt_list = glob.glob(os.path.join(options.dir_cnt, "*counts*.tsv"))
            header = cnt_list
            exon_t_gene = sp.loadtxt(cnt_file, delimiter='\t', dtype='string')[:, :-2]
            mycounts = sp.zeros((exon_t_gene.shape[0], count_files, 2))
            for i in xrange(count_files):
                mycounts[:, i, :] = sp.loadtxt(options.fn_out + 'counts' + str(i) + '.tsv', delimiter='\t', dtype='string')[:, -2:]
        else:
            print "No count files found in specified directory"
            sys.exit(2)
            


    if options.count_only:
        print "WARNING: Running only exon counts"
	pdb.set_trace()
        #MM CAVEAT: Order of exon-positions and counts might be switched (strand! --> see fct to get counts)
        for i in xrange(mycounts.shape[1]):
            exon_table = np.column_stack((exon_t_gene[:, :], mycounts[:, i, :]))
            sp.savetxt(options.fn_out + 'counts' + str(i) + '.tsv', exon_table, delimiter='\t', fmt='%s')
        sys.exit(0)


    logging.info("Calculate Bias and Find Median")
    vals = []
    for i in xrange(mycounts.shape[1]):
        if options.doPseudo:
            #AK: I had to filter for counts with some signal, otherwise the median is always 1.0
            i_ok = ((mycounts[:, i, 1] > 0) | (mycounts[:, i, 0] > 0))
            ratio = sp.percentile((mycounts[i_ok, i, 1] + 1) / (mycounts[i_ok, i, 0] + 1), 50)
        else:
            i_ok = ((mycounts[:, i, 1] > 0) & (mycounts[:, i, 0] > 0))
            ratio = sp.percentile(mycounts[i_ok, i, 1] / mycounts[i_ok, i, 0], 50)

        assert sp.sum(sp.isnan(ratio)) + sp.sum(sp.isinf(ratio)) == 0
        vals.append(ratio)

    vals = sp.array(vals)

    iqr = ((sp.percentile(vals, 75) - sp.percentile(vals, 25)) * 1.5)

    logging.info("Tukey Filter is estimated to be %f" % (iqr + sp.percentile(vals, 75)))
    if len(vals) > 1:
        print "Tukey Filter is estimated to be %f" % (iqr + sp.percentile(vals, 75))
        print "Tukey Filter is estimated to be %f" % (sp.percentile(vals, 25) - iqr)

    sp.savetxt('%s_sample_a_ratio_%s.tsv' % (options.fn_out, options.length),
               sp.vstack((header, vals.astype('string'))).T, delimiter='\t', fmt='%s')

    if options.doPlot:
        logging.info("Plot all samples")

        baseline_data = sp.loadtxt(options.fn_sample_ratio, delimiter='\t', dtype='string')
        baseline_data = baseline_data[:, 1].astype('float')

        base_p_val = sp.hstack((baseline_data, vals))
        midx = sp.hstack((sp.ones(baseline_data.shape[0]), sp.zeros(vals.shape[0]))).astype('bool')
        plotBias(base_p_val, '%s_bias_sorted_vline_%s.png' % (options.fn_out, options.length), midx)
        midx = sp.hstack((sp.ones(baseline_data.shape[0]), sp.zeros(vals.shape[0]))).astype('bool')
        plotBias(base_p_val, '%s_bias_sorted_vline_log_%s.png' % (options.fn_out, options.length), midx, logScale=True)


if __name__ == "__main__":
    main()
