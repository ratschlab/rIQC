import numpy as np
import pdb
import scipy as sp
from optparse import OptionParser, OptionGroup
import logging
import glob

from libs.annotation import *
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
    opt_gen.add_option('', '--bins', dest='nmb_bins', type='int',
                       help='Number of bins for different gene lengths', default=10)

    opt_gen.add_option('', '--log', dest='fn_log', metavar='FILE',
                       help='Log file', default='out.log')
    opt_gen.add_option('', '--verbose', dest='isVerbose', action="store_true",
                       help='Set Logger To Verbose', default=False)
    opt_gen.add_option('', '--sparse_bam', dest='sparse_bam', action="store_true",
                       help='Input BAM files are in sparse hdf5 format [off]', default=False)

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
              + int(options.dir_cnt != '-') + int(options.fastq_dir != '-')) != 1:
        print "Please specify exactly one type of input file(s) (e.g.: Exon quantification, Bam Files)"
        parser.print_help()
        sys.exit(2)
    if options.fastq_dir != '-' and options.fn_genome == '-':
        print >> sys.stderr, 'For usage on fastq files a genome file in fasta needs to be provided via -G/--genome'
        sys.exit(2)
    return options


def __get_counts_of_marginal_exons(exon_t_gene, data):
    mycounts = sp.zeros((exon_t_gene.shape[0], data.shape[1], 2))

    for i, rec in enumerate(exon_t_gene):

        istart = i * 2
        iend = i * 2 + 1

        if rec[0].split(':')[-1] == '-' and \
                int(rec[0].split(':')[1].split('-')[0]) \
                < int(rec[1].split(':')[1].split('-')[0]):
            istart, iend = iend, istart

        mycounts[i, :, 0] = data[istart, :]
        mycounts[i, :, 1] = data[iend, :]

    return mycounts


def main():
    # Parse options
    options = parse_options(sys.argv)

    # set up logger
    logging.basicConfig(filename=options.fn_log, level=0, format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger()
    if options.isVerbose:
        consoleHandler = logging.StreamHandler()
        log.addHandler(consoleHandler)

    # Read annotation from file
    logging.info("Reading Annotation from file")

    # MM type(exon_t_gene) is np.ndarray
    # MM rows:_genes
    # MM 6 columns: first_constitutive_exon last_constitutive_exon chromosome strand distance(basepairs) genname
    exon_t_gene = get_annotation_table(options, False) #no filtering for certain gene length (lengthFilter=False)

    if options.fastq_dir != '-':
        if(options.fn_pickle_filt != None and os.path.exists(options.fn_pickle_filt)):
            (kmers1, kmers2) = cPickle.load(open(options.fn_pickle_filt, 'r'))
        elif(os.path.exists('filt_kmers_k%i.pickle' % options.k)):
            (kmers1, kmers2) = cPickle.load(open(('filt_kmers_k%i.pickle' % options.k), 'r'))
        else:
            # MM type(kmers1) is list
            # MM kmers1 contains kmers in first consecutive exon
            # MM kmers2 contains kmers in last consecutive exon
            kmers1, kmers2 = prepare_kmers(options, exon_t_gene)
            kmers1, kmers2 = clean_kmers(options, kmers1, kmers2)

        fastq_list = glob.glob(os.path.join(options.fastq_dir, '*.fastq')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fastq.gz')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fq')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fq.gz'))

        # MM data is np.ndarray
        # MM data contains all counts for first and last exon of each gene
        # MM in case of multiple files in directory: each column contains counts for one file
        data = get_counts_from_multiple_fastq(fastq_list, kmers1, kmers2, options)

    elif options.dir_bam != '-':
        if options.sparse_bam:
            bam_list = glob.glob(os.path.join(options.dir_bam, '*.hdf5'))
            data = get_counts_from_multiple_bam_sparse(bam_list, exon_t_gene)
        else:
            bam_list = glob.glob(os.path.join(options.dir_bam, '*.bam'))
            data = get_counts_from_multiple_bam(bam_list, exon_t_gene)
    elif options.fn_bam != '-':
        if options.sparse_bam:
            data = get_counts_from_multiple_bam_sparse([options.fn_bam], exon_t_gene)
        else:
            data = get_counts_from_multiple_bam([options.fn_bam], exon_t_gene)


    ### normalize counts by exon length
    logging.info("Normalize counts by exon length")
    if (options.qmode == 'raw'):
        exonl = sp.array([int(x.split(':')[1].split('-')[1]) -
                          int(x.split(':')[1].split('-')[0]) + 1 for x in exon_t_gene[:, :2].ravel('C')],
                         dtype='float') / 1000.
        data /= sp.tile(exonl[:, sp.newaxis], data.shape[1])

    # get counts from both ends
    mycounts = __get_counts_of_marginal_exons(exon_t_gene, data)

    #MM for binning
    exonLengths = exon_t_gene[:, 4].astype(float)
    upperLengthBound = np.ceil(np.amax(exonLengths))

    scale = sp.zeros((mycounts.shape[0], mycounts.shape[1]))
    #MM average scales with interval and #genes that contribute
    avg_scale = sp.zeros((mycounts.shape[1], options.nmb_bins, 3))

    #MM for every file that was read in
    for i in xrange(mycounts.shape[1]):
        #MM only take genes that have some count
        #MM iOK is a boolean np.ndarray as long as number of genes
        iOK = sp.union1d(np.where(mycounts[:, i, 0] > 0)[0],
                         np.where(mycounts[:, i, 1] > 0)[0])

        #MM replace 0 with 1 to avoid division by 0
        i_avoid0 = np.where(mycounts[:, i, 0] == 0)[0]
        mycounts[i_avoid0, i, 0] = 1
        mycounts[i_avoid0, i, 1] = mycounts[i_avoid0, i, 1] + 1

        #MM scale-entrys stay 0 if they are not in iOK
        scale[iOK, i] = (mycounts[iOK, i, 1] / mycounts[iOK, i, 0])

        for j in range(options.nmb_bins):
            low_b = upperLengthBound / options.nmb_bins * j
            up_b = upperLengthBound / options.nmb_bins * (j + 1)
            idx_l = sp.intersect1d(np.where(low_b < exonLengths)[0], np.where(exonLengths <= up_b)[0])

            # indices of genes that have right length and a scale factor
            comb_idx = sp.intersect1d(iOK, idx_l)

            avg_scale[i, j, 1] = str(low_b) + "-" + str(up_b)
            if comb_idx.shape[0] != 0:
                avg_scale[i, j, 0] = sp.sum(scale[comb_idx, i]) / comb_idx.shape[0]
                avg_scale[i, j, 2] = comb_idx.shape[0]
            else:
                avg_scale[i, j, 0] = 0
                avg_scale[i, j, 2] = comb_idx.shape[0]

        sp.savetxt(options.fn_out + "_scalingFactors_" + str(i) + ".tsv", avg_scale[i, :, :], delimiter="\t", fmt="%s")


if __name__ == "__main__":
    main()
