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

    sampleinput.add_option('', '--bam_dir', dest='dir_bam', metavar='FILE', help='Directory of bam files',
                           default='-')
    sampleinput.add_option('', '--fastq_dir', dest='fastq_dir', metavar='FILE', help='Directory of fastq files',
                           default='-')
    sampleinput.add_option('', '--fn_bam', dest='fn_bam', metavar='FIlE', help='Specifies single bam file',
                           default='-')

    sampleinput.add_option('', '--fn_anno', dest='fn_anno', metavar='FILE', help='Annotation', default='-')

    sampleinput.add_option('', '--genome', dest='fn_genome', metavar='FILE', help='Path to genome file in fasta',
                           default='-')


    sampleinput.add_option('-i', '--genelist', dest='fn_genes', metavar='FILE', help='file with genenames to use',
                           default='-')
    sampleinput.add_option('', '--separate_files', dest='separate_files', action="store_true",
                           help='Consider all input files individually [off]', default=False)


    sampleoutput = OptionGroup(parser, 'Output')

    sampleoutput.add_option('', '--fn_out', dest='fn_out', metavar='FILE',
                            help='prefix for output', default='out')
    sampleoutput.add_option('', '--fn_anno_tmp', dest='fn_anno_tmp', metavar='FILE',
                            help='Temp file for storing anno info', default='anno.tmp')
    sampleoutput.add_option('', '--pickle_all', dest='fn_pickle_all', metavar='FILE',
                            help='Pickle file for storing all kmers', default=None)
    sampleoutput.add_option('', '--pickle_filt', dest='fn_pickle_filt', metavar='FILE',
                            help='Pickle file for storing filtered/cleaned kmers', default=None)


    opt_gen = OptionGroup(parser, 'General Options')

    opt_gen.add_option('', '--quant', dest='qmode', metavar='STRING',
                       help='What type of quantification to use [rpkm,raw]', default='raw')

    opt_gen.add_option('', '--count_only', dest='count_only', action="store_true",
                       help='Only do counting on given input [off]', default=False)

    opt_gen.add_option('', '--bins', dest='nmb_bins', type='int',
                       help='Number of bins for different gene lengths', default=10)

    opt_gen.add_option('', '--log', dest='fn_log', metavar='FILE', help='Log file', default='out.log')
    opt_gen.add_option('', '--verbose', dest='isVerbose', action="store_true", help='Set Logger To Verbose',
                       default=False)

    opt_gen.add_option('', '--sparse_bam', dest='sparse_bam', action="store_true",
                       help='Input BAM files are in sparse hdf5 format [off]', default=False)

    opt_gen.add_option('', '--protein-coding-filter_OFF', dest="protein_coding_filter", action="store_false",
                       help="Consider only genes that are protein-coding", default=True)


    opt_kmer = OptionGroup(parser, 'Options for k-mer counting')

    opt_kmer.add_option('', '--kmer_length', dest='k', type='int', help='Length of k-mer for alignmentfree counting [27]',
                        default=27)
    opt_kmer.add_option('', '--reads_kmer', dest='kmer_thresh', type='float',
                        help='Required active reads per sample [50000] / if btw 0 and 1 fraction of input reads considered',
                        default=50000)
    opt_kmer.add_option('', '--step_k', dest='step_k', type='int', help='Step-size for k-mer counting [4]', default=4)

    parser.add_option_group(sampleinput)
    parser.add_option_group(sampleoutput)
    parser.add_option_group(opt_gen)
    parser.add_option_group(opt_kmer)
    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)
    if sp.sum(int(options.dir_bam != '-') + int(options.fn_bam != '-') + int(options.fastq_dir != '-')) != 1:
        print "Please specify exactly one type of input file(s) (e.g.: Exon quantification, Bam Files)"
        parser.print_help()
        sys.exit(2)
    if options.fastq_dir != '-' and options.fn_genome == '-':
        print >> sys.stderr, 'For usage on fastq files a genome file in fasta needs to be provided via -G/--genome'
        sys.exit(2)
    return options



def main():
    # Parse options
    options = parse_options(sys.argv)

    #### set up logger
    logging.basicConfig(filename=options.fn_log, level=0, format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger()
    if options.isVerbose:
        consoleHandler = logging.StreamHandler()
        log.addHandler(consoleHandler)

    ### Read annotation from file
    logging.info("Reading Annotation from file")

    # MM type(exonTgene) is np.ndarray
    # MM rows:_genes
    # MM 6 columns: first_constitutive_exon last_constitutive_exon chromosome strand distance(basepairs) genname
    exonTgene = get_annotation_table(options, False) #no filtering for certain gene length (lengthFilter=False)

    if options.fastq_dir != '-':
        if(options.fn_pickle_filt != None and os.path.exists(options.fn_pickle_filt)):
            (kmers1, kmers2) = cPickle.load(open(options.fn_pickle_filt, 'r'))
        elif(os.path.exists('filt_kmers_k%i.pickle' % options.k)):
            (kmers1, kmers2) = cPickle.load(open(('filt_kmers_k%i.pickle' % options.k), 'r'))
        else:
            # MM type(kmers1) is list
            # MM kmers1 contains kmers in first consecutive exon
            # MM kmers2 contains kmers in last consecutive exon
            kmers1, kmers2 = prepare_kmers(options, exonTgene)
            kmers1, kmers2 = clean_kmers(options, kmers1, kmers2)

        fastq_list = glob.glob(os.path.join(options.fastq_dir, '*.fastq')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fastq.gz')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fq')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fq.gz'))

        # MM data is np.ndarray
        # MM data contains all counts for first and last exon of each gene
        # MM in case of multiple files in directory: each column contains counts for one file
        data = get_counts_from_multiple_fastq(fastq_list, kmers1, kmers2, options)
        # MM exonpos is np.ndarray
        # MM exonpos contains chr-position-strand information for first and last exon of each gene
        exonpos = exonTgene[:, :2].ravel('C')

    elif options.dir_bam != '-':
        if options.sparse_bam:
            bam_list = glob.glob(os.path.join(options.dir_bam, '*.hdf5'))
            data = get_counts_from_multiple_bam_sparse(bam_list, exonTgene)
        else:
            bam_list = glob.glob(os.path.join(options.dir_bam, '*.bam'))
            data = get_counts_from_multiple_bam(bam_list, exonTgene)
        exonpos = exonTgene[:, :2].ravel('C')
    elif options.fn_bam != '-':
        if options.count_only:
            print "WARNING: Running only gene counts"
            exonTable = sp.sort(exonTgene[:, [0, 1]].ravel())
            data = get_counts_from_single_bam(options.fn_bam, exonTable)
            sp.savetxt(options.fn_out + 'counts.tsv', sp.vstack((exonTable, data[::2])).T, delimiter='\t', fmt='%s')
            sys.exit(0)
        else:
            if options.sparse_bam:
                data = get_counts_from_multiple_bam_sparse([options.fn_bam], exonTgene)
            else:
                data = get_counts_from_multiple_bam([options.fn_bam], exonTgene)
            exonpos = exonTgene[:, :2].ravel('C')

    ### normalize counts by exon length
    logging.info("Normalize counts by exon length")
    if (options.qmode == 'raw'):
        exonl = sp.array([int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) + 1 for x in exonpos],
                         dtype='float') / 1000.
        data /= sp.tile(exonl[:, sp.newaxis], data.shape[1])

    # get counts from both ends
    mycounts = sp.zeros((exonTgene.shape[0], data.shape[1], 2))
    for i, rec in enumerate(exonTgene):

        ##MM istart is a boolean np.ndarray
        istart = exonpos == rec[0]
        iend = exonpos == rec[1]
        if sp.sum(istart) == 0:
            continue
        if sp.sum(iend) == 0:
           continue
        if exonpos[istart][0].split(':')[-1] == '-' \
                and int(exonpos[istart][0].split(':')[1].split('-')[0]) < int(exonpos[iend][0].split(':')[1].split('-')[0]):
            istart, iend = iend, istart

        # MM mycounts is np.ndarray
        mycounts[i, :, 0] = data[istart, :]
        mycounts[i, :, 1] = data[iend, :]

    #MM for binning
    exonLengths = exonTgene[:, 4].astype(float)
    upperLengthBound = np.ceil(np.amax(exonLengths))
    ##MM
    scale = sp.zeros((mycounts.shape[0], mycounts.shape[1]))
    #MM average scales with #genes that contribute
    avg_scale = sp.zeros((mycounts.shape[1], options.nmb_bins, 2))

    #MM for every file that was read in
    for i in xrange(mycounts.shape[1]):
        #MM avoid division by zero
        #MM iOK is a boolean np.ndarray as long as number of genes
        iOK = np.where(exonLengths > 0)[0]

        #MM scale-entrys stay 0 if they are not in iOK
        scale[iOK, i] = (mycounts[iOK, i, 1] / mycounts[iOK, i, 0])

        for j in range(options.nmb_bins):
            idx_l = sp.intersect1d(np.where(upperLengthBound / options.nmb_bins * j < exonLengths)[0],
                                   np.where(exonLengths <= upperLengthBound / options.nmb_bins * (j + 1))[0])

            # indices of genes that have right length and gradient is non-zero
            comb_idx = sp.intersect1d(iOK, idx_l)

            if(comb_idx.shape[0] != 0):
                avg_scale[i, j, 0] = sp.sum(scale[comb_idx, i]) / comb_idx.shape[0]
                avg_scale[i, j, 1] = comb_idx.shape[0]
            else:
                avg_scale[i, j, 0] = 0
                avg_scale[i, j, 1] = comb_idx.shape[0]











if __name__ == "__main__":
    main()
