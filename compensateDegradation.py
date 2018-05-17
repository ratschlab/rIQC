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

    opt_kmer.add_option('', '', dest='k', type='int', help='Length of k-mer for alignmentfree counting [27]',
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
    exonTgene = getAnnotationTable(options, False) #no filtering for certain gene length (lengthFilter=False)

    if options.fastq_dir != '-':
        if(options.fn_pickle_filt != None and os.path.exists(options.fn_pickle_filt)):
            (kmers1, kmers2) = cPickle.load(open(options.fn_pickle_filt, 'r'))
        elif(os.path.exists('filt_kmers_k%i.pickle' % options.k)):
            (kmers1, kmers2) = cPickle.load(open(('filt_kmers_k%i.pickle' % options.k), 'r'))
        else:
            kmers1, kmers2 = prepare_kmers(options, exonTgene)
            kmers1, kmers2 = clean_kmers(options, kmers1, kmers2)

        fastq_list = glob.glob(os.path.join(options.fastq_dir, '*.fastq')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fastq.gz')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fq')) \
                     + glob.glob(os.path.join(options.fastq_dir, '*.fq.gz'))

        data = get_counts_from_multiple_fastq(fastq_list, kmers1, kmers2, options)
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

        istart = exonpos == rec[0]
        iend = exonpos == rec[1]
        if sp.sum(istart) == 0:
            continue
        if sp.sum(iend) == 0:
           continue
        if exonpos[istart][0].split(':')[-1] == '-' \
                and int(exonpos[istart][0].split(':')[1].split('-')[0]) < int(exonpos[iend][0].split(':')[1].split('-')[0]):
            istart, iend = iend, istart

        mycounts[i, :, 0] = data[istart, :]
        mycounts[i, :, 1] = data[iend, :]

    exonLengths = exonTgene[:, 4].astype(float)
    upperLengthBound = np.ceil(np.amax(exonLengths))
    gradients = sp.zeros((mycounts.shape[0], 2))

    assert exonLengths.shape[0] == gradients.shape[0]

    for i in range(gradients.shape[1]):
        assert gradients.shape[1] == 2
        iOK = ((mycounts[:, i, 1] > 0))
        #set all zero-counts from the 5' end to 1 to avoid division by zero
        mycounts[np.where(mycounts[:, i, 0] == 0)[0], i, 0] = 1

        # gradients-entrys stay 0 if they are not in iOK
        gradients[iOK, i] = (mycounts[iOK, i, 1] - mycounts[iOK, i, 0]) / exonLengths[iOK]

    # average gradients with #genes that contribute
    avg_grad_1 = sp.zeros(options.nmb_bins, 2)
    avg_grad_2 = sp.zeros(options.nmb_bins, 2)


    idx_d_1 = np.where(gradients[:, 0] > 0) #we get back indices
    idx_d_2 = np.where(gradients[:, 1] > 0)
    for i in range(options.nmb_bins):
        idx_l = np.where(exonLengths > upperLengthBound / options.nmb_bins * i
                         and exonLengths <= upperLengthBound / options.nmb_bins * (i+1) )
        #for kmer1
        comb_idx = sp.intersect1d(idx_d_1, idx_l) #indices of genes that have right length and gradient is non-zero
        if(idx_d_1.shape[0] != 0):
            avg_grad_1[i, 0] = sp.sum(gradients[comb_idx, 0]) / idx_d_1.shape[0]
            avg_grad_1[i, 1] = idx_d_1.shape[0]
        else:
            avg_grad_1[i, 0] = 0
            avg_grad_1[i, 1] = idx_d_1.shape[0]

        #for kmer2
        comb_idx = sp.intersect1d(idx_d_2, idx_l)
        if (idx_d_2.shape[0] != 0):
            avg_grad_2[i, 0] = sp.sum(gradients[comb_idx, 0]) / idx_d_2.shape[0]
            avg_grad_2[i, 1] = idx_d_2.shape[0]
        else:
            avg_grad_2[i, 0] = 0
            avg_grad_2[i, 1] = idx_d_2.shape[0]

    logging.info("The average gradients per bin are: ")
    logging.info(avg_grad_1)
    logging.info(avg_grad_2)






if __name__ == "__main__":
    main()