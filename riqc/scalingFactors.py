import numpy as np
from optparse import OptionParser, OptionGroup
import logging
import glob

from .libs.annotation import *
from .libs.bam import *
from .libs.bam_sparse import *
from .libs.kmer import *


def parse_options(argv):
    parser = OptionParser()

    sample_input = OptionGroup(parser, 'Input')

    sample_input.add_option('', '--bam_dir',   dest='dir_bam', metavar='FILE', help='Directory of bam files', default='-')
    sample_input.add_option('', '--fastq_dir', dest='dir_fastq', metavar='FILE', help='Directory of fastq files', default='-')
    sample_input.add_option('', '--bam_fn',    dest='fn_bam', metavar='FIlE', help='Specifies single bam file', default='-')
    sample_input.add_option('', '--cnt_dir',   dest='dir_cnt', metavar='FILE', help='Directory of pre-produced tab delimited count files', default='-')

    sample_input.add_option('', '--anno_fn',   dest='fn_anno', metavar='FILE', help='Annotation', default='-')

    sample_input.add_option('', '--genome',    dest='fn_genome', metavar='FILE', help='Path to genome file in fasta', default='-')
    sample_input.add_option('', '--gene_list', dest='fn_genes', metavar='FILE', help='file with gene-names to use', default='-')
    sample_input.add_option('', '--separate_files_ON', dest='separateFiles', action="store_true", help='Consider all input files individually [off]', default=False)
    sample_input.add_option('', '--sparse_bam_ON',     dest='sparseBam', action="store_true", help='Input BAM files are in sparse hdf5 format [off]', default=False)

    sample_output = OptionGroup(parser, 'Output')

    sample_output.add_option('', '--out_dir', dest='dir_out', metavar='FILE', help='directory to store output in', default='.')
    sample_output.add_option('', '--out_fn',      dest='fn_out', metavar='FILE', help='prefix for output', default='out')
    sample_output.add_option('', '--anno_tmp_fn', dest='fn_anno_tmp', metavar='FILE', help='Temp file for storing anno info', default=''
                                                                                                                                      '')
    sample_output.add_option('', '--pickle_all',  dest='fn_pickle_all', metavar='FILE', help='Pickle file for storing all kmers', default=None)
    sample_output.add_option('', '--pickle_filt', dest='fn_pickle_filt', metavar='FILE', help='Pickle file for storing filtered/cleaned kmers', default=None)

    opt_gen = OptionGroup(parser, 'General Options')

    opt_gen.add_option('', '--quant',            dest='qMode', metavar='STRING', help='What type of quantification to use [rpkm,raw]', default='raw')
    opt_gen.add_option('', '--pseudo_count_OFF', dest='doPseudo', action="store_false", help='Add Pseudocounts to ratio', default=True)
    opt_gen.add_option('', '--length',           dest='readLength', metavar='STRING', help='Length filter [uq,mq,lq]', default='uq')

    opt_gen.add_option('', '--log',              dest='fn_log', metavar='FILE', help='Log file', default='out.log')
    opt_gen.add_option('', '--verbose_ON',       dest='isVerbose', action="store_true", help='Set Logger To Verbose', default=False)

    opt_gen.add_option('', '--protein_coding_filter_ON', dest="proteinCodingFilter", action="store_true", help="Consider only genes that are protein-coding", default=False)
    opt_gen.add_option('', '--length_filter_ON', dest="lengthFilter", action="store_true", help="Only consider genes of certain length", default=False)
    opt_gen.add_option('', '--bins', dest='nmb_bins', type='int', help='Number of bins for different gene lengths', default=10)
    opt_gen.add_option('', '--relative_binning_ON', dest='relativeBinning', action="store_true", help="Have relative (to number of genes) bin boundaries instead of absolute values", default=False)
    opt_gen.add_option('', '--average_factors_ON', dest='averageFactors', action='store_true', help="Compute factors by using average (instead of median)", default=False)

    opt_gen.add_option('', '--legacy', dest='legacy', action="store_true", help='Switch on some legacy behavior', default=False)

    opt_kmer = OptionGroup(parser, 'Options for k-mer counting')

    opt_kmer.add_option('', '--kmer_length', dest='k', type='int', help='Length of k-mer for alignmentfree counting', default=27)
    opt_kmer.add_option('', '--reads_kmer',  dest='kmerThresh', type='float', help='Required active reads per sample or if in [0, 1] then fraction of input reads considered', default=50000)
    opt_kmer.add_option('', '--step_k', dest='step_k', type='int', help='Step-size for k-mer counting', default=4)

    parser.add_option_group(sample_input)
    parser.add_option_group(sample_output)
    parser.add_option_group(opt_gen)
    parser.add_option_group(opt_kmer)
    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)
    if np.sum(int(options.dir_bam != '-') + int(options.fn_bam != '-')
              + int(options.dir_cnt != '-') + int(options.dir_fastq != '-')) != 1:
        print("Please specify exactly one type of input file(s) (e.g.: Bam, Fastq, Count file)")
        parser.print_help()
        sys.exit(2)
    if options.dir_fastq != '-' and options.fn_genome == '-':
        print('For usage on fastq files a genome file in fasta needs to be provided via -G/--genome', file=sys.stderr)
        sys.exit(2)
    return options


def __get_counts_of_marginal_exons(exon_t_gene, data):
    my_counts = np.zeros((exon_t_gene.shape[0], data.shape[1], 2))

    for i, rec in enumerate(exon_t_gene):

        i_start = i * 2
        i_end = i * 2 + 1

        if rec[0].split(':')[-1] == '-' and \
                int(rec[0].split(':')[1].split('-')[0]) \
                < int(rec[1].split(':')[1].split('-')[0]):
            i_start, i_end = i_end, i_start

        my_counts[i, :, 0] = data[i_start, :]
        my_counts[i, :, 1] = data[i_end, :]

    return my_counts


def get_scaling_factors(exon_t_gene, my_counts, nmb_bins=10,
                          doPseudo=True, relativeBinning=True, averageFactors=False):
    #MM for binning
    exon_lengths = exon_t_gene[:, 4].astype(float)
    upper_length_bound = np.ceil(np.amax(exon_lengths))
    nmb_exons = exon_lengths.shape[0]

    scale = np.zeros((my_counts.shape[0], my_counts.shape[1]))
    #MM average scales with interval and #genes that contribute
    avg_scale = np.zeros((my_counts.shape[1], nmb_bins, 4))

    #MM for every file that was read in
    for i in range(my_counts.shape[1]):
        #MM only take genes that have some count
        #MM i_ok is a boolean np.ndarray as long as number of genes
        if doPseudo:
            i_ok = np.union1d(np.where(my_counts[:, i, 0] > 0)[0],
                              np.where(my_counts[:, i, 1] > 0)[0])
            # MM add pseudocount to avoid division by 0
            my_counts[i_ok, i, :] = my_counts[i_ok, i, :] + 1
        else:
            i_ok = np.intersect1d(np.where(my_counts[:, i, 0] > 0)[0],
                                  np.where(my_counts[:, i, 1] > 0)[0])

        #MM scale-entrys stay 0 if they are not in i_ok
        scale[i_ok, i] = (my_counts[i_ok, i, 1] / my_counts[i_ok, i, 0])

        for j in range(nmb_bins):
            if relativeBinning:
                idx_s = np.argsort(exon_lengths)
                low_b = nmb_exons / nmb_bins * j
                up_b = nmb_exons / nmb_bins * (j + 1)
                idx_l = idx_s[low_b:up_b]

                avg_scale[i, j, 2] = exon_lengths[idx_s[low_b]]
                avg_scale[i, j, 3] = exon_lengths[idx_s[up_b]]
            else:
                low_b = upper_length_bound / nmb_bins * j
                up_b = upper_length_bound / nmb_bins * (j + 1)
                idx_l = np.intersect1d(np.where(low_b < exon_lengths)[0], np.where(exon_lengths <= up_b)[0])

                avg_scale[i, j, 2] = low_b
                avg_scale[i, j, 3] = up_b

            # indices of genes that have right length and a scale factor
            comb_idx = np.intersect1d(i_ok, idx_l)
            

            if comb_idx.shape[0] != 0 and averageFactors:
                avg_scale[i, j, 0] = np.sum(scale[comb_idx, i]) / comb_idx.shape[0]
            elif comb_idx.shape[0] != 0:  # Take median
                avg_scale[i, j, 0] = np.percentile(scale[comb_idx, i], 50)
            else:
                avg_scale[i, j, 0] = 0
                assert comb_idx.shape[0] == 0
            avg_scale[i, j, 1] = comb_idx.shape[0]

    description = np.array([['scaling_factor_for_genes_with_length', 'number_of_genes_with_length', 'length_lower_bound', 'length_upper_bound']])
    assert description.shape[1] == avg_scale.shape[2]

    return description, avg_scale


def main():
    # Parse options
    options = parse_options(sys.argv)

    # Set up logger
    logging.basicConfig(filename=options.dir_out + "/" + options.fn_log, level=0, format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger()
    if options.isVerbose:
        console_handler = logging.StreamHandler()
        log.addHandler(console_handler)

    if options.dir_cnt != '-':
        count_files = 0
        cnt_file = None
        print(glob.glob(options.dir_cnt + "/counts_*.npy"))
        for cnt_file in glob.glob(options.dir_cnt + "/counts_*.npy"):
            count_files = count_files + 1

        if cnt_file is not None and count_files > 0:
            header = np.loadtxt(options.dir_cnt + "/counts_header.tsv", delimiter="\t", dtype="string")
            exon_t_gene = np.load(cnt_file)[:, :-2]
            my_counts = np.zeros((exon_t_gene.shape[0], count_files, 2))
            for i in range(count_files):
                my_counts[:, i, :] = np.load(options.dir_cnt + '/counts_' + str(i) + '.npy')[:, -2:]
        else:
            print("No count files found in specified directory")
            sys.exit(1)

    else:
        # Read annotation from file
        logging.info("Reading Annotation from file")

        # MM type(exon_t_gene) is np.ndarray
        # MM rows:_genes
        # MM 6 columns: first_constitutive_exon last_constitutive_exon chromosome strand distance(basepairs) genname
        exon_t_gene = get_annotation_table(
            options.fn_genes,
            options.fn_anno_tmp,
            options.fn_anno,
            options.proteinCodingFilter,
            options.lengthFilter,
            options.readLength,
            options.legacy)

        if options.dir_fastq != '-':
            if options.fn_pickle_filt is not None \
                    and os.path.exists(options.fn_pickle_filt):
                (kmers1, kmers2) = cPickle.load(open(options.fn_pickle_filt, 'r'))
            elif os.path.exists('filt_kmers_k%i.pickle' % options.k):
                (kmers1, kmers2) = cPickle.load(open(('filt_kmers_k%i.pickle' % options.k), 'r'))
            else:
                # MM type(kmers1) is list
                # MM kmers1 contains kmers in first consecutive exon
                # MM kmers2 contains kmers in last consecutive exon
                kmers1, kmers2 = prepare_kmers(
                    exon_t_gene,
                    options.fn_genome,
                    options.k)
                kmers1, kmers2 = clean_kmers(
                    kmers1, kmers2,
                    options.fn_pickle_all,
                    options.fn_pickle_filt,
                    options.k,
                    options.fn_genome)

            fastq_list = glob.glob(os.path.join(options.dir_fastq, '*.fastq')) \
                + glob.glob(os.path.join(options.dir_fastq, '*.fastq.gz')) \
                + glob.glob(os.path.join(options.dir_fastq, '*.fq')) \
                + glob.glob(os.path.join(options.dir_fastq, '*.fq.gz'))

            if options.separateFiles:
                header = fastq_list
            else:
                header = ','.join(fastq_list)

            # MM data is np.ndarray
            # MM data contains all counts for first and last exon of each gene
            # MM in case of multiple files in directory: each column contains counts for one file
            data = get_counts_from_multiple_fastq(fastq_list, kmers1, kmers2, options)

        elif options.dir_bam != '-':
            if options.sparseBam:
                bam_list = glob.glob(os.path.join(options.dir_bam, '*.hdf5'))
                data = get_counts_from_multiple_bam_sparse(bam_list, exon_t_gene)
            else:
                bam_list = glob.glob(os.path.join(options.dir_bam, '*.bam'))
                data = get_counts_from_multiple_bam(bam_list, exon_t_gene)
            header = bam_list
        elif options.fn_bam != '-':
            header = [options.fn_bam]
            if options.sparseBam:
                data = get_counts_from_multiple_bam_sparse([options.fn_bam], exon_t_gene)
            else:
                data = get_counts_from_multiple_bam([options.fn_bam], exon_t_gene)


        ### normalize counts by exon length
        logging.info("Normalize counts by exon length")
        if options.qMode == 'raw':
            exonl = np.array([int(x.split(':')[1].split('-')[1]) -
                              int(x.split(':')[1].split('-')[0]) + 1 for x in exon_t_gene[:, :2].ravel('C')],
                             dtype='float') / 1000.
            data /= np.tile(exonl[:, np.newaxis], data.shape[1])

        # get counts from both ends
        my_counts = __get_counts_of_marginal_exons(exon_t_gene, data)

    description, avg_scale = get_scaling_factors(exon_t_gene, my_counts, options.nmb_bins,
                                                 options.doPseudo, options.relativeBinning, options.averageFactors)
    np.save(options.dir_out + "/scalingFactors.npy", avg_scale)
    np.savetxt(options.dir_out + "/scalingFactors_header.tsv", header, delimiter="\t", fmt="%s")
    for i in range(my_counts.shape[1]):
        np.savetxt(options.dir_out + "/scalingFactors_" + str(i) + ".tsv", np.concatenate((description, avg_scale[i, :, :])), delimiter="\t", fmt="%s")


if __name__ == "__main__":
    main()
