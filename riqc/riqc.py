import scipy.stats as spst
import numpy as np

import sys
import os
import glob
import fnmatch
import pickle
from optparse import OptionParser, OptionGroup
import logging

from .libs.annotation import *
from .libs.viz import *
from .libs.bam import *
from .libs.bam_sparse import *
from .libs.kmer import *
from .scalingFactors import get_scaling_factors


def parse_options(argv):
    parser = OptionParser()

    sample_input = OptionGroup(parser, 'Input')

    sample_input.add_option('', '--bam_dir',          dest='dir_bam', metavar='FILE', help='Directory of bam files', default='-')
    sample_input.add_option('', '--bam_fn',           dest='fn_bam', metavar='FIlE', help='Specifies single bam file', default='-')
    sample_input.add_option('', '--fastq_dir',        dest='dir_fastq', metavar='FILE', help='Directory of fastq files', default='-')
    sample_input.add_option('', '--fastq_fn',        dest='fn_fastq', metavar='FILE', help='Specifies single fastq file', default='-')
    sample_input.add_option('', '--cnt_dir',          dest='dir_cnt', metavar='FILE', help='Directory of pre-produced tab delimited count files', default='-')

    sample_input.add_option('', '--scale_factors_dir', dest='dir_scale_factors', metavar='FILE', help='Directory of scaling factor files', default='-')

    sample_input.add_option('', '--anno_fn',           dest='fn_anno', metavar='FILE', help='Annotation', default='-')

    sample_input.add_option('', '--genome',            dest='fn_genome', metavar='FILE', help='Path to genome file in fasta', default='-')
    sample_input.add_option('', '--gene_list',         dest='fn_genes', metavar='FILE', help='File with gene-names to use', default='-')

    sample_input.add_option('', '--separate_files_ON', dest='separateFiles', action="store_true", help='Consider all input files individually [off]', default=False)
    sample_input.add_option('', '--sparse_bam_ON',     dest='sparseBam', action="store_true", help='Input BAM files are in sparse hdf5 format [off]', default=False)

    sample_output = OptionGroup(parser, 'Output')

    sample_output.add_option('', '--out_dir',     dest='dir_out', metavar='FILE', help='directory to store output in', default='.')
    sample_output.add_option('', '--out_fn',      dest='fn_out', metavar='FILE', help='prefix for output', default='out')
    sample_output.add_option('', '--anno_tmp_fn', dest='fn_anno_tmp', metavar='FILE', help='Temp file for storing anno info', default='')
    sample_output.add_option('', '--pickle_all',  dest='fn_pickle_all', metavar='FILE', help='Pickle file for storing all kmers', default=None)
    sample_output.add_option('', '--pickle_filt', dest='fn_pickle_filt', metavar='FILE', help='Pickle file for storing filtered/cleaned kmers', default=None)

    opt_gen = OptionGroup(parser, 'General Options')
    
    opt_gen.add_option('', '--quant',           dest='qMode', metavar='STRING', help='What type of quantification to use [rpkm,raw]', default='raw')
    opt_gen.add_option('', '--pseudo_count_ON', dest='doPseudo', action="store_true", help='Add Pseudocounts to ratio', default=False)
    opt_gen.add_option('', '--length',          dest='readLength', metavar='STRING', help='Length filter [uq,mq,lq]', default='uq')

    opt_gen.add_option('', '--score_on_bases_ON',    dest='baseScore', action="store_true",
                       help="Calculate score not from last and first exon but from certain amount of bases at beginning and end that can be set via "
                            "--base_number (normalized length)", default=False)
    opt_gen.add_option('', '--base_number',     dest="baseNumber", type='int',
                       help="Number of bases at beginning/end for calculating score (only relevant if --score_on_bases_ON is set)", default=100)

    opt_gen.add_option('', '--log',             dest='fn_log', metavar='FILE', help='Log file', default='out.log')
    opt_gen.add_option('', '--verbose_ON',      dest='isVerbose', action="store_true", help='Set Logger To Verbose', default=False)
    opt_gen.add_option('', '--plot_ON',         dest='doPlot', action="store_true", help='Plot figures', default=False)
    opt_gen.add_option('', '--fn_sample_ratio', dest='fn_sample_ratio', metavar='FILE', help='Sample Ratios in relation to yours', default=os.path.join(os.path.realpath(__file__).rsplit('/', 1)[:-1][0], 'data', 'sampleRatios/TCGA_sample_a_ratio_uq.tsv'))

    opt_gen.add_option('', '--mask_filter',               dest='filt', help='Mask all read-counts below this integer', default='0')
    opt_gen.add_option('', '--protein_coding_filter_OFF', dest="proteinCodingFilter", action="store_false", help="Consider only genes that are protein-coding", default=True)
    opt_gen.add_option('', '--length_filter_OFF',         dest="lengthFilter", action="store_false", help="Only consider genes of certain length", default=True)

    opt_gen.add_option('', '--save_counts_ON',  dest='saveCounts', action='store_true', help='Store the exon counts in .npy and .tsv files for later use', default=False)
    opt_gen.add_option('', '--scale_counts_ON', dest='scaleCounts', action='store_true', help='Scale counts with pre-computed scaling factors for degradation compensation', default=False)

    opt_gen.add_option('', '--scale_mode', dest='scaleMode', metavar='STRING', help='How to scale counts [first, last, pseudoFirst]', default='first')

    opt_gen.add_option('', '--legacy', dest='legacy', action="store_true", help='Switch on some legacy behavior', default=False)
    
    opt_kmer = OptionGroup(parser, 'Options for k-mer counting')

    opt_kmer.add_option('', '--kmer_length', dest='k', type='int', help='Length of k-mer for alignmentfree counting [27]', default=27)
    opt_kmer.add_option('', '--reads_kmer',  dest='kmerThresh', type='float', help='Required active reads per sample or if in [0, 1] then fraction of input reads considered', default=50000)
    opt_kmer.add_option('', '--step_k',      dest='step_k', type='int', help='Step-size for k-mer counting', default=4)

    parser.add_option_group(sample_input)
    parser.add_option_group(sample_output)
    parser.add_option_group(opt_gen)
    parser.add_option_group(opt_kmer)
    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)
    if np.sum(int(options.dir_bam != '-') + int(options.fn_bam != '-')
              + int(options.dir_cnt != '-') + int(options.dir_fastq != '-') + int(options.fn_fastq != '-')) != 1:
        print("Please specify exactly one type of input file(s) (e.g.: Bam, Fastq, Count file)")
        parser.print_help()
        sys.exit(2)
    if (options.dir_fastq != '-' or options.fn_fastq != '-') and (options.fn_genome == '-'):
        print('For usage on fastq files a genome file in fasta needs to be provided via -G/--genome', file=sys.stderr)
        sys.exit(2)
    if options.dir_scale_factors == '-' and options.scaleCounts:
        print('CAVEAT: You have not provided pre-computed scaling factors ' \
              '--> they will be computed on the fly using default settings (see scalingFactors.py)')
    return options


def __get_counts_from_marginal_exons(exon_t_gene, data):

    my_counts = np.zeros((exon_t_gene.shape[0], data.shape[1], 2))

    for i, rec in enumerate(exon_t_gene):

        i_start = i*2
        i_end = i*2 + 1

        if rec[0].split(':')[-1] == '-' and \
                int(rec[0].split(':')[1].split('-')[0]) \
                < int(rec[1].split(':')[1].split('-')[0]):
            i_start, i_end = i_end, i_start

        my_counts[i, :, 0] = data[i_start, :]
        my_counts[i, :, 1] = data[i_end, :]

    return my_counts


def main():
    # Parse options
    options = parse_options(sys.argv)

    # Set up logger
    logging.basicConfig(filename=(options.dir_out + "/" + options.fn_log), level=0, format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger()
    if options.isVerbose:
        console_handler = logging.StreamHandler()
        log.addHandler(console_handler)

    if options.dir_cnt != '-':
        count_files = 0
        cnt_file = None
        for cnt_file in glob.glob(options.dir_cnt + "/counts_*.npy"):
            count_files = count_files + 1

        if cnt_file is not None and count_files > 0:
            header = np.loadtxt(options.dir_cnt + "/counts_header.tsv", delimiter="\t", dtype="str")
            exon_t_gene = np.load(cnt_file)[:, :-2]
            my_counts = np.zeros((exon_t_gene.shape[0], count_files, 2))
            for i in range(count_files):
                my_counts[:, i, :] = np.load(options.dir_cnt + '/counts_' + str(i) + '.npy')[:, -2:]
        else:
            print(("No count files found in specified directory %s" % options.dir_cnt))
            sys.exit(1)

    else:
        # Read annotation from file
        logging.info("Reading Annotation from file")
        
        exon_t_gene = get_annotation_table(
            options.fn_genes,
            options.fn_anno_tmp,
            options.fn_anno,
            options.proteinCodingFilter,
            options.lengthFilter,
            options.readLength,
            options.legacy)

        if options.dir_fastq != '-' or options.fn_fastq != '-':
            if options.fn_pickle_filt is not None \
                    and os.path.exists(options.fn_pickle_filt):
                (kmers1, kmers2) = pickle.load(open(options.fn_pickle_filt, 'rb'), encoding='latin1')
            elif os.path.exists('filt_kmers_k%i.pickle' % options.k):
                (kmers1, kmers2) = pickle.load(open(('filt_kmers_k%i.pickle' % options.k), 'rb'), encoding='latin1')
            else:
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

            if options.fn_fastq != '-':
                fastq_list = [options.fn_fastq]
            elif options.dir_fastq != '-':
                fastq_list = glob.glob(os.path.join(options.dir_fastq, '*.fastq')) \
                    + glob.glob(os.path.join(options.dir_fastq, '*.fastq.gz')) \
                    + glob.glob(os.path.join(options.dir_fastq, '*.fq')) \
                    + glob.glob(os.path.join(options.dir_fastq, '*.fq.gz'))
            if options.separateFiles:
                header = fastq_list
            else:
                header = ','.join(fastq_list)

            data = get_counts_from_multiple_fastq(fastq_list, kmers1, kmers2,
                                                  options.separateFiles,
                                                  options.kmerThresh,
                                                  options.k,
                                                  options.step_k)

        elif options.dir_bam != '-':
            if options.sparseBam:
                bam_list = glob.glob(os.path.join(options.dir_bam, '*.hdf5'))
                header = bam_list  ### change this TODO
                data = get_counts_from_multiple_bam_sparse(bam_list, exon_t_gene)  ### REMOVE
            else:
                bam_list = glob.glob(os.path.join(options.dir_bam, '*.bam'))
                header = bam_list  ### change this TODO
                data = get_counts_from_multiple_bam(bam_list, exon_t_gene)  ### REMOVE
        elif options.fn_bam != '-':
            header = [options.fn_bam]  ### change this TODO
            if options.sparseBam:
                data = get_counts_from_multiple_bam_sparse([options.fn_bam], exon_t_gene)  ### REMOVE
            else:
                data = get_counts_from_multiple_bam([options.fn_bam], exon_t_gene)  ### REMOVE


        ### normalize counts by exon length
        logging.info("Normalize counts by exon length")
        if options.qMode == 'raw':
            exonl = np.array([int(x.split(':')[1].split('-')[1])
                             - int(x.split(':')[1].split('-')[0]) + 1 for x in exon_t_gene[:, :2].ravel('C')],
                             dtype='float') / 1000.
            data /= np.tile(exonl[:, np.newaxis], data.shape[1])

        # Get counts from first and last exon
        logging.info("Get counts from marginal exons")
        my_counts = __get_counts_from_marginal_exons(exon_t_gene, data)

        if options.saveCounts:
            # MM CAVEAT: Order of exon-positions and counts might be switched (strand! --> see fct to get counts)
            np.savetxt(options.dir_out + "/counts_header.tsv", header, delimiter="\t", fmt="%s")
            for i in range(my_counts.shape[1]):
                exon_table = np.column_stack((exon_t_gene[:, :], my_counts[:, i, :]))
                np.save(options.dir_out + '/counts_' + str(i) + '.npy', exon_table)
                np.savetxt(options.dir_out + '/counts_' + str(i) + '.tsv', exon_table, delimiter='\t', fmt='%s')

    if options.scaleCounts:
        if not (os.path.exists(options.dir_scale_factors + "/scalingFactors.npy")):
            description, scaling_factors = get_scaling_factors(exon_t_gene, my_counts)
        else:
            scaling_factors = np.load(options.dir_scale_factors + "/scalingFactors.npy")
        exon_lengths = exon_t_gene[:, 4].astype(float)
        for i in range(my_counts.shape[1]):
            # MM j corresponds to number of bins
            for j in range(scaling_factors.shape[1]):
                low_b = scaling_factors[i, j, 2]
                up_b = scaling_factors[i, j, 3]
                factor = scaling_factors[i, j, 0]
                i_ok = np.intersect1d(np.where(low_b < exon_lengths)[0], np.where(exon_lengths <= up_b)[0])

                if options.scaleMode == 'first':
                    my_counts[i_ok, i, 0] = my_counts[i_ok, i, 0] * factor
                elif options.scaleMode == 'last':
                    if factor != 0:
                        my_counts[i_ok, i, 1] = my_counts[i_ok, i, 1] / factor
                else:
                    assert options.scaleMode == 'pseudoFirst'
                    i_ok_scale = np.intersect1d(i_ok, np.where(my_counts[:, i, 1] > 0)[0])
                    my_counts[i_ok, i, :] = my_counts[i_ok, i, :] + 1
                    my_counts[i_ok_scale, i, 0] = my_counts[i_ok_scale, i, 0] * factor

        # MM Save scaled counts for experimental purposes - can be removed later
        np.savetxt(options.dir_out + "/scaledCounts_header.tsv", header, delimiter="\t", fmt="%s")
        for i in range(my_counts.shape[1]):
            exon_table = np.column_stack((exon_t_gene[:, :], my_counts[:, i, :]))
            np.save(options.dir_out + '/scaledCounts_' + str(i) + '.npy', exon_table)
            np.savetxt(options.dir_out + '/scaledCounts_' + str(i) + '.tsv', exon_table, delimiter='\t', fmt='%s')

    logging.info("Calculate Bias and Find Median")

    vals = []
    for i in range(my_counts.shape[1]):
        if options.doPseudo:
            #AK: I had to filter for counts with some signal, otherwise the median is always 1.0
            i_ok = ((my_counts[:, i, 1] > 0) | (my_counts[:, i, 0] > 0))
            ratio = np.percentile((my_counts[i_ok, i, 1] + 1) / (my_counts[i_ok, i, 0] + 1), 50)
        else:
            i_ok = ((my_counts[:, i, 1] > 0) & (my_counts[:, i, 0] > 0))
            ratio = np.percentile(my_counts[i_ok, i, 1] / my_counts[i_ok, i, 0], 50)

        assert np.sum(np.isnan(ratio)) + np.sum(np.isinf(ratio)) == 0
        vals.append(ratio)

    vals = np.array(vals)

    iqr = ((np.percentile(vals, 75) - np.percentile(vals, 25)) * 1.5)

    logging.info("Tukey Filter is estimated to be %f" % (iqr + np.percentile(vals, 75)))
    if len(vals) > 1:
        print(("Tukey Filter is estimated to be %f" % (iqr + np.percentile(vals, 75))))
        print(("Tukey Filter is estimated to be %f" % (np.percentile(vals, 25) - iqr)))

    np.savetxt(options.dir_out + '/%s_sample_a_ratio_%s.tsv' % (options.fn_out, options.readLength),
               np.vstack((header, vals.astype('str'))).T, delimiter='\t', fmt='%s')

    if options.doPlot:
        logging.info("Plot all samples")

        baseline_data = np.loadtxt(options.fn_sample_ratio, delimiter='\t', dtype='str')
        baseline_data = baseline_data[:, 1].astype('float')

        base_p_val = np.hstack((baseline_data, vals))
        midx = np.hstack((np.ones(baseline_data.shape[0]), np.zeros(vals.shape[0]))).astype('bool')
        plotBias(base_p_val, '%s_bias_sorted_vline_%s.png' % (options.fn_out, options.readLength), midx)
        midx = np.hstack((np.ones(baseline_data.shape[0]), np.zeros(vals.shape[0]))).astype('bool')
        plotBias(base_p_val, '%s_bias_sorted_vline_log_%s.png' % (options.fn_out, options.readLength), midx, logScale=True)

    return 0

if __name__ == "__main__":
    sys.exit(main())

