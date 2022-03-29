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


def parse_options(argv):
    parser = OptionParser()

    sample_input = OptionGroup(parser, 'Input')

    sample_input.add_option('', '--bam_dir',           dest='dir_bam', metavar='FILE', help='Directory of bam files', default='-')
    sample_input.add_option('', '--fastq_dir',         dest='dir_fastq', metavar='FILE', help='Directory of fastq files', default='-')
    sample_input.add_option('', '--bam_fn',            dest='fn_bam', metavar='FIlE', help='Specifies single bam file', default='-')
    sample_input.add_option('', '--cnt_dir',           dest='dir_cnt', metavar='FILE', help='Directory of pre-produced tab delimited count files', default='-')

    sample_input.add_option('', '--anno_fn',           dest='fn_anno', metavar='FILE', help='Annotation', default='-')

    sample_input.add_option('', '--genome',            dest='fn_genome', metavar='FILE', help='Path to genome file in fasta', default='-')
    sample_input.add_option('', '--gene_list',         dest='fn_genes', metavar='FILE', help='File with gene-names to use', default='-')
    sample_input.add_option('', '--separate_files_ON', dest='separateFiles', action="store_true", help='Consider all input files individually [off]', default=False)
    sample_input.add_option('', '--sparse_bam_ON',     dest='sparseBam', action="store_true", help='Input BAM files are in sparse hdf5 format [off]', default=False)

    sample_output = OptionGroup(parser, 'Output')

    sample_output.add_option('', '--out_fn',      dest='fn_out', metavar='FILE', help='prefix for output', default='out')
    sample_output.add_option('', '--anno_tmp_fn', dest='fn_anno_tmp', metavar='FILE', help='Temp file for storing anno info', default='anno.tmp')
    sample_output.add_option('', '--pickle_all',  dest='fn_pickle_all', metavar='FILE', help='Pickle file for storing all kmers', default=None)
    sample_output.add_option('', '--pickle_filt', dest='fn_pickle_filt', metavar='FILE', help='Pickle file for storing filtered/cleaned kmers', default=None)

    opt_gen = OptionGroup(parser, 'General Options')

    opt_gen.add_option('', '--quant',                     dest='qMode', metavar='STRING', help='What type of quantification to use [rpkm,raw]', default='raw')
    opt_gen.add_option('', '--length',                    dest='readLength', metavar='STRING', help='Length filter [uq,mq,lq]', default='uq')

    opt_gen.add_option('', '--log',                       dest='fn_log', metavar='FILE', help='Log file', default='out.log')
    opt_gen.add_option('', '--verbose_ON',                dest='isVerbose', action="store_true", help='Set Logger To Verbose', default=False)

    opt_gen.add_option('', '--mask_filter',               dest='filt', help='Mask all read-counts below this integer', default='0')
    opt_gen.add_option('', '--protein_coding_filter_OFF', dest="proteinCodingFilter", action="store_false", help="Consider only genes that are protein-coding", default=True)
    opt_gen.add_option('', '--length_filter_OFF',         dest="lengthFilter", action="store_false", help="Only consider genes of certain length", default=True)

    opt_kmer = OptionGroup(parser, 'Options for k-mer counting')

    opt_kmer.add_option('', '--kmer_length', dest='k', type='int', help='Length of k-mer for alignmentfree counting', default=27)
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
              + int(options.dir_cnt != '-') + int(options.dir_fastq != '-')) != 1:
        print("Please specify exactly one type of input file(s) (e.g.: Exon quantification, Bam Files)")
        parser.print_help()
        sys.exit(2)
    if options.dir_fastq != '-' and options.fn_genome == '-':
        print('For usage on fastq files a genome file in fasta needs to be provided via -G/--genome', file=sys.stderr)
        sys.exit(2)
    return options


def __get_counts_from_marginal_exons(exon_t_gene, data):
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


def main():
    # Parse options
    options = parse_options(sys.argv)

    # Set up logger
    logging.basicConfig(filename=options.fn_log, level=0, format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger()
    if options.isVerbose:
        console_handler = logging.StreamHandler()
        log.addHandler(console_handler)

    # Read annotation from file
    logging.info("Reading Annotation from file")

    exon_t_gene = get_annotation_table(
            options.fn_genes,
            options.fn_anno_tmp,
            options.fn_anno,
            options.proteinCodingFilter,
            options.lengthFilter,
            options.readLength)

    if options.dir_fastq != '-':
        if options.fn_pickle_filt is not None and os.path.exists(options.fn_pickle_filt):
            (kmers1, kmers2) = pickle.load(open(options.fn_pickle_filt, 'r'))
        elif os.path.exists('filt_kmers_k%i.pickle' % options.k):
            (kmers1, kmers2) = pickle.load(open(('filt_kmers_k%i.pickle' % options.k), 'r'))
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

        fastq_list = glob.glob(os.path.join(options.dir_fastq, '*.fastq')) \
            + glob.glob(os.path.join(options.dir_fastq, '*.fastq.gz')) \
            + glob.glob(os.path.join(options.dir_fastq, '*.fq')) \
            + glob.glob(os.path.join(options.dir_fastq, '*.fq.gz'))
        if options.separateFiles:
            header = fastq_list
        else:
            header = ','.join(fastq_list)

        data = get_counts_from_multiple_fastq(fastq_list, kmers1, kmers2, options)

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

    # normalize counts by exon length
    logging.info("Normalize counts by exon length")
    if options.qMode == 'raw':
        exonl = np.array([int(x.split(':')[1].split('-')[1])
                          - int(x.split(':')[1].split('-')[0]) + 1 for x in exon_t_gene[:, :2].ravel('C')],
                         dtype='float') / 1000.
        data /= np.tile(exonl[:, np.newaxis], data.shape[1])

    # Get counts from the first an last exon
    logging.info("Get counts from marginal exons")
    my_counts = __get_counts_from_marginal_exons(exon_t_gene, data)

    #MM CAVEAT: Order of exon-positions and counts might be switched (strand! --> see fct to get counts)
    np.savetxt(options.fn_out + "_header.tsv", header, delimiter="\t", fmt="%s")
    for i in range(my_counts.shape[1]):
        exon_table = np.column_stack((exon_t_gene[:, :], my_counts[:, i, :]))
        np.save(options.fn_out + '_counts_' + str(i) + '.npy', exon_table)
        np.savetxt(options.fn_out + '_counts_' + str(i) + '.tsv', exon_table, delimiter='\t', fmt='%s')


if __name__ == "__main__":
    main()
