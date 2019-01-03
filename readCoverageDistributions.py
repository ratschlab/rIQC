import os
import sys
import glob
from optparse import OptionParser, OptionGroup
import pysam
import time
import numpy as np
import scipy as sp
from libs.usefulTools import *

# Some Numbers
NMB_CHR = 23

#  GTF column numbers
SEQ_NAME = 0    # name of chromosome or scaffold
SOURCE = 1      # name of program that generated this feature
FEATURE = 2     # feature type name (e.g. "gene", "transcript", "exon") / type of feature (term or accession from SOFA sequence ontology)
START = 3       # start position of feature (seq numbering starting at 1)
END = 4         # end position of feature (seq numbering starting at 1)
SCORE = 5       # a floating point value
STRAND = 6      # + (forward) or - (reverse)
FRAME = 7       # 0/1/2 : position in codon
ATTRIBUTE = 8   # semicolon-separated list of tag-value pairs


# Average count per exon (histogram) (depending on location); for both all exons and only constitutive ones
def avg_count_per_exon():
    """
    Plots one histogram for all genes: Average count of first, second, ..., last exon with variance as error metric
    (positions at the end have most likely less samples than at the beginning)
    NOT considering length of exons or length of gene, also different genes have last exon at different positions (different amounts of exons)
    """
    print "Hi"


# Expression distribution over normalized gene length (in different length-bins-> same as already used; constitutive and not)
# # Generate average coverage distribution across all genes within that bin after projecting all genes in bin onto one length
def distr_over_gene_lengths():
    print "Hi"


# Distribution of length of last-constitutive exons (joint and in different length bins)
def last_const_exon_length():
    print "Hi"


# Position of last constitutive exon in the list of all exons (based on normalized gene length)
def last_const_exon_pos():
    print "Hi"

###############################################################################################


def get_counts_from_marginal_exons(exon_t_gene, data):
    my_counts = sp.zeros((exon_t_gene.shape[0], data.shape[1], 2))

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


def get_overlap_genes(fn_anno):
    """
    Returns a list of gene names which are overlapping
    """

    data = []
    for l in open(fn_anno, 'r'):
        if l[SEQ_NAME] == '#':
            continue  # comments
        l_spl = l.strip('\n').split('\t')
        if l_spl[2].lower() != 'gene':
            continue
        tags = dict()
        for t in l_spl[ATTRIBUTE].strip(';').split(';'):
            tt = t.strip(' ').split(' ')
            tags[tt[0]] = tt[1].strip('"')
        data.append([tags['gene_id'], '%s:%s-%s' % (l_spl[SEQ_NAME], l_spl[START], l_spl[END])])

    data = sp.array(data)

    # fix positions
    pos = data[:, 1]
    pos = sp.array([x.split(':')[0] + '-' + x.split(':')[1] for x in pos])
    pos = sp.array([x.strip('chr') for x in pos])
    pos = sp.array([x.split('-') for x in pos])
    pos[pos[:, 0] == 'X', 0] = '23'
    pos[pos[:, 0] == 'Y', 0] = '24'

    # filter weird things like mitochondria etc.
    i_ok = np.core.defchararray.isdigit(pos[:, 0])
    pos = pos[i_ok, :]
    data = data[i_ok, :]
    pos = pos.astype('int')

    # sort everything nicely
    sidx = sp.lexsort((pos[:, 2], pos[:, 1], pos[:, 0]))
    pos = pos[sidx, :]
    data = data[sidx, :]

    # find genes with overlapping annotations
    my_overlap_genes = []
    for i in xrange(pos.shape[0]):
        my_pos = pos[i, :]

        # same chr
        i_chr = my_pos[0] == pos[:, 0]

        # end is in something else
        i_lb_end = my_pos[2] >= pos[:, 1]
        i_ub_end = my_pos[2] <= pos[:, 2]

        # st is in something else
        i_lb_st = my_pos[1] <= pos[:, 2]
        i_ub_st = my_pos[1] >= pos[:, 1]

        # on both ends the only entry that overlaps to i is i itself --> continue
        if (sp.sum(i_chr & i_lb_end & i_ub_end) == 1) and (sp.sum(i_chr & i_lb_st & i_ub_st) == 1):
            continue

        # extract IDs of overlapping genes
        overlap_genes_st = data[i_chr & i_ub_st & i_lb_st, 0]
        overlap_genes_end = data[i_chr & i_ub_end & i_lb_end, 0]

        overlap_genes_st = sp.array([x.split('|')[0] for x in overlap_genes_st])
        overlap_genes_end = sp.array([x.split('|')[0] for x in overlap_genes_end])

        # this should actually never happen ...
        if (sp.unique(overlap_genes_st).shape[0] == 1) and (sp.unique(overlap_genes_end).shape[0] == 1):
            continue
        if sp.unique(overlap_genes_st).shape[0] > 1:
            my_overlap_genes.extend(overlap_genes_st.tolist())
        if sp.unique(overlap_genes_end).shape[0] > 1:
            my_overlap_genes.extend(overlap_genes_end.tolist())
    return sp.unique(my_overlap_genes)


def reading_anno(fn_anno, overlap_genes, protein_coding_filter):
    """
    Reads in all transcript annotations,
    removes overlapping genes and
    removes non-chr-contigs and
    eventually filters for protein-coding genes on the fly
    """
    # for removing non-chr-contigs later
    chr_whitelist = [str(x) for x in range(NMB_CHR)]
    chr_whitelist.extend(['chr%i' % i for i in range(NMB_CHR)])
    chr_whitelist.extend(['chrx', 'chry', 'chrm', 'x', 'y', 'm', 'mt'])

    data = dict()
    # collect transcript information
    transcripts = dict()
    for l in open(fn_anno, 'r'):
        if l[SEQ_NAME] == '#':
            continue

        l_spl = l.strip('\n').split('\t')

        if l_spl[FEATURE].lower() != 'exon':
            continue
        tags = dict()
        for t in l_spl[ATTRIBUTE].strip(';').split(';'):
            tt = t.strip(' ').split(' ')
            tags[tt[0]] = tt[1].strip('"')
        try:
            transcripts[tags['transcript_id']].append('-'.join([l_spl[START], l_spl[END]]))
        except KeyError:
            transcripts[tags['transcript_id']] = ['-'.join([l_spl[START], l_spl[END]])]

    # read transcript annotation
    for l in open(fn_anno, 'r'):
        if l[SEQ_NAME] == '#':
            continue
        l_spl = l.strip('\n').split('\t')

        if l_spl[FEATURE].lower() != 'transcript':
            continue
        tags = dict()
        for t in l_spl[ATTRIBUTE].strip(';').split(';'):
            tt = t.strip(' ').split(' ')
            tags[tt[0]] = tt[1].strip('"')

        key = tags['gene_id']
        gene_type = tags['gene_type']

        # remove overlapping genes
        if key in overlap_genes:
            continue
        # remove non-chr-contigs
        if l_spl[SEQ_NAME].lower() in chr_whitelist:
            continue
        # filter for protein-coding genes only (if filter turned on)
        if protein_coding_filter and (gene_type != "protein_coding"):
            continue

        value = '%s:%s:%s' % (l_spl[SEQ_NAME], ','.join(transcripts[tags['transcript_id']]), l_spl[STRAND])

        try:
            data[key].append(value)
        except KeyError:
            data[key] = [value]

    return data


def get_transcript_length(rec):
    """
        Returns transcript length defined as sum of length of exons
    """
    ex_pieces = sp.array(rec.split(':')[1].split(','))
    lgt = 0

    for i, x in enumerate(ex_pieces[0:]):
        start, end = x.split('-')
        lgt += int(end) - int(start) + 1
    assert lgt != 0, "Outch, there are transcripts with no length"
    return lgt


def process_single_transcript_genes(tcrpt):
    assert len(tcrpt) == 1, "Too many transcripts to process"

    # checking that we have at least two exons
    tcrpt = tcrpt[0]
    if tcrpt.find(',') == -1:
        return None

    # reformat to somewhat convenient reading
    # format # ID : exon1_positions,exon2_positions,...,exonN_positions : STRAND

    first_ex = tcrpt.split(':')[0] + ':' + tcrpt.split(':')[1].split(',')[0] + ':' + tcrpt.split(':')[2]
    last_ex = tcrpt.split(':')[0] + ':' + tcrpt.split(':')[1].split(',')[-1] + ':' + tcrpt.split(':')[2]
    return [first_ex, last_ex, tcrpt.split(':')[0], tcrpt.split(':')[2], get_transcript_length(tcrpt)]


def process_multi_transcript_genes(tcrpt):
    # We only use transcript isoforms that have at least two exons
    if sp.sum(np.core.defchararray.find(tcrpt, ',') != -1) != len(tcrpt):
        # TODO: so we not just exclude those transcripts but ignore genes where one or more transcripts have less than 2 exons?!
        return None

    # make matrix of transcript struct and length
    my_exons = [x.split(':')[1].split(',') for x in tcrpt]
    # unravel exons into one list of exons
    my_exons = sp.array([reduce(lambda x, y: x + y, my_exons)]).ravel()
    my_exons_int = sp.array([x.split('-') for x in my_exons]).astype('int')

    # sort this
    sidx_int = sp.lexsort((my_exons_int[:, 1], my_exons_int[:, 0]))
    my_exons = my_exons[sidx_int]
    my_exons_int = my_exons_int[sidx_int, :]

    # see how often we got each item
    dummy, u_idx, dists = unique_rows(my_exons_int, index=True, counts=True)
    n_match = sp.sum(dists == len(tcrpt))

    # make sure we have at least 3 constitutive exons
    if n_match < 3:
        return None

    # get constitutive exons
    i_const = dists == len(tcrpt)
    uq_const_ex = my_exons[u_idx][i_const]

    first_ex = uq_const_ex[0]
    last_ex = uq_const_ex[-1]

    # get length of all transcripts
    my_ex_struct_l = []
    for i, rec in enumerate(tcrpt):
        my_ex_struct_l.append(get_transcript_length(rec))

    first_ex = tcrpt[0].split(':')[0] + ':' + first_ex + ':' + tcrpt[0].split(':')[2]
    last_ex = tcrpt[0].split(':')[0] + ':' + last_ex + ':' + tcrpt[0].split(':')[2]
    return [first_ex, last_ex, tcrpt[0].split(':')[0], tcrpt[0].split(':')[2], str(sp.median(my_ex_struct_l))]


def read_annotation_file(fn_anno, protein_coding_filter):
    # get list of overlapping genes
    overlap_genes = get_overlap_genes(fn_anno)

    # reading annotation file in
    data = reading_anno(fn_anno, overlap_genes, protein_coding_filter)

    uq_g_id = data.keys()  # unique gene ids
    new_data = []
    for gid in uq_g_id:
        # process transcripts
        if len(data[gid]) == 1:
            temp = process_single_transcript_genes(data[gid])
        else:
            temp = process_multi_transcript_genes(data[gid])

        # make sure it has been processed correctly
        if temp is None:
            continue
        else:
            temp.extend([gid])
            new_data.append(temp)
    new_data = sp.array(new_data)
    print new_data.shape
    print new_data[0:5,:]
    s_idx = sp.argsort(new_data[:, 5])
    new_data = new_data[s_idx, :]
    # filter gene with no name
    return sp.array(new_data)


def get_annotation_table(fn_anno, protein_coding_filter):

    if fn_anno.lower().endswith('gtf'):
        exon_t_gene = read_annotation_file(fn_anno, protein_coding_filter)
    else:
        raise Exception(
            "Only annotation files in format gtf are supported. File name must end accordingly")

    return exon_t_gene


def get_counts_from_single_bam(fn_bam, regions):
    """This function extracts read counts from a given bam file spanning
       a set of given intervals."""

    if not os.path.exists(fn_bam + '.bai'):
        warnings.warn('WARNING: alignment file %s seems not to be indexed and will be skipped! \n' % fn_bam)
        dummy = sp.zeros(regions.shape[0] * 2)
        dummy[:] = sp.nan
        return dummy
    if not os.stat(fn_bam).st_size > 0:
        warnings.warn('WARNING: alignment file %s seems to be empty and will be skipped! \n' % fn_bam)
        dummy = sp.zeros(regions.shape[0] * 2)
        dummy[:] = sp.nan
        return dummy

    bam_file = pysam.Samfile(fn_bam, 'rb')
    ref_seqs = bam_file.references
    cnts = sp.zeros((regions.shape[0], 2), dtype='float')
    t0 = time.time()

    if len(regions.shape) > 1:
        sidx = sp.argsort(regions[:, 0])
    else:
        sidx = sp.argsort(regions)

    for i, ii in enumerate(sidx):
        rec = regions[ii]
        if i > 0 and i % 100 == 0:
            print '%i rounds to go. ETA %.0f seconds' % (regions.shape[0] - i, (time.time() - t0) / i * (regions.shape[0] - i))
        if len(regions.shape) == 1:
            chrm = rec.split(':')[0]
            if chrm not in ref_seqs:
                chrm = chrm.strip('chr')
            start1 = int(rec.split(':')[1].split('-')[0])
            end1 = int(rec.split(':')[1].split('-')[1])
            start2 = None
            end2 = None
        else:
            chrm = rec[0].split(':')[0]
            if chrm not in ref_seqs:
                chrm = chrm.strip('chr')
            start1 = int(rec[0].split(':')[1].split('-')[0])
            end1 = int(rec[0].split(':')[1].split('-')[1])
            start2 = int(rec[1].split(':')[1].split('-')[0])
            end2 = int(rec[1].split(':')[1].split('-')[1])
        try:
            cnt1 = int(sp.ceil(sp.sum(
                [sp.sum((sp.array(read.positions) >= start1) & (sp.array(read.positions) < end1)) for read in
                 bam_file.fetch(chrm, start1, end1) if not read.is_secondary]) / 50.0))
            if start2 is None:
                cnt2 = cnt1
            else:
                cnt2 = int(sp.ceil(sp.sum(
                    [sp.sum((sp.array(read.positions) >= start2) & (sp.array(read.positions) < end2)) for read in
                     bam_file.fetch(chrm, start2, end2) if not read.is_secondary]) / 50.0))
        except ValueError:
            print >> sys.stderr, 'Ignored %s' % chrm
            cnt1 = 1
            cnt2 = 1
        finally:
            cnts[ii, :] = [cnt1, cnt2]
    bam_file.close()

    return cnts.ravel('C')


def get_counts_from_multiple_bam(fn_bams, regions):
    """ This is a wrapper to concatenate counts for a given list of bam
        files"""

    if len(fn_bams) == 1:
        return get_counts_from_single_bam(fn_bams[0], regions)[:, sp.newaxis]
    else:
        return sp.hstack([get_counts_from_single_bam(fn_bams[i], regions)[:, sp.newaxis] for i in range(len(fn_bams))])


def parse_options(argv):
    parser = OptionParser()

    sample_input = OptionGroup(parser, 'Input')
    sample_input.add_option('', '--bam_dir', dest='dir_bam', metavar='FILE', help='Directory of bam files', default='-')
    sample_input.add_option('', '--bam_fn', dest='fn_bam', metavar='FIlE', help='Specifies single bam file', default='-')
    sample_input.add_option('', '--anno_fn', dest='fn_anno', metavar='FILE', help='Annotation', default='-')

    opt_gen = OptionGroup(parser, 'General Options')
    opt_gen.add_option('', '--protein_coding_filter_OFF', dest="proteinCodingFilter", action="store_false",
                       help="Consider only genes that are protein-coding", default=True)

    parser.add_option_group(sample_input)
    parser.add_option_group(opt_gen)
    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)
    if sp.sum(int(options.dir_bam != '-') + int(options.fn_bam != '-')) != 1:
        print "Please specify either bam file or directory"
        parser.print_help()
        sys.exit(2)
    return options


def main():
    options = parse_options(sys.argv)
    exon_t_gene = get_annotation_table(options.fn_anno, options.proteinCodingFilter)

    if options.dir_bam != '-':
        file_names = glob.glob(os.path.join(options.dir_bam, '*.bam'))
        data = get_counts_from_multiple_bam(file_names, exon_t_gene)
    elif options.fn_bam != '-':
        file_names = [options.fn_bam]
        data = get_counts_from_multiple_bam(file_names, exon_t_gene)

    # Normalize counts by exon length
    exon_l = sp.array([int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) + 1 for x in exon_t_gene[:, :2].ravel('C')],
                      dtype='float') / 1000.
    data /= sp.tile(exon_l[:, sp.newaxis], data.shape[1])

    # Get counts from first and last exon
    my_counts = get_counts_from_marginal_exons(exon_t_gene, data)

    avg_count_per_exon()


if __name__ == "__main__":
    main()
