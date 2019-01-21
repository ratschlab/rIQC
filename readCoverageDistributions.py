import os
import sys
import glob
import pickle
from optparse import OptionParser, OptionGroup
import pysam
import time
import numpy as np
from libs.usefulTools import *  # brings packages: warnings, scipy as sp
import matplotlib.pyplot as plt

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
def avg_count_per_exon(counts, regions, file_name):
    for gene in regions:
        if not gene[0] in counts:
            continue
        val = counts[gene[0]]
        if not np.sum(val[:, 3] > 0.000):
            continue
        if gene[2] == "+":
            start = val[0][0]  # first position
            end = val[-1][1]  # last position
        else:
            end = val[0][1]  # first position
            start = val[-1][0]  # last position
        interval = (end - start) / 100
        for ex in val:
            rel_pos = (ex[0]-start)/interval + ((ex[1]-start)/interval - (ex[0]-start)/interval) / 2  # middle position of exon normalized to [0,100]
            if rel_pos < 0 or rel_pos > 100:
                print "Oh"
            if start > end:
                print "Oha"
            plt.plot(rel_pos, ex[3], ".", color="firebrick")

    plt.title("Constitutive Exons - all lengths (%s)" % file_name)
    plt.ylabel("Normalized Count")
    plt.xlabel("Relative Location")

    plt.savefig("../2018_degradationPaper/Coverage/constitutive_all_%s_filtered" %file_name)
    plt.show()


# Expression distribution over normalized gene length (in different length-bins-> same as already used; constitutive and not)
# # Generate average coverage distribution across all genes within that bin after projecting all genes in bin onto one length
def distr_over_gene_lengths(counts, regions, file_name):
    gene_lengths = regions[:, 3]
    nmb_genes = gene_lengths.shape[0]

    nmb_bins = 10
    for b in range(nmb_bins):
        idx_s = np.argsort(gene_lengths)
        low_b = nmb_genes / nmb_bins * b
        up_b = nmb_genes / nmb_bins * (b + 1)

        low_l = gene_lengths[idx_s[low_b]]
        up_l = gene_lengths[idx_s[up_b]]

        for gene in regions:
            if not gene[0] in counts:
                continue
            if not gene[3] >= low_l and gene[3] < up_l:
                continue
            val = counts[gene[0]]
            if gene[2] == "+":
                start = val[0][0]  # first position
                end = val[-1][1]  # last position
            else:
                end = val[0][1]  # first position
                start = val[-1][0]  # last position
            interval = (end - start) / 100
            for ex in val:
                rel_pos = (ex[0]-start)/interval + ((ex[1]-start)/interval - (ex[0]-start)/interval) / 2  # middle position of exon normalized to [0,100]
                if rel_pos < 0 or rel_pos > 100:
                    print "Oh"
                if start > end:
                    print "Oha"
                plt.plot(rel_pos, ex[3], ".", color="blue")

        plt.title("Constitutive Exons - bin %i (%s)" % (b+1, file_name))
        plt.ylabel("Normalized Count")
        plt.xlabel("Relative Location")

        plt.savefig("../2018_degradationPaper/Coverage/constitutive_bin%i_%s" % (b+1, file_name))
        plt.show()


# Distribution of length of last-constitutive exons (joint and in different length bins)
def last_const_exon_length():
    print "Hi"


# Position of last constitutive exon in the list of all exons (based on normalized gene length)
def last_const_exon_pos():
    print "Hi"

###############################################################################################


def get_overlap_genes(fn_anno):
    """
    Returns a list of gene names which are overlapping
    """

    data = []
    # For all gene entries
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

    data = dict()  # a dictionary with gene IDs as keys and a list of transcripts with format CHR:listOfExons(e1,e2,...):STRAND as values
    # collect transcript information
    transcripts = dict()  # a dictionary with transcript IDs as keys and a list of exons (start-end) as value
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
        if l_spl[SEQ_NAME].lower() not in chr_whitelist:
            continue
        # filter for protein-coding genes only (if filter turned on)
        if protein_coding_filter and (gene_type != "protein_coding"):
            continue

        value = '%s:%s:%s' % (l_spl[SEQ_NAME], ','.join(transcripts[tags['transcript_id']]), l_spl[STRAND])

        try:
            data[key].append(value)
        except KeyError:
            data[key] = [value]

    return data  # a dictionary with gene-ID as key and a list of transcripts, each in the form
                 # "chr_name:exon1_start-exon1_end,...,exonN_start-exonN_end:strand" as value


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
        return None, None

    # format is ID:exon1_positions,exon2_positions,...,exonN_positions:STRAND

    exons = tcrpt.split(':')[1].split(',')  # as list

    return [tcrpt.split(':')[0], tcrpt.split(':')[2], get_transcript_length(tcrpt)], exons


def process_multi_transcript_genes(tcrpt):
    # We only use transcript isoforms that have at least two exons
    if sp.sum(np.core.defchararray.find(tcrpt, ',') != -1) != len(tcrpt):
        return None, None, None

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

    # make sure we have at least 2 constitutive exons
    if n_match < 2:
        return None, None, None

    # get constitutive exons
    i_const = dists == len(tcrpt)
    uq_const_ex = my_exons[u_idx][i_const]

    # get length of all transcripts
    my_ex_struct_l = []
    for i, rec in enumerate(tcrpt):
        my_ex_struct_l.append(get_transcript_length(rec))

    all_exons = []  # TODO

    return [tcrpt[0].split(':')[0], tcrpt[0].split(':')[2], str(sp.median(my_ex_struct_l))], uq_const_ex, all_exons


def read_annotation_file(fn_anno, protein_coding_filter):
    # get list of overlapping genes
    overlap_genes = get_overlap_genes(fn_anno)

    # reading in annotation file
    data = reading_anno(fn_anno, overlap_genes, protein_coding_filter)
    # data is a dictionary with gene IDs as keys and a list of transcripts with format CHR:listOfExons(e1,e2,...):STRAND as values

    uq_g_id = data.keys()  # unique gene ids
    new_data = []

    const_exons = dict()
    # all_exons = dict()  # TODO
    for gid in uq_g_id:
        # process transcripts
        if len(data[gid]) == 1:
            temp, c_e = process_single_transcript_genes(data[gid])
            # a_e = c_e  # TODO
        else:
            temp, c_e, a_e = process_multi_transcript_genes(data[gid])

        # make sure it has been processed correctly
        if temp is None:
            continue
        else:
            const_exons[gid] = c_e
            # all_exons[gid] = a_e  # TODO
            new_data.append([gid] + temp)
    new_data = sp.array(new_data)
    s_idx = sp.argsort(new_data[:, -1])
    new_data = new_data[s_idx, :]
    # filter gene with no name

    # new_data is array with entries: gene_ID, seq_name, strand, (median) transcript-length
    return sp.array(new_data), const_exons  # , all_exons TODO


def get_annotation_table(fn_anno, protein_coding_filter):

    if fn_anno.lower().endswith('gtf'):
        exon_t_gene, const_exons = read_annotation_file(fn_anno, protein_coding_filter)
    else:
        raise Exception(
            "Only annotation files in format gtf are supported. File name must end accordingly")

    return exon_t_gene, const_exons


def get_counts_from_single_bam(fn_bam, regions, exons):
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
    cnts = dict()
    t0 = time.time()

    # Sort regions by chr
    if len(regions.shape) > 1:
        sidx = sp.argsort(regions[:, 1])  # seq_name
    else:
        exit("regions.shape was not > 1")
        # sidx = sp.argsort(np.vstack((regions[1], regions[0])))

    for i, ii in enumerate(sidx):
        if i > 0 and i % 100 == 0:
            print '%i rounds to go. ETA %.0f seconds' % (regions.shape[0] - i, (time.time() - t0) / i * (regions.shape[0] - i))

        rec = regions[ii]  # e.g. ENSG00000233493.3_2	chr19	-	1000
        rec_exons = exons[rec[0]]  # rec[0] is unique gene ID

        if rec[2] == "-" and int(rec_exons[0].split("-")[0]) < int(rec_exons[-1].split('-')[0]):
            rec_exons = np.flipud(rec_exons)

        exon_counts = np.zeros((len(rec_exons), 4), dtype=float)  # store start, end, length, count

        if len(regions.shape) == 1:
            exit("regions.shape == 1")
        # else
        chrm = rec[1]
        if chrm not in ref_seqs:
            chrm = chrm.strip('chr')
            if chrm not in ref_seqs:
                exit("%s is not in bam-references" % chrm)

        for e in range(len(rec_exons)):
                start = int(rec_exons[e].split("-")[0])
                end = int(rec_exons[e].split("-")[1])
                exon_counts[e, 0] = start
                exon_counts[e, 1] = end
                exon_counts[e, 2] = end - start
                cnt = 1
                try:
                    cnt = int(sp.ceil(sp.sum(
                        [sp.sum((sp.array(read.positions) >= start) & (sp.array(read.positions) < end)) for read in
                         bam_file.fetch(str(chrm), start, end) if not read.is_secondary]) / 50.0))
                except ValueError:
                    print >> sys.stderr, 'Ignored %s' % chrm
                    cnt = 1
                finally:
                    exon_counts[e, 3] = cnt / exon_counts[e, 2]
        cnts[rec[0]] = exon_counts
    bam_file.close()

    return cnts


def get_counts_from_multiple_bam(fn_bams, regions, exons):
    """ This is a wrapper to concatenate counts for a given list of bam
        files"""

    if len(fn_bams) == 1:
        return [get_counts_from_single_bam(fn_bams[0], regions, exons)]
    else:
        li = []
        for i in range(len(fn_bams)):
            li.append(get_counts_from_single_bam(fn_bams[i], regions, exons))
        return li


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
    if os.path.exists("./anno.tmp") and os.path.exists("./const_ex.pkl"):
        exon_t_gene = sp.loadtxt("./anno.tmp", delimiter='\t', dtype='string')
        const_exons = pickle.load(open("./const_ex.pkl", "rb"))
    else:
        exon_t_gene, const_exons = get_annotation_table(options.fn_anno, options.proteinCodingFilter)
        sp.savetxt("./anno.tmp", exon_t_gene, delimiter='\t', fmt='%s')
        f = open("./const_ex.pkl", "wb")
        pickle.dump(const_exons, f)
        f.close()

    if not os.path.exists("./count_data.pkl"):
        if options.dir_bam != '-':
            file_names = glob.glob(os.path.join(options.dir_bam, '*.bam'))
            data = get_counts_from_multiple_bam(file_names, exon_t_gene, const_exons)
        else:
            assert options.fn_bam != '-'
            file_names = [options.fn_bam]
            data = get_counts_from_multiple_bam(file_names, exon_t_gene, const_exons)
        f = open("./count_data.pkl", "wb")
        pickle.dump(data, f)
        f.close()
    else:
        file_names = ['FFPE_1', 'FFPE_2', 'FFPE_3', 'FFPE_4', 'FF_1', 'FF_2', 'FF_3', 'FF_4']
        data = pickle.load(open("./count_data.pkl", "rb"))
        # data is a dictionary with unique gene_IDs as keys and
        # a list of (constitutive) exons (start, end, (normalized) count) as value

    for i in range(len(file_names)):
        # distr_over_gene_lengths(data[i], exon_t_gene, file_names[i])
        avg_count_per_exon(data[i], exon_t_gene, file_names[i])


if __name__ == "__main__":
    main()
