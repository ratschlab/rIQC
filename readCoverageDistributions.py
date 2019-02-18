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
ANNO_SEQ_NAME = 0    # name of chromosome or scaffold
ANNO_SOURCE = 1      # name of program that generated this feature
ANNO_FEATURE = 2     # feature type name (e.g. "gene", "transcript", "exon") / type of feature (term or accession from SOFA sequence ontology)
ANNO_START = 3       # start position of feature (seq numbering starting at 1)
ANNO_END = 4         # end position of feature (seq numbering starting at 1)
ANNO_SCORE = 5       # a floating point value
ANNO_STRAND = 6      # + (forward) or - (reverse)
ANNO_FRAME = 7       # 0/1/2 : position in codon
ANNO_ATTRIBUTE = 8   # semicolon-separated list of tag-value pairs

#  COUNTS
CNT_ST = 0
CNT_END = 1
CNT_LEN = 2
CNT_CNT = 3

#  REGIONS / EXON_T_GENE
REG_ID = 0
REG_CHR = 1
REG_STR = 2
REG_LEN = 3

DEG = ['TCGA-G9-6498-01A-12R-A311-07', 'TCGA-EJ-A46F-01A-31R-A250-07', 'TCGA-KK-A59V-01A-11R-A29R-07', 'TCGA-G9-6362-01A-11R-1789-07', 'TCGA-J4-A67O-01A-11R-A30B-07']
NON_DEG = ['TCGA-EJ-7789-01A-11R-2118-07', 'TCGA-J4-AATZ-01A-11R-A41O-07', 'TCGA-EJ-7781-01A-11R-2118-07', 'TCGA-EJ-5511-01A-01R-1580-07', 'TCGA-EJ-7321-11A-01R-2263-07']


def prepare_outlier_filter(counts, regions):
    all_counts = []
    for gene in regions:
        if not gene[REG_ID] in counts:
            continue  # we have no counts for that gene
        val = counts[gene[REG_ID]]
        for ex in val:
            all_counts.append(ex[CNT_CNT])
    return np.percentile(all_counts, q=99.5)


def avg_count_per_exon(counts, regions, file_name, cut_off, histo, gn_version):

    coords = dict()

    for gene in regions:
        if not gene[REG_ID] in counts:
            continue  # we have no counts for that gene
        val = counts[gene[REG_ID]]
        if not np.sum(val[:, CNT_CNT] > 0.000):
            continue  # none of the exons is longer than 0

        last_end = val[0][CNT_ST] - 1  # need this for "cutting together" exons

        take_gene = True
        for ex in val:
            if ex[CNT_CNT] > cut_off:
                take_gene = False
                break
            if ex[CNT_ST] > last_end:
                ex[CNT_END] = ex[CNT_END] - (ex[CNT_ST] - (last_end + 1))
                ex[CNT_ST] = last_end + 1
                last_end = ex[CNT_END]
            # if exons overlap:
            else:
                print "CAUTION: Exons overlap"
                last_end = ex[CNT_END]

        if not take_gene:
            continue

        start = val[0][CNT_ST]
        end = val[-1][CNT_END]
        interval = (end - start) / 1000.0

        xs = []
        ys = []
        bases = np.zeros(1001)
        for ex in val:
            if ex[CNT_ST] < start:
                print "Error: start was not smallest index"
            if ex[CNT_END] > end:
                print "Error: end was not biggest index"

            rel_pos_start = int(np.ceil((ex[CNT_ST] - start) / interval))
            rel_pos_end = int(np.floor((ex[CNT_END] - start) / interval))
            if gene[REG_STR] == "-":
                rel_pos_start, rel_pos_end = 1000 - rel_pos_end, 1000 - rel_pos_start
            if rel_pos_start < 0 or rel_pos_start > 1000 or rel_pos_end < 0 or rel_pos_end > 1000:
                print "Start-Position in unexpected range " + str(rel_pos_start) + " or End-Position " + str(rel_pos_end)
            if start > end or rel_pos_start > rel_pos_end:
                print "Made a mistake at determining start and end"

            if histo:
                for i in range(rel_pos_start, rel_pos_end+1, 1):
                    bases[i] = ex[CNT_CNT]
            else:
                # we need this distinction for the graph - otherwise we have "zick-zack" pattern (--> we need to go from lower to higher positions)
                if gene[REG_STR] == "-":
                    xs.insert(0, rel_pos_end)
                    ys.insert(0, ex[CNT_CNT])
                    xs.insert(0, rel_pos_start)
                    ys.insert(0, ex[CNT_CNT])
                else:
                    xs.append(rel_pos_start)
                    ys.append(ex[CNT_CNT])
                    xs.append(rel_pos_end)
                    ys.append(ex[CNT_CNT])
        if histo:
            for i in range(0, 10, 1):
                try:
                    coords[i].append(np.median(bases[100*i:100*(i+1)]))
                except KeyError:
                    coords[i] = [np.median(bases[100*i:100*(i+1)])]

        else:
            plt.plot(xs, ys, "--", color="firebrick")

    if histo:

        for key in coords:
            p1 = plt.bar((100 * key), len(filter(lambda ex_l: (ex_l > 0.5), coords[key])), width=25, color='royalblue', align='edge')
            p2 = plt.bar((100 * key) + 25, len(filter(lambda ex_l: (ex_l > 1), coords[key])), width=25, color='firebrick', align='edge')
            p3 = plt.bar((100 * key) + 50, len(filter(lambda ex_l: (ex_l > 2), coords[key])), width=25, color='orange', align='edge')
            p4 = plt.bar((100 * key) + 75, len(filter(lambda ex_l: (ex_l > 4), coords[key])), width=25, color='olivedrab', align='edge')

        plt.xlim(0, 1000)
        plt.ylim(0, 380)  # FFPE_VS_FF
        # plt.ylim(0, 550)
        plt.legend([p1, p2, p3, p4], ["> 0.5", "> 1", "> 2", "> 4"])
        plt.title("Constitutive Exons - all lengths (%s)" % file_name)
        plt.ylabel("Number of exons longer than threshold")
        plt.xlabel("Relative Location (10 areas)")

        plt.savefig("../2018_degradationPaper/ReadCoverage/" + gn_version + "/histo_median_constExons_allLengths_%s" % file_name)
        plt.show()

    else:
        plt.xlim(0, 1000)
        plt.ylim(0, max(cut_off, 40))
        plt.title("Constitutive Exons - all lengths (%s)" % file_name)
        plt.ylabel("Normalized Count")
        plt.xlabel("Relative Location")

        plt.savefig("../2018_degradationPaper/ReadCoverage/" + gn_version + "/constExons_allLengths_%s" % file_name)
        plt.show()


def distr_over_gene_lengths(counts, regions, file_name, cut_off, histo, gn_version):
    gene_lengths = regions[:, 3]
    nmb_genes = gene_lengths.shape[0]
    idx_s = np.argsort(gene_lengths)

    nmb_bins = 4
    for b in range(nmb_bins):
        low_b = nmb_genes / nmb_bins * b
        up_b = nmb_genes / nmb_bins * (b + 1) - 1

        low_l = gene_lengths[idx_s[low_b]]
        up_l = gene_lengths[idx_s[up_b]]

        coords = dict()

        for gene in regions:
            if not gene[REG_ID] in counts:
                continue  # we have no counts for that gene
            if not (low_l <= gene[REG_LEN] < up_l):
                continue
            val = counts[gene[REG_ID]]
            if not np.sum(val[:, CNT_CNT] > 0.000):
                continue  # none of the exons is longer than 0

            last_end = val[0][CNT_ST] - 1  # need this for "cutting together" exons

            take_gene = True
            for ex in val:
                if ex[CNT_CNT] > cut_off:
                    take_gene = False
                    break
                if ex[CNT_ST] > last_end:
                    ex[CNT_END] = ex[CNT_END] - (ex[CNT_ST] - (last_end + 1))
                    ex[CNT_ST] = last_end + 1
                    last_end = ex[CNT_END]
                # if exons overlap:
                else:
                    print "CAUTION: Exons overlap"
                    last_end = ex[CNT_END]

            if not take_gene:
                continue

            start = val[0][CNT_ST]
            end = val[-1][CNT_END]
            interval = (end - start) / 1000.0

            xs = []
            ys = []
            bases = np.zeros(1001)
            for ex in val:
                if ex[CNT_ST] < start:
                    print "Error: start was not smallest index"
                if ex[CNT_END] > end:
                    print "Error: end was not biggest index"

                rel_pos_start = int(np.ceil((ex[CNT_ST] - start) / interval))
                rel_pos_end = int(np.floor((ex[CNT_END] - start) / interval))
                if gene[REG_STR] == "-":
                    rel_pos_start, rel_pos_end = 1000 - rel_pos_end, 1000 - rel_pos_start
                if rel_pos_start < 0 or rel_pos_start > 1000 or rel_pos_end < 0 or rel_pos_end > 1000:
                    print "Start-Position in unexpected range " + str(rel_pos_start) + " or End-Position " + str(rel_pos_end)
                if start > end or rel_pos_start > rel_pos_end:
                    print "Made a mistake at determining start and end"

                if histo:
                    for i in range(rel_pos_start, rel_pos_end+1, 1):
                        bases[i] = ex[CNT_CNT]
                else:
                    if gene[REG_STR] == "-":
                        xs.insert(0, rel_pos_end)
                        ys.insert(0, ex[CNT_CNT])
                        xs.insert(0, rel_pos_start)
                        ys.insert(0, ex[CNT_CNT])
                    else:
                        xs.append(rel_pos_start)
                        ys.append(ex[CNT_CNT])
                        xs.append(rel_pos_end)
                        ys.append(ex[CNT_CNT])
            if histo:
                for i in range(0, 10, 1):
                    try:
                        coords[i].append(np.median(bases[100 * i:100 * (i + 1)]))
                    except KeyError:
                        coords[i] = [np.median(bases[100 * i:100 * (i + 1)])]

            else:
                plt.plot(xs, ys, "--", color="firebrick")

        if histo:

            for key in coords:
                p1 = plt.bar((100 * key), len(filter(lambda ex_l: (ex_l > 0.5), coords[key])), width=25, color='royalblue', align='edge')
                p2 = plt.bar((100 * key) + 25, len(filter(lambda ex_l: (ex_l > 1), coords[key])), width=25, color='firebrick', align='edge')
                p3 = plt.bar((100 * key) + 50, len(filter(lambda ex_l: (ex_l > 2), coords[key])), width=25, color='orange', align='edge')
                p4 = plt.bar((100 * key) + 75, len(filter(lambda ex_l: (ex_l > 4), coords[key])), width=25, color='olivedrab', align='edge')

            plt.xlim(0, 1000)
            plt.ylim(0, 130)
            # plt.ylim(0, 170)
            plt.legend([p1, p2, p3, p4], ["> 0.5", "> 1", "> 2", "> 4"])
            plt.title("Constitutive Exons - bin %i (%s)" % (b + 1, file_name))
            plt.ylabel("Number of exons longer than threshold")
            plt.xlabel("Relative Location (10 areas)")

            plt.savefig("../2018_degradationPaper/ReadCoverage/" + gn_version + "/histo_median_constitutive_%s_bin%i" % (file_name, b + 1))
            plt.show()

        else:
            plt.xlim(0, 1000)
            plt.ylim(0, max(cut_off, 40))
            plt.title("Constitutive Exons - bin %i (%s)" % (b + 1, file_name))
            plt.ylabel("Normalized Count")
            plt.xlabel("Relative Location")

            plt.savefig("../2018_degradationPaper/ReadCoverage/" + gn_version + "/constitutive_%s_bin%i" % (file_name, b + 1))
            plt.show()


# Distribution of length of last-constitutive exons (joint and in different length bins)
def last_const_exon_length(counts, regions, file_name, gn_version):
    end_counts = []
    for gene in regions:
        if not gene[REG_ID] in counts:
            continue  # we have no counts for that gene
        val = counts[gene[REG_ID]]
        if not np.sum(val[:, CNT_CNT] > 0.000):
            continue  # none of the exons is longer than 0

        if gene[REG_STR] == "-":
            end_counts.append(min(val[0][CNT_LEN], 2000))
        else:
            end_counts.append(min(val[-1][CNT_LEN], 2000))
    n, bins, patches = plt.hist(end_counts, bins=100, color="slateblue")
    plt.title("Length of last exon (%s)" % file_name)
    plt.savefig("../2018_degradationPaper/ReadCoverage/" + gn_version + "/lastExonLength_%s" % file_name)

    plt.show()


# Position of last constitutive exon in the list of all exons (based on normalized gene length)
def last_const_exon_pos(counts, regions, file_name, exon_lookup, gn_version):
    # exon_lookup is a dictionary with gene-IDs as key and a list as value that contains for each transcript a string of the form:
    # EX1_START-EX1_END,EX2_START-EX2_END,.....
    end_pos = []
    excl_count = 0
    for gene in regions:
        if not gene[REG_ID] in counts:
            continue  # we have no counts for that gene
        val = counts[gene[REG_ID]]
        if not np.sum(val[:, CNT_CNT] > 0.000):
            continue  # none of the exons is longer than 0

        # if we only have 1 transcript
        if len(exon_lookup[gene[REG_ID]]) == 1:
            trcpt = exon_lookup[gene[REG_ID]][0].split(",")
            trcpt = np.asarray([x.split("-") for x in trcpt]).astype(int)
            # at this point, "trcpt" is an array of arrays of start- and end-position for each exon in the transcript
            last_ex_ind = 0
            last_ex = np.zeros(2)
            for i in range(len(trcpt)):
                if (gene[REG_STR] == "-" and val[0][CNT_ST] == trcpt[i, 0] and val[0][CNT_END] == trcpt[i, 1]) or \
                        (gene[REG_STR] == "+" and val[-1][CNT_ST] == trcpt[i, 0] and val[-1][CNT_END] == trcpt[i, 1]):
                    last_ex_ind = i
                    break

            if gene[REG_STR] == '-':
                loop_range = range(len(trcpt)-1, -1, -1)
                last_end = trcpt[-1][0] - 1
            else:
                loop_range = range(0, len(trcpt), 1)
                last_end = trcpt[0][0] - 1

            take_gene = True
            for i in loop_range:
                if trcpt[i, 0] > last_end:
                    trcpt[i, 1] = trcpt[i, 1] - (trcpt[i, 0] - (last_end + 1))
                    trcpt[i, 0] = last_end + 1
                    last_end = trcpt[i, 1]
                # if exons overlap:
                else:
                    print "CAUTION: Exons overlap"
                    last_end = trcpt[i, 1]
                if i == last_ex_ind:
                    last_ex = trcpt[i]

            if not take_gene:
                continue

            if gene[REG_STR] == '-':
                start = trcpt[-1, 0]
                end = trcpt[0, 1]
                interval = (end - start) / 1000.0
                rel_pos_start = 1000 - int(np.floor((last_ex[1] - start) / interval))
                rel_pos_end = 1000 - int(np.ceil((last_ex[0] - start) / interval))
            else:
                start = trcpt[0, 0]
                end = trcpt[-1, 1]
                interval = (end - start) / 1000.0
                rel_pos_start = int(np.ceil((last_ex[0] - start) / interval))
                rel_pos_end = int(np.floor((last_ex[1] - start) / interval))

            end_pos.append(rel_pos_start + (rel_pos_end - rel_pos_start) / 2)

        else:
            last_ex_pos = []
            for line in exon_lookup[gene[REG_ID]]:
                trcpt = line.split(",")
                trcpt = np.asarray([x.split("-") for x in trcpt]).astype(int)
                # at this point, "trcpt" is an array of arrays of start- and end-position for each exon in the transcript
                last_ex_ind = 0
                last_ex = np.zeros(2)
                for i in range(len(trcpt)):
                    if (gene[REG_STR] == "-" and val[0][CNT_ST] == trcpt[i, 0] and val[0][CNT_END] == trcpt[i, 1]) or \
                            (gene[REG_STR] == "+" and val[-1][CNT_ST] == trcpt[i, 0] and val[-1][CNT_END] == trcpt[i, 1]):
                        last_ex_ind = i
                        break

                if gene[REG_STR] == '-':
                    loop_range = range(len(trcpt) - 1, -1, -1)
                    last_end = trcpt[-1][0] - 1
                else:
                    loop_range = range(0, len(trcpt), 1)
                    last_end = trcpt[0][0] - 1

                take_gene = True
                for i in loop_range:
                    if trcpt[i, 0] > last_end:
                        trcpt[i, 1] = trcpt[i, 1] - (trcpt[i, 0] - (last_end + 1))
                        trcpt[i, 0] = last_end + 1
                        last_end = trcpt[i, 1]
                    # if exons overlap:
                    else:
                        print "CAUTION: Exons overlap"
                        last_end = trcpt[i, 1]
                    if i == last_ex_ind:
                        last_ex = trcpt[i]

                if not take_gene:
                    continue

                if gene[REG_STR] == '-':
                    start = trcpt[-1, 0]
                    end = trcpt[0, 1]
                    interval = (end - start) / 1000.0
                    rel_pos_start = 1000 - int(np.floor((last_ex[1] - start) / interval))
                    rel_pos_end = 1000 - int(np.ceil((last_ex[0] - start) / interval))
                else:
                    start = trcpt[0, 0]
                    end = trcpt[-1, 1]
                    interval = (end - start) / 1000.0
                    rel_pos_start = int(np.ceil((last_ex[0] - start) / interval))
                    rel_pos_end = int(np.floor((last_ex[1] - start) / interval))

                last_ex_pos.append(rel_pos_start + (rel_pos_end - rel_pos_start) / 2)
            if max(last_ex_pos) - min(last_ex_pos) > 200:
                excl_count += 1
                continue
            else:
                end_pos.append(min(last_ex_pos) + (max(last_ex_pos) - min(last_ex_pos)) / 2)
    n, bins, patches = plt.hist(end_pos, bins=50, color="forestgreen")
    print str(len(end_pos)) + " genes excluding " + str(excl_count)
    plt.title("Position of last exon (%s)" % file_name)
    plt.savefig("../2018_degradationPaper/ReadCoverage/" + gn_version + "/lastExonPos_%s" % file_name)

    plt.show()

###############################################################################################


def get_tags(tagline):
    """Extract tags from given tagline in a gtf file"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags


def get_overlap_genes(fn_anno):
    """
    Returns a list of gene names which are overlapping
    """

    data = []
    # For all gene entries - extract unique ID
    for l in open(fn_anno, 'r'):
        if l[ANNO_SEQ_NAME] == '#':
            continue  # comments
        l_spl = l.strip('\n').split('\t')
        if l_spl[2].lower() != 'gene':
            continue
        tags = get_tags(l_spl[ANNO_ATTRIBUTE])
        data.append([tags['gene_id'], '%s:%s-%s' % (l_spl[ANNO_SEQ_NAME], l_spl[ANNO_START], l_spl[ANNO_END])])

    data = sp.array(data)

    # fix positions
    pos = data[:, 1]  # example-entry:  chr1:11869-14409
    pos = sp.array([x.split(':')[0] + '-' + x.split(':')[1] for x in pos])  # example-entry:  chr1-11869-14409
    pos = sp.array([x.strip('chr') for x in pos])  # example-entry:  1-11869-14409
    pos = sp.array([x.split('-') for x in pos])  # example-entry:  ['1' '11869' '14409']
    pos[pos[:, 0] == 'X', 0] = '23'
    pos[pos[:, 0] == 'Y', 0] = '24'

    # filter weird things like mitochondria etc.
    i_ok = np.core.defchararray.isdigit(pos[:, 0])  # boolean array of same size
    pos = pos[i_ok, :]
    data = data[i_ok, :]
    pos = pos.astype('int')  # example-entry array([    1, 11869, 14409])

    # sort everything nicely
    sidx = sp.lexsort((pos[:, 2], pos[:, 1], pos[:, 0]))  # last is primary sorting-key
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
    return sp.unique(my_overlap_genes)  # np-array of gene-IDs


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
        if l[ANNO_SEQ_NAME] == '#':
            continue

        l_spl = l.strip('\n').split('\t')

        if l_spl[ANNO_FEATURE].lower() != 'exon':
            continue
        tags = get_tags(l_spl[ANNO_ATTRIBUTE])
        try:
            transcripts[tags['transcript_id']].append('-'.join([l_spl[ANNO_START], l_spl[ANNO_END]]))
        except KeyError:
            transcripts[tags['transcript_id']] = ['-'.join([l_spl[ANNO_START], l_spl[ANNO_END]])]

    # read transcript annotation
    for l in open(fn_anno, 'r'):
        if l[ANNO_SEQ_NAME] == '#':
            continue
        l_spl = l.strip('\n').split('\t')

        if l_spl[ANNO_FEATURE].lower() != 'transcript':
            continue
        tags = get_tags(l_spl[ANNO_ATTRIBUTE])

        key = tags['gene_id']
        gene_type = tags['gene_type']

        # remove overlapping genes
        if key in overlap_genes:
            continue
        # remove non-chr-contigs
        if l_spl[ANNO_SEQ_NAME].lower() not in chr_whitelist:
            continue
        # filter for protein-coding genes only (if filter turned on)
        if protein_coding_filter and (gene_type != "protein_coding"):
            continue

        value = '%s:%s:%s' % (l_spl[ANNO_SEQ_NAME], ','.join(transcripts[tags['transcript_id']]), l_spl[ANNO_STRAND])

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
    exons.sort()

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

    all_exons = my_exons[u_idx]

    # get length of all transcripts
    my_ex_struct_l = []
    for i, rec in enumerate(tcrpt):
        my_ex_struct_l.append(get_transcript_length(rec))

    return [tcrpt[0].split(':')[0], tcrpt[0].split(':')[2], str(sp.median(my_ex_struct_l))], uq_const_ex, all_exons


def read_annotation_file(fn_anno, protein_coding_filter):
    # get np-array of IDs of overlapping genes
    overlap_genes = get_overlap_genes(fn_anno)

    # reading in annotation file
    data = reading_anno(fn_anno, overlap_genes, protein_coding_filter)
    # data is a dictionary with gene IDs as keys and a list of transcripts with format CHR:listOfExons(e1,e2,...):STRAND as values

    exon_lookup = dict()
    for gid in data.keys():
        for t in data[gid]:
            try:
                exon_lookup[gid].append(t.split(":")[1])
            except KeyError:
                exon_lookup[gid] = [t.split(":")[1]]

    uq_g_id = data.keys()  # unique gene ids
    new_data = []

    const_exons = dict()
    all_exons = dict()
    for gid in uq_g_id:
        # process transcripts
        if len(data[gid]) == 1:
            temp, c_e = process_single_transcript_genes(data[gid])
            a_e = c_e
        else:
            temp, c_e, a_e = process_multi_transcript_genes(data[gid])

        # make sure it has been processed correctly
        if temp is None:
            continue
        else:
            const_exons[gid] = c_e
            all_exons[gid] = a_e
            new_data.append([gid] + temp)
    new_data = sp.array(new_data)
    s_idx = sp.argsort(new_data[:, 0])  # filter gene with no name
    new_data = new_data[s_idx, :]

    # new_data is array with entries: gene_ID, seq_name, strand, (median) transcript-length
    return sp.array(new_data), const_exons, all_exons, exon_lookup


def get_annotation_table(fn_anno, protein_coding_filter):

    if fn_anno.lower().endswith('gtf'):
        exon_t_gene, const_exons, all_exons, exon_lookup = read_annotation_file(fn_anno, protein_coding_filter)
    else:
        raise Exception(
            "Only annotation files in format gtf are supported. File name must end accordingly")

    return exon_t_gene, const_exons, all_exons, exon_lookup


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
    if len(regions.shape) <= 1:
        exit("regions.shape was not > 1")
    sidx = sp.argsort(regions[:, 1])  # seq_name

    for i, ii in enumerate(sidx):
        if i > 0 and i % 100 == 0:
            print '%i rounds to go. ETA %.0f seconds' % (regions.shape[0] - i, (time.time() - t0) / i * (regions.shape[0] - i))

        rec = regions[ii]  # e.g. ENSG00000233493.3_2	chr19	-	1000
        rec_exons = exons[rec[0]]  # rec[0] is unique gene ID

        exon_counts = np.zeros((len(rec_exons), 4), dtype=float)  # store start, end, length, count

        if len(regions.shape) == 1:
            exit("regions.shape == 1")
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
    gn_version = "hg38_ffpe_vs_ff"
    options = parse_options(sys.argv)
    if os.path.exists("./" + gn_version + "/anno.tmp") \
            and os.path.exists("./" + gn_version + "/const_ex.pkl") \
            and os.path.exists("./" + gn_version + "/exon_lookup.pkl"):
            # and os.path.exists("./" + gn_version + "/all_ex.pkl") \
        exon_t_gene = sp.loadtxt("./" + gn_version + "/anno.tmp", delimiter='\t', dtype='string')
        const_exons = pickle.load(open("./" + gn_version + "/const_ex.pkl", "rb"))
        # #all_exons = pickle.load(open("./" + gn_version + "/all_ex.pkl", "rb"))
        exon_lookup = pickle.load(open("./" + gn_version + "/exon_lookup.pkl", "rb"))
    else:
        exon_t_gene, const_exons, all_exons, exon_lookup = get_annotation_table(options.fn_anno, options.proteinCodingFilter)
        sp.savetxt("./" + gn_version + "/anno.tmp", exon_t_gene, delimiter='\t', fmt='%s')
        f = open("./" + gn_version + "/const_ex.pkl", "wb")
        pickle.dump(const_exons, f)
        f.close()
        f = open("./" + gn_version + "/all_ex.pkl", "wb")
        pickle.dump(all_exons, f)
        f.close()
        f = open("./" + gn_version + "/exon_lookup.pkl", "wb")
        pickle.dump(exon_lookup, f)
        f.close()

    if os.path.exists("./" + gn_version + "/const_count_data.pkl") \
            and os.path.exists("./" + gn_version + "/file_names.pkl"):
            # and os.path.exists("./" + gn_version + "/all_count_data.pkl"):
        # file_names = ['FFPE_1', 'FFPE_2', 'FFPE_3', 'FFPE_4', 'FF_1', 'FF_2', 'FF_3', 'FF_4']
        file_names = pickle.load(open("./" + gn_version + "/file_names.pkl", "rb"))
        const_data = pickle.load(open("./" + gn_version + "/const_count_data.pkl", "rb"))
        # #all_data = pickle.load(open("./" + gn_version + "/all_count_data.pkl", "rb"))
        # data is a dictionary with unique gene_IDs as keys and
        # a list of (constitutive) exons (start, end, (normalized) count) as value
    else:
        if options.dir_bam != '-':
            file_names = glob.glob(os.path.join(options.dir_bam, '*.bam'))
            const_data = get_counts_from_multiple_bam(file_names, exon_t_gene, const_exons)
            # #all_data = get_counts_from_multiple_bam(file_names, exon_t_gene, all_exons)
        else:
            assert options.fn_bam != '-'
            file_names = [options.fn_bam]
            const_data = get_counts_from_multiple_bam(file_names, exon_t_gene, const_exons)
            # #all_data = get_counts_from_multiple_bam(file_names, exon_t_gene, all_exons)
        f = open("./" + gn_version + "/const_count_data.pkl", "wb")
        pickle.dump(const_data, f)
        f.close()
        # #f = open("./" + gn_version + "/all_count_data.pkl", "wb")
        # #pickle.dump(all_data, f)
        # #f.close()
        f = open("./" + gn_version + "/file_names.pkl", "wb")
        pickle.dump(file_names, f)
        f.close()

    for i in range(len(file_names)):
        f_name = file_names[i].split("/")[-1].strip(".bam")  # FOR FFPE_VS_FF
        # f_name = file_names[i].split("/")[-1].split(".")[0]
        # if f_name in DEG:
        #     f_name = "Deg-" + f_name
        # elif f_name in NON_DEG:
        #     f_name = "Non-" + f_name
        # else:
        #     print "ERROR"
        cut_off = prepare_outlier_filter(const_data[i], exon_t_gene)
        # avg_count_per_exon(const_data[i], exon_t_gene, f_name, cut_off, False, gn_version)
        # distr_over_gene_lengths(const_data[i], exon_t_gene, f_name, cut_off, False, gn_version)
        # last_const_exon_length(const_data[i], exon_t_gene, f_name, gn_version)
        last_const_exon_pos(const_data[i], exon_t_gene, f_name, exon_lookup, gn_version)


if __name__ == "__main__":
    main()
