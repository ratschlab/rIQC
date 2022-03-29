import pysam
import scipy.stats as spst
import numpy as np
import time
from . import utils as ut
import os
from functools import reduce

# Some Numbers
NMB_CHR = 23

# GFF and GTF column numbers
SEQ_NAME = 0    # name of chromosome or scaffold
SOURCE = 1      # name of program that generated this feature
FEATURE = 2     # feature type name (e.g. "gene", "transcript", "exon")
START = 3       # start position of feature (seq numbering starting at 1)
END = 4         # end position of feature (seq numbering starting at 1)
SCORE = 5       # a floating point value
STRAND = 6      # + (forward) or - (reverse)
FRAME = 7       # 0/1/2 : position in codon
ATTRIBUTE = 8   # semicolon-separated list of tag-value pairs


def __filter_non_chr_contigs(exon_t_gene):
    chr_whitelist = [str(x) for x in range(NMB_CHR)]
    chr_whitelist.extend(['chr%i' % i for i in range(NMB_CHR)])
    chr_whitelist.extend(['chrx', 'chry', 'chrm', 'x', 'y', 'm', 'mt'])
    k_idx = np.array([x.lower() in chr_whitelist for x in exon_t_gene[:, 2]], dtype='bool')
    
    return exon_t_gene[k_idx, :]


def __filter_gene_length(exon_t_gene, length):
    t_25 = spst.scoreatpercentile(exon_t_gene[:, 4].astype('float'), 25)
    t_75 = spst.scoreatpercentile(exon_t_gene[:, 4].astype('float'), 75)

    if length == 'uq':
        k_idx = np.where(exon_t_gene[:, 4].astype('float') > t_75)[0]
    elif length == 'mq':
        k_idx = np.where((exon_t_gene[:, 4].astype('float') > t_25) & (exon_t_gene[:, 4].astype('float') < t_75))[0]
    elif length == 'lq':
        k_idx = np.where(exon_t_gene[:, 4].astype('float') < t_25)[0]
    else:
        raise Exception('--length should be one of: uq, mq, lq -- currently is: %s' % length)
    
    return exon_t_gene[k_idx, :]


def __name_anno_tmp_file(proteinCodingFilter, legacy):
    fn_anno_tmp = 'anno_tmp'

    if proteinCodingFilter:
        fn_anno_tmp = fn_anno_tmp + '_pcfON'
    else:
        fn_anno_tmp = fn_anno_tmp + '_pcfOFF'

    if legacy:
        fn_anno_tmp = fn_anno_tmp + '_legacy'

    return fn_anno_tmp + '.tsv'


def get_annotation_table(fn_genes, fn_anno_tmp, fn_anno, proteinCodingFilter, lengthFilter, length, legacy):
    if fn_genes == '-':
        if fn_anno_tmp == '':
            fn_anno_tmp = __name_anno_tmp_file(proteinCodingFilter, legacy)
        if os.path.exists(fn_anno_tmp):
            exon_t_gene = np.loadtxt(fn_anno_tmp, delimiter='\t', dtype='str')
        else:
            if fn_anno.lower().endswith('gff') or fn_anno.lower().endswith('gff3'):
                exon_t_gene = __read_annotation_file(fn_anno, proteinCodingFilter, file_format='gff', legacy=legacy)
            elif fn_anno.lower().endswith('gtf'):
                exon_t_gene = __read_annotation_file(fn_anno, proteinCodingFilter, file_format='gtf', legacy=legacy)
            else:
                raise Exception(
                    "Only annotation files in formats: gff and gtf are supported. File name must end accordingly")
            # the temporary anno_*.tsv is saved without being filtered for "interesting" genes
            np.savetxt(fn_anno_tmp, exon_t_gene, delimiter='\t', fmt='%s')

        # Filtering
        exon_t_gene = __filter_non_chr_contigs(exon_t_gene)
        if lengthFilter:
            exon_t_gene = __filter_gene_length(exon_t_gene, length)

    else:
        exon_t_gene = np.loadtxt(fn_genes, delimiter=' ', dtype='str')

    return exon_t_gene


def __get_transcript_length(rec, legacy=False):
    # list of all exon-intervals
    ex_pieces = np.array(rec.split(':')[1].split(','))
    lgt = 0

    for i, x in enumerate(ex_pieces[0:]):
        start, end = x.split('-')
        if legacy and i == 0:
            lgt += (int(end) - int(start) + 1) / 2.
        else:
            lgt += int(end) - int(start) + 1
    return lgt


def __get_transcript_length_bex(rec, firstEx, lastEx, legacy=False):
    """
    Returns transcript length defined as sum of length of exons
    """
    expieces = np.array(rec.split(':')[1].split(','))
    foundFirst = False
    lgt = 0
    for i, x in enumerate(expieces):
        start, end = x.split('-')
        if legacy and x == firstEx:
            foundFirst = True
            lgt += 0.5 * (int(end) - int(start) + 1)
            continue
        if legacy and x == lastEx:
            foundFirst = False
            lgt += 0.5 * (int(end) - int(start) + 1)
            break
        if foundFirst or not legacy:
            lgt += int(end) - int(start) + 1
    assert lgt != 0, "Outch, there are transcripts with no length"
    return lgt


def __get_overlap_genes(fn, format):
    """
    Returns a list of gene names which are overlapping
    """

    data = []
    if format == 'gtf':
        for l in open(fn, 'r'):
            # comments
            if l[SEQ_NAME] == '#':
                continue
            lSpl = l.strip('\n').split('\t')
            if lSpl[2].lower() != 'gene':
                continue
            tags = __get_tags_gtf(lSpl[ATTRIBUTE])
            data.append([tags['gene_id'], '%s:%s-%s' % (lSpl[SEQ_NAME], lSpl[START], lSpl[END])])
    elif format in ['gff', 'gff3']:
        for l in open(fn, 'r'):
            # comments
            if l[SEQ_NAME] == '#':
                continue
            lSpl = l.strip('\n').split('\t')
            if not lSpl[2].lower() in ['gene', 'lincrna_gene', 'mirna_gene', 'processed_transcript', 'rrna_gene',
                                       'snrna_gene', 'snorna_gene']:
                continue
            tags = __get_tags_gff(lSpl[ATTRIBUTE])
            try:
                data.append([tags['ID'], '%s:%s-%s' % (lSpl[SEQ_NAME], lSpl[START], lSpl[END])])
            except KeyError:
                data.append([tags['Parent'], '%s:%s-%s' % (lSpl[SEQ_NAME], lSpl[START], lSpl[END])])

    # data contains two columns: gene_ID, GeneLocus (e.g., chr7:130020290-130027948:+)
    data = np.array(data)

    # fix positions
    pos = data[:, 1]
    pos = np.array([x.split(':')[0] + '-' + x.split(':')[1] for x in pos])
    pos = np.array([x.strip('chr') for x in pos])
    pos = np.array([x.split('-') for x in pos])
    pos[pos[:, 0] == 'X', 0] = '23'
    pos[pos[:, 0] == 'Y', 0] = '24'

    # filter weird things like mitochondria etc.
    iOK = np.core.defchararray.isdigit(pos[:, 0])
    pos = pos[iOK, :]
    data = data[iOK, :]
    pos = pos.astype('int')

    ### sort everything nicely
    sidx = np.lexsort((pos[:, 2], pos[:, 1], pos[:, 0]))
    pos = pos[sidx, :]
    data = data[sidx, :]

    ### find genes with overlapping annotations
    myOverlapGenes = []
    for i in range(pos.shape[0]):
        mypos = pos[i, :]

        ## same chr
        iChr = mypos[0] == pos[:, 0]

        ## end is in something else
        iLBEnd = mypos[2] >= pos[:, 1]
        iUBEnd = mypos[2] <= pos[:, 2]

        ## st is in something else
        iLBSt = mypos[1] <= pos[:, 2]
        iUBSt = mypos[1] >= pos[:, 1]

        ## on both ends the only entry that overlaps to i is i itself --> continue
        if (np.sum(iChr & iLBEnd & iUBEnd) == 1) and (np.sum(iChr & iLBSt & iUBSt) == 1):
            continue

        ### extract IDs of overlapping genes
        overlapgenesSt = data[iChr & iUBSt & iLBSt, 0]
        overlapgenesEnd = data[iChr & iUBEnd & iLBEnd, 0]

        overlapgenesSt = np.array([x.split('|')[0] for x in overlapgenesSt])
        overlapgenesEnd = np.array([x.split('|')[0] for x in overlapgenesEnd])

        ### this shoudl actually never happen ...
        if (np.unique(overlapgenesSt).shape[0] == 1) and (np.unique(overlapgenesEnd).shape[0] == 1):
            continue
        if np.unique(overlapgenesSt).shape[0] > 1:
            myOverlapGenes.extend(overlapgenesSt.tolist())
        if np.unique(overlapgenesEnd).shape[0] > 1:
            myOverlapGenes.extend(overlapgenesEnd.tolist())
    return np.unique(myOverlapGenes)


def __reading_anno(fn, overlapgenes, protein_coding_filter, format):
    """
    Reads in all transcript annotations,
    removes overlapping genes and eventually filters for protein-coding genes on the fly
    """
    data = dict()
    # collect transcript information
    transcripts = dict()
    for l in open(fn, 'r'):
        if l[SEQ_NAME] == '#':
            continue

        l_spl = l.strip('\n').split('\t')

        if l_spl[FEATURE].lower() != 'exon':
            continue
        if format == 'gtf':
            tags = __get_tags_gtf(l_spl[ATTRIBUTE])
            try:
                transcripts[tags['transcript_id']].append('-'.join([l_spl[START], l_spl[END]]))
            except KeyError:
                transcripts[tags['transcript_id']] = ['-'.join([l_spl[START], l_spl[END]])]
        elif format == 'gff':
            tags = __get_tags_gff(l_spl[ATTRIBUTE])
            try:
                transcripts[tags['Parent']].append('-'.join([l_spl[START], l_spl[END]]))
            except KeyError:
                transcripts[tags['Parent']] = ['-'.join([l_spl[START], l_spl[END]])]

    # read transcript annotation
    for l in open(fn, 'r'):
        if l[SEQ_NAME] == '#':
            continue
        l_spl = l.strip('\n').split('\t')

        if format == 'gtf':
            if l_spl[FEATURE].lower() != 'transcript':
                continue
            tags = __get_tags_gtf(l_spl[ATTRIBUTE])

            key = tags['gene_id']
            gene_type = tags['gene_type']
            if key in overlapgenes:
                continue
            if protein_coding_filter and (gene_type != "protein_coding"):
                continue
            value = '%s:%s:%s' % (l_spl[SEQ_NAME], ','.join(transcripts[tags['transcript_id']]), l_spl[STRAND])
        elif format == 'gff':
            if not l_spl[FEATURE].lower() in \
                   ['transcript', 'pseudogenic_transcript', 'snrna', 'snorna', 'rrna', 'pseudogene',
                    'processed_transcript', 'processed_pseudogene', 'lincrna', 'mirna']:
                continue
            tags = __get_tags_gff(l_spl[ATTRIBUTE])

            gene_type = tags['gene_type']

            key_flag = True
            try:
                key = tags['Parent']
            except KeyError:
                try:
                    key = tags['geneID']
                except KeyError:
                    key_flag = False

            if key_flag and key in overlapgenes:
                continue
            if protein_coding_filter and (gene_type != "protein_coding"):
                continue
            try:
                value = '%s:%s:%s' % (l_spl[SEQ_NAME], ','.join(transcripts[tags['ID']]), l_spl[STRAND])
            except KeyError:
                continue

        try:
            data[key].append(value)
        except KeyError:
            data[key] = [value]

    return data


def __process_single_transcript_genes(tcrpt, legacy=False):
    assert len(tcrpt) == 1, "Too many transcripts to process"

    # checking that we have at least two exons
    tcrpt = tcrpt[0]
    if tcrpt.find(',') == -1:
        return None

    #format# ID : exon1positions,exon2positions,...,exonNpositions : STRAND

    firstEx = tcrpt.split(':')[0] + ':' + tcrpt.split(':')[1].split(',')[0] + ':' + tcrpt.split(':')[2]
    lastEx = tcrpt.split(':')[0] + ':' + tcrpt.split(':')[1].split(',')[-1] + ':' + tcrpt.split(':')[2]

    return [firstEx, lastEx, tcrpt.split(':')[0], tcrpt.split(':')[2], __get_transcript_length(tcrpt, legacy)]


def __process_multi_transcript_genes(tcrpts, legacy=False):
    # CAVEAT: We only use transcript isoforms that have at least two exons
    if np.sum(np.core.defchararray.find(tcrpts, ',') != -1) != len(tcrpts):
        return None

    # make matrix of transcript struct and length
    myExons = [x.split(':')[1].split(',') for x in tcrpts]
    # unravel exons into one list of exons
    myExons = np.array([reduce(lambda x, y: x + y, myExons)]).ravel()
    myExonsInt = np.array([x.split('-') for x in myExons]).astype('int')

    # sort this
    sidxInt = np.lexsort((myExonsInt[:, 1], myExonsInt[:, 0]))
    myExons = myExons[sidxInt]
    myExonsInt = myExonsInt[sidxInt, :]

    # see how often we got each item
    dummy, uidx, dists = ut.unique_rows(myExonsInt, index=True, counts=True)
    N_match = np.sum(dists == len(tcrpts))

    # make sure we have at least 3 constitutive exons
    if N_match < 3:
        return None

    # get constitutitve exons
    iConst = dists == len(tcrpts)
    uqConstEx = myExons[uidx][iConst]

    firstEx = uqConstEx[0]
    lastEx = uqConstEx[-1]

    # get length of all transcripts
    myExStrucL = []
    for i, rec in enumerate(tcrpts):
        myExStrucL.append(__get_transcript_length_bex(rec, firstEx, lastEx, legacy))

    firstEx = tcrpts[0].split(':')[0] + ':' + firstEx + ':' + tcrpts[0].split(':')[2]
    lastEx = tcrpts[0].split(':')[0] + ':' + lastEx + ':' + tcrpts[0].split(':')[2]

    return [firstEx, lastEx, tcrpts[0].split(':')[0], tcrpts[0].split(':')[2], str(np.median(myExStrucL))]


def __read_annotation_file(fn, protein_coding_filter, file_format, legacy=False):
    # get list of overlapping genes
    overlapgenes = __get_overlap_genes(fn, file_format)

    # reading annotation file in
    data = __reading_anno(fn, overlapgenes, protein_coding_filter, file_format)

    new_data = []
    for gid in data:  # iterate over unique gene ids
        # process transcripts
        if len(data[gid]) == 1:
            temp = __process_single_transcript_genes(data[gid], legacy)
        else:
            temp = __process_multi_transcript_genes(data[gid], legacy)

        # make sure it has been processed correctly
        if temp is None:
            continue
        else:
            temp.extend([gid])
            new_data.append(temp)
    new_data = np.array(new_data)
    s_idx = np.argsort(new_data[:, 5])
    new_data = new_data[s_idx, :]
    # filter gene with no name
    return np.array(new_data)


def __get_tags_gff(tagline):
    """Extract tags from given tagline in a gff or gff3 file"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.split('=')
        tags[tt[0]] = tt[1]
    return tags


def __get_tags_gtf(tagline):
    """Extract tags from given tagline in a gtf file"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags



