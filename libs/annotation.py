import pysam
import pdb
import scipy as sp
import scipy.stats as spst
import numpy as np
import time
import usefulTools as ut
import os

# Some Numbers
NMB_CHR = 23

# GFF and GTF column numbers
SEQ_NAME = 0    # name of chromosome or scaffold
SOURCE = 1      # name of program that generated this feature
FEATURE = 2     # feature type name (e.g. "gene", "transcript", "exon")
                #  type of feature (term or accession from SOFA sequence ontology)
START = 3       # start position of feature (seq numbering starting at 1)
END = 4         # end position of feature (seq numbering starting at 1)
SCORE = 5       # a floating point value
STRAND = 6      # + (forward) or - (reverse)
FRAME = 7       # 0/1/2 : position in codon
ATTRIBUTE = 8   # semicolon-separated list of tag-value pairs

def removeNonChrContigs(exonTgene):
    chr_whitelist = [str(x) for x in range(NMB_CHR)]
    chr_whitelist.extend(['chr%i' % i for i in range(NMB_CHR)])
    chr_whitelist.extend(['chrx', 'chry', 'chrm', 'x', 'y', 'm', 'mt'])
    k_idx = sp.array([x.lower() in chr_whitelist for x in exonTgene[:, 2]], dtype='bool')
    
    return exonTgene[k_idx, :]

def filterToInterestingGenes(exonTgene, length):
    t_25 = spst.scoreatpercentile(exonTgene[:, 4].astype('float'), 25)
    t_75 = spst.scoreatpercentile(exonTgene[:, 4].astype('float'), 75)

    if length == 'uq':
        k_idx = sp.where(exonTgene[:, 4].astype('float') > t_75)[0]
    elif length == 'mq':
        k_idx = sp.where((exonTgene[:, 4].astype('float') > t_25) & (exonTgene[:, 4].astype('float') < t_75))[0]
    elif length == 'lq':
        k_idx = sp.where(exonTgene[:, 4].astype('float') < t_25)[0]
    else:
        raise Exception('--length should be one of: uq, mq, lq -- currently is: %s' % length)
    
    return exonTgene[k_idx, :]



def getAnnotationTable(options, lengthFilter=True):
    if options.fn_genes == '-':
        if os.path.exists(options.fn_anno_tmp):
            exonTgene = sp.loadtxt(options.fn_anno_tmp, delimiter='\t', dtype='string')
        else:
            if options.fn_anno.lower().endswith('gff') or options.fn_anno.lower().endswith('gff3'):
                exonTgene = readAnnotationFile(options.fn_anno, options.protein_coding_filter, format='gff')
            elif options.fn_anno.lower().endswith('gtf'):
                exonTgene = readAnnotationFile(options.fn_anno, options.protein_coding_filter, format='gtf')
            else:
                raise Exception(
                    "Only annotation files in formats: gff and gtf are supported. File name must end accordingly")
            sp.savetxt(options.fn_anno_tmp, exonTgene, delimiter='\t', fmt='%s')

        ### remove non chr contigs
        exonTgene = removeNonChrContigs(exonTgene)

        ## filter exonTgene to only retain genes we are interested in (length, splicing, etc)
        if(lengthFilter):
            exonTgene = filterToInterestingGenes(exonTgene, options.length)

    else:
        exonTgene = sp.loadtxt(options.fn_genes, delimiter=' ', dtype='string')

    return exonTgene


def getTranscriptLength(rec):
    expieces = sp.array(rec.split(':')[1].split(',')) # list of all exonintervals
    iEnd = expieces.shape[0]
    lgt = 0

    for i, x in enumerate(expieces[0:]):
        start, end = x.split('-')
        if i == 0:
            lgt += (int(end) - int(start) + 1) / 2.
        elif i == (iEnd - 1):
            lgt += (int(end) - int(start) + 1) / 2.
        else:
            lgt += int(end) - int(start) + 1

    return lgt


def getTranscriptLengthBex(rec, firstEx, lastEx):
    """
    Returns transcript legnth defined as 0.5 first exon
    everything in between and 0.5 last exon
    """
    expieces = sp.array(rec.split(':')[1].split(','))
    foundFirst = False
    lgt = 0
    for i, x in enumerate(expieces):
        start, end = x.split('-')
        if x == firstEx:
            foundFirst = True
            lgt += 0.5 * (int(end) - int(start) + 1)
            continue
        if x == lastEx:
            foundFirst = False
            lgt += 0.5 * (int(end) - int(start) + 1)
            break
        if foundFirst:
            lgt += int(end) - int(start) + 1
    assert lgt != 0, "Outch, there are transcripts with no length"
    return lgt


def getOverlapGenes(fn, format):
    """
    Returns a list of gene names which are overlapping
    """

    data = []
    if format == 'gtf':
        for l in open(fn, 'r'):
            ## comments
            if l[SEQ_NAME] == '#':
                continue
            lSpl = l.strip('\n').split('\t')
            if lSpl[FEATURE].lower() != 'gene':
                continue
            tags = get_tags_gtf(lSpl[ATTRIBUTE])
            data.append([tags['gene_id'], '%s:%s-%s' % (lSpl[SEQ_NAME], lSpl[START], lSpl[END])])
    elif format in ['gff', 'gff3']:
        for l in open(fn, 'r'):
            ## comments
            if l[SEQ_NAME] == '#':
                continue
            lSpl = l.strip('\n').split('\t')
            if not lSpl[FEATURE].lower() in ['gene', 'lincrna_gene', 'mirna_gene', 'processed_transcript', 'rrna_gene',
                                       'snrna_gene', 'snorna_gene']:
                continue
            tags = get_tags_gff3(lSpl[ATTRIBUTE])
            try:
                data.append([tags['ID'], '%s:%s-%s' % (lSpl[SEQ_NAME], lSpl[START], lSpl[END])])
            except KeyError:
                data.append([tags['Parent'], '%s:%s-%s' % (lSpl[SEQ_NAME], lSpl[START], lSpl[END])])

                ### data contains two   o columns: gene_ID, GeneLocus (e.g., chr7:130020290-130027948:+)
    data = sp.array(data)

    ### fix positions
    pos = data[:, 1]
    pos = sp.array([x.split(':')[0] + '-' + x.split(':')[1] for x in pos])
    pos = sp.array([x.strip('chr') for x in pos])
    pos = sp.array([x.split('-') for x in pos])
    pos[pos[:, 0] == 'X', 0] = '23'
    pos[pos[:, 0] == 'Y', 0] = '24'

    ### filter weird things like mitochondria etc.
    iOK = np.core.defchararray.isdigit(pos[:, 0])
    pos = pos[iOK, :]
    data = data[iOK, :]
    pos = pos.astype('int')

    ### sort everything nicely
    sidx = sp.lexsort((pos[:, 2], pos[:, 1], pos[:, 0]))
    pos = pos[sidx, :]
    data = data[sidx, :]

    ### find genes with overlapping annotations
    myOverlapGenes = []
    for i in xrange(pos.shape[0]):
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
        if (sp.sum(iChr & iLBEnd & iUBEnd) == 1) and (sp.sum(iChr & iLBSt & iUBSt) == 1):
            continue

        ### extract IDs of overlapping genes
        overlapgenesSt = data[iChr & iUBSt & iLBSt, 0]
        overlapgenesEnd = data[iChr & iUBEnd & iLBEnd, 0]

        overlapgenesSt = sp.array([x.split('|')[0] for x in overlapgenesSt])
        overlapgenesEnd = sp.array([x.split('|')[0] for x in overlapgenesEnd])

        ### this shoudl actually never happen ...
        if (sp.unique(overlapgenesSt).shape[0] == 1) and (sp.unique(overlapgenesEnd).shape[0] == 1):
            continue
        if sp.unique(overlapgenesSt).shape[0] > 1:
            myOverlapGenes.extend(overlapgenesSt.tolist())
        if sp.unique(overlapgenesEnd).shape[0] > 1:
            myOverlapGenes.extend(overlapgenesEnd.tolist())
    return sp.unique(myOverlapGenes)


def readinganno(fn, overlapgenes, proteinCodingFilter, format):
    """
    Reads in all transcript annotations,
    removes overlapping genes and eventually filters for protein-coding genes on the fly
    """
    data = dict()
    ### collect transcript information
    transcripts = dict()
    for l in open(fn, 'r'):
        if l[SEQ_NAME] == '#':
            continue

        lSpl = l.strip('\n').split('\t')

        if lSpl[FEATURE].lower() != 'exon':
            continue
        if format == 'gtf':
            tags = get_tags_gtf(lSpl[ATTRIBUTE])
            try:
                transcripts[tags['transcript_id']].append('-'.join([lSpl[START], lSpl[END]]))
            except KeyError:
                transcripts[tags['transcript_id']] = ['-'.join([lSpl[START], lSpl[END]])]
        elif format == 'gff':
            tags = get_tags_gff3(lSpl[ATTRIBUTE])
            try:
                transcripts[tags['Parent']].append('-'.join([lSpl[START], lSpl[END]]))
            except KeyError:
                transcripts[tags['Parent']] = ['-'.join([lSpl[START], lSpl[END]])]

    #### read transcript annotation
    for l in open(fn, 'r'):
        if l[SEQ_NAME] == '#':
            continue
        lSpl = l.strip('\n').split('\t')

        if format == 'gtf':
            if lSpl[FEATURE].lower() != 'transcript':
                continue
            tags = get_tags_gtf(lSpl[ATTRIBUTE])

            key = tags['gene_id']
            gene_type = tags['gene_type']
            if key in overlapgenes:
                continue

            if (proteinCodingFilter) and (gene_type != "protein_coding"):
                continue
            value = '%s:%s:%s' % (lSpl[SEQ_NAME], ','.join(transcripts[tags['transcript_id']]), lSpl[STRAND])
        elif format == 'gff':
            if not lSpl[FEATURE].lower() in ['transcript', 'pseudogenic_transcript', 'snrna', 'snorna', 'rrna', 'pseudogene',
                                       'processed_transcript', 'processed_pseudogene', 'lincrna', 'mirna']:
                continue
            tags = get_tags_gff3(lSpl[ATTRIBUTE])

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
            if (proteinCodingFilter) and (gene_type != "protein_coding"):
                continue
            try:
                value = '%s:%s:%s' % (lSpl[SEQ_NAME], ','.join(transcripts[tags['ID']]), lSpl[STRAND])
            except KeyError:
                continue

        try:
            data[key].append(value)
        except KeyError:
            data[key] = [value]

    return data


def processSingleTranscriptGenes(tcrpt):
    assert len(tcrpt) == 1, "Too many transcripts to process"

    ### checking that we have at least two exons
    tcrpt = tcrpt[0]
    if tcrpt.find(',') == -1:
        return None

    ### reformat to somewhat convenient reading
    #format# ID : exon1positions,exon2positions,...,exonNpositions : STRAND

    firstEx = tcrpt.split(':')[0] + ':' + tcrpt.split(':')[1].split(',')[0] + ':' + tcrpt.split(':')[2]
    lastEx = tcrpt.split(':')[0] + ':' + tcrpt.split(':')[1].split(',')[-1] + ':' + tcrpt.split(':')[2]
    return [firstEx, lastEx, tcrpt.split(':')[0], tcrpt.split(':')[2], getTranscriptLength(tcrpt)]


def processMultiTranscriptGenes(tcrpts):
    ### all transcript isoforms have at least two exons
    if sp.sum(np.core.defchararray.find(tcrpts, ',') != -1) != len(tcrpts):
        return None

    #### make matrix of transcript struc and length
    myExons = [x.split(':')[1].split(',') for x in tcrpts]
    myExons = sp.array([reduce(lambda x, y: x + y, myExons)]).ravel()  ### unravel exons into one list of exons
    myExonsInt = sp.array([x.split('-') for x in myExons]).astype('int')

    ### sort this
    sidxInt = sp.lexsort((myExonsInt[:, 1], myExonsInt[:, 0]))
    myExons = myExons[sidxInt]
    myExonsInt = myExonsInt[sidxInt, :]

    ### see how often i got each item
    dummy, uidx, dists = ut.unique_rows(myExonsInt, index=True, counts=True)
    N_match = sp.sum(dists == len(tcrpts))

    if N_match < 3:  ### i want at least 3 constitutive exons
        return None

    ### get constitutitve exons
    iConst = dists == len(tcrpts)
    uqConstEx = myExons[uidx][iConst]

    firstEx = uqConstEx[0]
    lastEx = uqConstEx[-1]

    ## get length of all transcripts
    myExStrucL = []
    for i, rec in enumerate(tcrpts):
        myExStrucL.append(getTranscriptLengthBex(rec, firstEx, lastEx))

    firstEx = tcrpts[0].split(':')[0] + ':' + firstEx + ':' + tcrpts[0].split(':')[2]
    lastEx = tcrpts[0].split(':')[0] + ':' + lastEx + ':' + tcrpts[0].split(':')[2]
    return [firstEx, lastEx, tcrpts[0].split(':')[0], tcrpts[0].split(':')[2], str(sp.median(myExStrucL))]


def readAnnotationFile(fn, proteinCodingFilter, format):
    ### get list of overlapping genes
    overlapgenes = getOverlapGenes(fn, format)

    ### reading in
    data = readinganno(fn, overlapgenes, proteinCodingFilter, format)

    uqgid = data.keys()  ###  unique gene ids
    newdata = []
    for gid in uqgid:
        ### process transcripts
        if len(data[gid]) == 1:
            temp = processSingleTranscriptGenes(data[gid])
        else:
            temp = processMultiTranscriptGenes(data[gid])

        ### make sure it has been processed correctly
        if temp is None:
            continue
        else:
            temp.extend([gid])
            newdata.append(temp)

    newdata = sp.array(newdata)
    sidx = sp.argsort(newdata[:, 5])
    newdata = newdata[sidx, :]
    ### filter gene with no name
    return sp.array(newdata)


def get_tags_gff3(tagline):
    """Extract tags from given tagline in a gff or gff3 file"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.split('=')
        tags[tt[0]] = tt[1]
    return tags


def get_tags_gtf(tagline):
    """Extract tags from given tagline in a gtf file"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags



