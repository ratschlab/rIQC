import pysam
import pdb
import time

import scipy as sp
import scipy.stats as spst
import numpy as np
import usefulTools as ut
import os


### SOME NUMBERS
NMB_CHR = 23

### INDICES FOR GTF
GTF_SEQ_NAME = 0 # name of chromosome or scaffold
GTF_SOURCE = 1 # name of program that generated this feature
GTF_FEATURE = 2 # feature type name (e.g. "gene", "transcript", "exon")
GTF_START = 3 # start position of feature (seq numbering starting at 1)
GTF_END = 4 # end position of feature (seq numbering starting at 1)
GTF_SCORE = 5 # a floating point value
GTF_STRAND = 6 # + (forward) or - (reverse)
GTF_FRAME = 7 # 0/1/2 : position in codon
GTF_ATTRIBUTE = 8 # semicolon-separated list of tag-value pairs

### INDICES FOR GFF
GFF_SEQ_ID = 0 # name of chromosome or scaffold
GFF_SOURCE = 1 # name of program that generated this feature
GFF_TYPE = 2 # type of feature (term or accession from SOFA sequence ontology)
GFF_START = 3 # start position of feature (seq numbering starting at 1)
GFF_END = 4 # end position of feature (seq numbering starting at 1)
GFF_SCORE = 5 # a floating point value
GFF_STRAND = 6 # + (forward) or - (reverse)
GFF_PHASE = 7 # 0/1/2 : position in codon
GFF_ATTRIBUTE = 8 # semicolon-separated list of tag-value pairs

### INDICES FOR GAF
GAF_DB = 0 # database from which identifier in "DB_OBJECT_ID" is drawn
GAF_DB_OBJECT_ID = 1 # unique identifier from database in "DB"
GAF_DB_OBJECT_SYMBOL = 2 # (unique and valid) symbol to which "DB_OBJECT_ID" is matched
GAF_QUALIFIER = 3 # flags that modify interpretation of an annotation
GAF_GO_ID = 4 # GO identifier for term attributed to "DB_OBJECT_ID"
GAF_DB_REFERENCE = 5 # one or more unique identifiers for a single source cited as an authority
    # for the attribution of the "GO_ID" to the "DB_OBJECT_ID"
GAF_EVIDENCE_CODE = 6 # see GO evidence code guide for list of valid evidence codes
GAF_WITH_FROM = 7 # one of DB:gene_symbol ; DB:gene_symbol[allele_symbol] ; DB:gene_id ;
    # DB:protein_name ; DB:sequence_id ; GO:GO_id ; CHEBI:CHEBI_id
GAF_ASPECT = 8 # refers to namespace or ontology to which "GO_ID" belongs
GAF_DB_OBJECT_NAME = 9 # name of gene or gene product
GAF_DB_OBJECT_SYNONYM = 10 # gene symbol [or other text]
GAF_DB_OBJECT_TYPE = 11 # description of the type of gene product being annotated
GAF_TAXON = 12 # taxonomic identifier(s)
GAF_DATE = 13 # date on which annotation was made
GAF_ASSIGNED_BY = 14 # database which made the annotation
GAF_ANNOTATION_EXTENSION = 15 # one of DB:gene_id ; DB:sequence_id ; CHEBI:CHEBI_id ;
    # Cell Type Ontology:CL_id ; GO:GO_id
GAF_GENE_PRODUCT_FORM_ID = 16 # allows the annotation of specific variants of that gene or gene product



def remove_non_chr_contigs(exonTgene):
    """Filter exonTgene to remove contigs that are not chromosomes"""
    chr_whitelist = [str(x) for x in range(NMB_CHR)]
    chr_whitelist.extend(['chr%i' % i for i in range(NMB_CHR)])
    chr_whitelist.extend(['chrx', 'chry', 'chrm', 'x', 'y', 'm', 'mt'])
    k_idx = sp.array([x.lower() in chr_whitelist for x in exonTgene[:, 2]], dtype='bool')
    return exonTgene[k_idx, :]



def filter_genes(exonTgene, options):
    """filter exonTgene to only retain genes we are interested in (length, splicing, etc)"""
    t_25 = spst.scoreatpercentile(exonTgene[:, 4].astype('float'), 25)
    t_75 = spst.scoreatpercentile(exonTgene[:, 4].astype('float'), 75)

    if options.length == 'uq':
        k_idx = sp.where(exonTgene[:, 4].astype('float') > t_75)[0]
    elif options.length == 'mq':
        k_idx = sp.where((exonTgene[:, 4].astype('float') > t_25) & (exonTgene[:, 4].astype('float') < t_75))[0]
    elif options.length == 'lq':
        k_idx = sp.where(exonTgene[:, 4].astype('float') < t_25)[0]
    else:
        raise Exception('--length should be one of: uq, mq, lq -- currently is: %s' % options.length)

    return exonTgene[k_idx, :]



def getAnnotationTable(options):
    print "DEBUG : annotation.py : getAnnotationTable()"

    if options.fn_genes == '-':
        if os.path.exists(options.fn_anno_tmp):
            exonTgene = sp.loadtxt(options.fn_anno_tmp, delimiter = '\t', dtype = 'string') 
        else:
            if options.fn_anno.lower().endswith('gaf'):
                exonTgene = readAnnotationFile(options.fn_anno, format='gaf')
            elif options.fn_anno.lower().endswith('gff') or options.fn_anno.lower().endswith('gff3'):
                exonTgene = readAnnotationFile(options.fn_anno, format='gff')
            elif options.fn_anno.lower().endswith('gtf'):
                exonTgene = readAnnotationFile(options.fn_anno, format='gtf')
            else:
                raise Exception("Only annotation files in formats: gaf, gff and gtf are supported. File name must end accordingly")
            sp.savetxt(options.fn_anno_tmp, exonTgene, delimiter = '\t', fmt = '%s')

        remove_non_chr_contigs(exonTgene)

    else:
        exonTgene = sp.loadtxt(options.fn_genes, delimiter = ' ', dtype = 'string')
    return exonTgene


def getTranscriptLength(rec,iFirst, iLast):
    print "DEBUG : annotation.py : getTranscriptLength()"
    expieces = sp.array(rec.split(':')[1].split(','))
    iEnd     = expieces.shape[0] - iLast
    lgt      = 0

    for i,x in enumerate(expieces[iFirst:]):
        start, end = x.split('-')
        if i == 0:
            lgt += ( int(end) - int(start) + 1) / 2.
        elif i == iLast - 1:
            lgt += ( int(end) - int(start) + 1) / 2.
        else:
            lgt += int(end) - int(start) + 1
    return lgt

def getTranscriptLengthBex(rec, firstEx, lastEx):
    """
    Returns transcript length defined as 0.5 first exon
    everything in between and 0.5 last exon
    """
    print "DEBUG : annotation.py : getTranscriptLengthBex()"
    expieces = sp.array(rec.split(':')[1].split(','))
    foundFirst = False
    lgt = 0
    for i,x in enumerate(expieces):
        start, end = x.split('-')
        if x == firstEx:
            foundFirst = True
            lgt += 0.5 * ( int(end) - int(start) + 1)
            continue
        if x == lastEx:
            foundFirst = False
            lgt += 0.5 * ( int(end) - int(start) + 1)
            break
        if foundFirst:
            lgt += int(end) - int(start) + 1
    assert lgt != 0, "Outch, there are transcripts with no length"
    return lgt







def readinganno(fn, overlapgenes, format='gaf'):
    """
    Reads in all transcript annotations
    and removes overlapping genes on the fly
    and keeps only protein-coding genes
    """
    print "DEBUG : annotation.py : readinganno()"

    data = dict()
    ### collect transcript information for certain formats
    if format in ['gtf', 'gff', 'gff3']:
        transcripts = dict()
        for l in open(fn, 'r'):
            if l[0] == '#':
                continue
            lSpl = l.strip('\n').split('\t')
            if lSpl[GTF_FEATURE].lower() != 'exon':
                continue
            if format == 'gtf':
                tags = get_tag_value_pairs_gtf(lSpl[GTF_ATTRIBUTE])
                try:
                    transcripts[tags['transcript_id']].append('-'.join([lSpl[GTF_START], lSpl[GTF_END]]))
                except KeyError:
                    transcripts[tags['transcript_id']] = ['-'.join([lSpl[GTF_START], lSpl[GTF_END]])]
            else:
                tags = get_tag_value_pairs_gff3(lSpl[GFF_ATTRIBUTE])
                try:
                    transcripts[tags['Parent']].append('-'.join([lSpl[GFF_START], lSpl[GFF_END]]))
                except KeyError:
                    transcripts[tags['Parent']] = ['-'.join([lSpl[GFF_START], lSpl[GFF_END]])]
        
    #### read transcript annotation
    for l in open(fn, 'r'):
        if l[0] == '#':
            continue
        lSpl = l.strip('\n').split('\t')
        if format == 'gaf':
            if lSpl[GAF_DB_OBJECT_SYMBOL] != 'transcript':
                continue
            if lSpl[GAF_ASPECT] != 'genome':
                continue
            if lSpl[GAF_ANNOTATION_EXTENSION] == '':
                continue
            if lSpl[GAF_ANNOTATION_EXTENSION].split('|')[0] in overlapgenes:
                continue
            key = lSpl[GAF_ANNOTATION_EXTENSION]
            value = lSpl[GAF_ASSIGNED_BY]
        elif format == 'gtf':
            if lSpl[GTF_FEATURE] != 'transcript':
                continue
            tags = get_tag_value_pairs_gtf(lSpl[GTF_ATTRIBUTE])
            key = tags['gene_id']
            type = tags['gene_type']
            if key in overlapgenes:
                continue
            if type != "protein_coding":
                continue
            value = '%s:%s:%s' % (lSpl[GTF_SEQ_NAME],
                                  ','.join(transcripts[tags['transcript_id']]), lSpl[GTF_STRAND])
        elif format in ['gff', 'gff3']:
            if not lSpl[2].lower() in ['transcript', 'pseudogenic_transcript', 'snrna', 'snorna', 'rrna', 'pseudogene', 'processed_transcript', 'processed_pseudogene', 'lincrna', 'mirna']:
                continue
            tags = get_tag_value_pairs_gff3(lSpl[GFF_ATTRIBUTE])
            try:
                key = tags['Parent']
            except KeyError:
                try:
                    key = tags['geneID']
                except KeyError:
                    continue
            if key in overlapgenes:
                continue
            try:
                value = '%s:%s:%s' % (lSpl[GFF_SEQ_ID], ','.join(transcripts[tags['ID']]), lSpl[GFF_STRAND])
            except KeyError:
                continue
        try:
            data[key].append(value)
        except KeyError:
            data[key] = [value]
    return data


def processSingleTranscriptGenes(tcrpt):
    print "DEBUG : annotation.py : processSingleTranscriptGenes()"

    assert len(tcrpt) == 1, "Too many transcripts to process"

    ### checking that we have at least two exons
    tcrpt = tcrpt[0] 
    if tcrpt.find(',') == -1:
        return None

    ### reformat to somewhat convenient reading
    firstEx = tcrpt.split(':')[1].split(',')[0]
    lastEx  = tcrpt.split(':')[1].split(',')[-1]

    firstEx = tcrpt.split(':')[0]+':' + firstEx + ':' + tcrpt.split(':')[2]
    lastEx  = tcrpt.split(':')[0]+':' + lastEx + ':'  + tcrpt.split(':')[2]
    return [firstEx, lastEx, tcrpt.split(':')[0], tcrpt.split(':')[2], getTranscriptLength(tcrpt, 0, 0)]

def processMultiTranscriptGenes(tcrpts):
    print "DEBUG : annotation.py : processMultiTranscriptGenes()"

    ### all transcript isoforms have at least two exons
    if sp.sum(np.core.defchararray.find(tcrpts,',') != -1) != len(tcrpts):
        return None

    #### make matrix of transcript struc and length
    myExons       = [x.split(':')[1].split(',') for x in tcrpts]
    myExons       = sp.array([reduce(lambda x, y: x + y, myExons)]).ravel() ### unravel exons into one list of exons
    myExonsInt    = sp.array([x.split('-') for x in myExons]).astype('int')

    ### sort this
    sidxInt       = sp.lexsort((myExonsInt[:,1], myExonsInt[:,0]))
    myExons       = myExons[sidxInt]
    myExonsInt    = myExonsInt[sidxInt,:]

    ### see how often i got each item
    dummy, uidx, dists = ut.unique_rows(myExonsInt, index=True, counts = True)
    N_match            = sp.sum(dists == len(tcrpts))

    if N_match < 3: ### i want at lest 3 constitutive exons
        return None

    ### get constitutitve exons
    iConst    = dists == len(tcrpts)
    uqConstEx = myExons[uidx][iConst]

    firstEx   = uqConstEx[0]
    lastEx    = uqConstEx[-1]

    ## get length of all transcripts
    myExStrucL = []
    for i,rec in enumerate(tcrpts):
        myExStrucL.append(getTranscriptLengthBex(rec, firstEx, lastEx))        

    firstEx   = tcrpts[0].split(':')[0] + ':' + firstEx + ':' + tcrpts[0].split(':')[2]
    lastEx    = tcrpts[0].split(':')[0] + ':' + lastEx  + ':' + tcrpts[0].split(':')[2]
    return [firstEx, lastEx, tcrpts[0].split(':')[0],tcrpts[0].split(':')[2], str(sp.median(myExStrucL))]


def getOverlapGenes(fn, format):
    """
    Returns a list of gene names which are overlapping TODO:What does overlapping mean exactly?
    """

    print "DEBUG : annotation.py : getOverlapGenes()"
    ### Read gene annotation
    data = []
    if format == 'gaf':
        for l in open(fn, 'r'):
            lSpl = l.strip('\n').split('\t')
            if lSpl[GAF_DB_OBJECT_SYMBOL] != 'gene':
                continue
            if lSpl[GAF_ASPECT] != 'genome':
                continue
            if lSpl[GAF_ANNOTATION_EXTENSION] == '':
                continue
            data.append([lSpl[GAF_DB_OBJECT_ID], lSpl[GAF_GENE_PRODUCT_FORM_ID]])
    elif format == 'gtf':
        for l in open(fn, 'r'):
            ## comments
            if l[0] == '#':
                continue
            lSpl = l.strip('\n').split('\t')
            if lSpl[GTF_FEATURE].lower() != 'gene':
                continue
            tags = get_tag_value_pairs_gtf(lSpl[GTF_ATTRIBUTE])
            data.append([tags['gene_id'], '%s:%s-%s' %
                         (lSpl[GTF_SEQ_NAME], lSpl[GTF_START], lSpl[GTF_END])])  # 0:seqname, 3:startpos, 4:endpos
    elif format in ['gff', 'gff3']:
        for l in open(fn, 'r'):
            ## comments
            if l[0] == '#':
                continue
            lSpl = l.strip('\n').split('\t')
            if not lSpl[GFF_TYPE].lower() in ['gene', 'lincrna_gene', 'mirna_gene', 'processed_transcript', 'rrna_gene',
                                       'snrna_gene', 'snorna_gene']:
                continue
            tags = get_tag_value_pairs_gff3(lSpl[GFF_ATTRIBUTE])
            try:
                data.append([tags['ID'], '%s:%s-%s' %
                             (lSpl[GFF_SEQ_ID], lSpl[GFF_START], lSpl[GFF_END])])
            except KeyError:
                data.append([tags['Parent'], '%s:%s-%s' %
                             (lSpl[GFF_SEQ_ID], lSpl[GFF_START], lSpl[GFF_END])])

    # data contains two columns: gene_ID, gene_locus (e.g. chr7:140020290-130027948)
    data = sp.array(data)

    ### fix positions
    pos = data[:, 1]
    pos = sp.array([x.split(':')[0] + '-' + x.split(':')[1] for x in pos])
    pos = sp.array([x.strip('chr') for x in pos])
    pos = sp.array([x.split('-') for x in pos])
    pos[pos[:, 0] == 'X', 0] = '23'
    pos[pos[:, 0] == 'Y', 0] = '24'

    ### filter weird things like mitochondria etc.
    # (everything that has no digit as seqname)
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
        mypos = pos[i, :]  # MM mypos = entire row

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

def readAnnotationFile(fn, format='gaf'):
    print "DEBUG : annotation.py : readAnnotationFile()"

    overlapgenes = getOverlapGenes(fn, format)

    ### reading in
    data   = readinganno(fn, overlapgenes, format)

    uqgid   = data.keys() ###  unique gene ids
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
    sidx    = sp.argsort(newdata[:,5])
    newdata = newdata[sidx,:]
    ### filter gene with no name
    return sp.array(newdata)

def get_tag_value_pairs_gff3(tagline):
    """Extract tag-value pairs from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.split('=')
        tags[tt[0]] = tt[1]
    return tags

def get_tag_value_pairs_gtf(tagline):
    """Extract tag-value pairs from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags



