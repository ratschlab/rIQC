import sys
import scipy as sp
import os
import time
import cPickle
import gzip
import numpy.random as npr

npr.seed(23)

REV_DIC = {'A': 'T',
           'C': 'G',
           'G': 'C',
           'T': 'A',
           'N': 'N'}


def __read_genome(fname):
    seq = []
    curr_chrm = ''

    if fname.endswith('.gz'):
        fh = gzip.open(fname, 'r')
    else:
        fh = open(fname, 'r')

    for l, line in enumerate(fh):
        if line[0] == '>':
            if curr_chrm != '':
                yield (curr_chrm, ''.join(seq))
            seq = []
            curr_chrm = line.strip().split(' ')[0].strip('>').strip('chr')
        else:
            seq.append(line.strip().upper())
    if curr_chrm != '':
        yield (curr_chrm, ''.join(seq))
    fh.close()


def __reverse_complement(seq):
    return ''.join([REV_DIC[_] for _ in seq][::-1])


def prepare_kmers(options, regions):
    print 'Preparing genomic kmers'
    cnt = 0
    # MM: creates array of empty sets (one set for each entry in annotation file)
    kmers1 = [set() for _ in regions]
    kmers2 = [set() for _ in regions]
    t0 = time.time()
    # MM: creates array of chromosome names in the order they occur in annotation file (same length as kmers)
    chrms = sp.array([_.strip('chr') for _ in regions[:, 2]])
    # MM: sequence for each chrm in .fasta file
    for chrm, seq in __read_genome(options.fn_genome):
        print 'processing %s' % chrm
        # MM: array of all indices in chrms-array that match chrm name from .fasta
        idx = sp.where(chrms == chrm)[0]
        for i in idx:
            rec = regions[i, :]
            if cnt > 0 and cnt % 100 == 0:
                # MM: regions.shape[0] is number of genes from annotation file
                print '%i rounds to go. ETA %.0f seconds' \
                      % (regions.shape[0] - cnt, (time.time() - t0) / cnt * (regions.shape[0] - cnt))
            cnt += 1

            # MM: TODO: when would that happen?
            if len(regions.shape) == 1:
                start1 = int(rec.split(':')[1].split('-')[0])
                end1 = int(rec.split(':')[1].split('-')[1])
                start2 = None
                end2 = None
            else:
                # MM: first consecutive exon
                start1 = int(rec[0].split(':')[1].split('-')[0])
                end1 = int(rec[0].split(':')[1].split('-')[1])
                # MM: last consecutive exon
                start2 = int(rec[1].split(':')[1].split('-')[0])
                end2 = int(rec[1].split(':')[1].split('-')[1])

            if end1 - start1 > options.k:
                for s in range(start1, end1 - options.k + 1):
                    kmers1[i].add(seq[s:s + options.k])
            if not start2 is None and end2 - start2 > options.k:
                for s in range(start2, end2 - options.k + 1):
                    kmers2[i].add(seq[s:s + options.k])
    # MM: kmers1/2 now consist of sets of kmers of length k, extraced from the fasta file
    #     (kmers that match one entry in annotation are in one set)
    return (kmers1, kmers2)


def clean_kmers(options, kmers1, kmers2):
    if(options.fn_pickle_all != None):
        kmer_pickle = options.fn_pickle_all
    else:
        kmer_pickle = 'all_kmers_k%i.pickle' % options.k

    print 'Making kmers unique'
    if os.path.exists(kmer_pickle):
        (all_kmers1, all_kmers2) = cPickle.load(open(kmer_pickle, 'r'))
    else:
        # MM: creates dictionaries with values: 0 and keys: all existing kmers (from fasta)
        all_kmers1 = dict([[_, 0] for s in kmers1 for _ in s])
        all_kmers2 = dict([[_, 0] for s in kmers2 for _ in s])
        for chrm, seq in __read_genome(options.fn_genome):
            print '\nprocessing %s' % chrm
            for s in range(0, len(seq) - options.k + 1):
                if s > 0 and s % 100000 == 0:
                    sys.stdout.write('.')
                    if s % 1000000 == 0:
                        sys.stdout.write('%i/%i\n' % (s, len(seq)))
                    sys.stdout.flush()
                # MM: only takes kmers that we picked based on annotation file
                try:
                    all_kmers1[seq[s:s + options.k]] += 1
                    all_kmers1[__reverse_complement(seq[s:s + options.k])] += 1
                except KeyError:
                    pass
                try:
                    all_kmers2[seq[s:s + options.k]] += 1
                    all_kmers2[__reverse_complement(seq[s:s + options.k])] += 1
                except KeyError:
                    pass
        cPickle.dump((all_kmers1, all_kmers2), open(kmer_pickle, 'w'), -1)

    ### remove all non-unique entries
    removed = 0
    total = 0
    # MM: rec is a set of kmers
    for i, rec in enumerate(kmers1):
        size_old = len(rec)
        total += size_old
        kmers1[i] = filter(lambda x: all_kmers1[x] == 1 and not x in all_kmers2, rec)
        removed += (size_old - len(kmers1[i]))
    for i, rec in enumerate(kmers2):
        size_old = len(rec)
        total += size_old
        kmers2[i] = filter(lambda x: all_kmers2[x] == 1 and not x in all_kmers1, rec)
        removed += (size_old - len(kmers2[i]))
    if(float(total) != 0):
        print 'Removed %i non-unique kmers (%.2f percent)' % (removed, removed / float(total) * 100)

    if(options.fn_pickle_filt != None):
        cPickle.dump((kmers1, kmers2), open(options.fn_pickle_filt, 'w'), -1)
    else:
        cPickle.dump((kmers1, kmers2), open(('filt_kmers_k%i.pickle' % options.k), 'w'), -1)
    return (kmers1, kmers2)


def get_counts_from_multiple_fastq(fn_fastq, kmers1, kmers2, options):
    """ This is a wrapper to concatenate counts for a given list of fastq
        files"""

    if not options.separate_files:
        return get_counts_from_single_fastq(fn_fastq, kmers1, kmers2, options)[:, sp.newaxis]
    else:
        return sp.hstack([get_counts_from_single_fastq(fn_fastq[i], kmers1, kmers2, options)[:, sp.newaxis] for i in
                          range(len(fn_fastq))])


def get_counts_from_single_fastq(fn_fastqs, kmers1, kmers2, options):
    all_kmers1 = dict([[_, 0] for s in kmers1 for _ in s])
    all_kmers2 = dict([[_, 0] for s in kmers2 for _ in s])

    use_fraction = False
    if options.kmer_thresh <= 1:
        use_fraction = True

    for fn_fastq in fn_fastqs:
        print 'Processing %s' % fn_fastq
        cnt = 0
        cnt1 = 0
        if fn_fastq.endswith('gz'):
            fh = gzip.open(fn_fastq, 'r')
        else:
            fh = open(fn_fastq, 'r')
        for l, line in enumerate(fh):
            if l % 4 != 1:
                continue
            cnt += 1
            if not use_fraction and cnt1 > options.kmer_thresh:
                break
            if cnt % 10000 == 0:
                sys.stdout.write('.')
                if cnt % 100000 == 0:
                    sys.stdout.write(' processed %i reads - %i (%.2f%%) used for quantification\n' % (
                    cnt, cnt1, cnt1 / float(cnt) * 100))
                sys.stdout.flush()
            if use_fraction and npr.random() > options.kmer_thresh:
                continue
            sl = line.strip()
            slr = __reverse_complement(sl)
            for s in range(0, len(sl) - options.k + 1, options.step_k):
                try:
                    all_kmers1[sl[s:s + options.k]] += 1
                    cnt1 += 1
                    break
                except KeyError:
                    pass
                try:
                    all_kmers1[slr[s:s + options.k]] += 1
                    cnt1 += 1
                    break
                except KeyError:
                    pass
                try:
                    all_kmers2[sl[s:s + options.k]] += 1
                    cnt1 += 1
                    break
                except KeyError:
                    pass
                try:
                    all_kmers2[slr[s:s + options.k]] += 1
                    cnt1 += 1
                    break
                except KeyError:
                    pass
        fh.close()

    return sp.array([[sp.sum([all_kmers1[x] for x in kmers1[y]]), sp.sum([all_kmers2[x] for x in kmers2[y]])] for y in
                     range(len(kmers1))], dtype='float').ravel('C')
