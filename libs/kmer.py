import sys
import scipy as sp
import os
import time
import cPickle
import numpy.random as npr
npr.seed(23)

REV_DIC = {'A':'T',
           'C':'G',
           'G':'C',
           'T':'A',
           'N':'N'}

def __read_genome(fname):

    seq = []
    curr_chrm = ''
    for l, line in enumerate(open(fname, 'r')):
        if line[0] == '>':
            if curr_chrm != '':
                yield (curr_chrm, ''.join(seq))
            seq = []
            curr_chrm = line.strip().split(' ')[0][1:]
        else:
            seq.append(line.strip())
    if curr_chrm != '':
        yield (curr_chrm, ''.join(seq))

def __reverse_complement(seq):

    return ''.join([REV_DIC[_] for _ in seq][::-1])
             
def prepare_kmers(options, regions):

    print 'Preparing genomic kmers'
    cnt = 0
    kmers1 = [set() for _ in regions]
    kmers2 = [set() for _ in regions]
    t0 = time.time()
    chrms = sp.array([_.strip('chr') for _ in regions[:, 2]])
    for chrm, seq in __read_genome(options.fn_genome):
        print 'processing %s' % chrm
        idx = sp.where(chrms == chrm)[0]
        for i in idx:
            rec = regions[i, :]
            if cnt > 0 and cnt % 100 == 0:
                print '%i rounds to go. ETA %.0f seconds' % (regions.shape[0] - cnt, (time.time() - t0) / cnt * (regions.shape[0] - cnt))
            cnt += 1

            if len(regions.shape) == 1:
                start1 = int(rec.split(':')[1].split('-')[0])
                end1   = int(rec.split(':')[1].split('-')[1])
                start2 = None
                end2   = None
            else:
                start1 = int(rec[0].split(':')[1].split('-')[0])
                end1   = int(rec[0].split(':')[1].split('-')[1])
                start2 = int(rec[1].split(':')[1].split('-')[0])
                end2   = int(rec[1].split(':')[1].split('-')[1])

            if end1 - start1 > options.k:
                for s in range(start1, end1 - options.k + 1):
                    kmers1[i].add(seq[s:s+options.k])
            if not start2 is None and end2 - start2 > options.k:        
                for s in range(start2, end2 - options.k + 1):
                    kmers2[i].add(seq[s:s+options.k])
    return (kmers1, kmers2)

def clean_kmers(options, kmers1, kmers2):

    kmer_pickle = 'all_kmers_k%i.pickle' % options.k
    print 'Making kmers unique'
    if os.path.exists(kmer_pickle):
        (all_kmers1, all_kmers2) = cPickle.load(open(kmer_pickle, 'r'))
    else:
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
                try:
                    all_kmers1[seq[s:s+options.k]] += 1
                    all_kmers1[__reverse_complement(seq[s:s+options.k])] += 1
                except KeyError:
                    pass
                try:
                    all_kmers2[seq[s:s+options.k]] += 1
                    all_kmers2[__reverse_complement(seq[s:s+options.k])] += 1
                except KeyError:
                    pass
        cPickle.dump((all_kmers1, all_kmers2), open(kmer_pickle, 'w'), -1)

    ### remove all non-unique entries
    removed = 0
    total = 0
    for i, rec in enumerate(kmers1):
        size_old = len(rec)
        total += size_old
        kmers1[i] = filter(lambda x: all_kmers1[x] == 1, rec)
        removed += (size_old - len(kmers1[i]))
    for i, rec in enumerate(kmers2):
        size_old = len(rec)
        total += size_old
        kmers2[i] = filter(lambda x: all_kmers2[x] == 1, rec)
        removed += (size_old - len(kmers2[i]))
    print 'Removed %i non-unique kmers (%.2f percent)' % (removed, removed / float(total) * 100)

    return (kmers1, kmers2)

def get_counts_from_multiple_fastq(fn_fastq, kmers1, kmers2, options):
    """ This is a wrapper to concatenate counts for a given list of fastq
        files"""
    
    if len(fn_fastq) == 1:
        return get_counts_from_single_fastq(fn_fastq[0], kmers1, kmers2, options)[:, sp.newaxis]
    else:
        return sp.hstack([get_counts_from_single_fastq(fn_fastq[i], kmers1, kmers2, options)[:,sp.newaxis] for i in range(len(fn_fastq))])

def get_counts_from_single_fastq(fn_fastq, kmers1, kmers2, options):

    all_kmers1 = dict([[_, 0] for s in kmers1 for _ in s])
    all_kmers2 = dict([[_, 0] for s in kmers2 for _ in s])

    use_fraction = False
    if options.kmer_thresh <= 1:
        use_fraction = True

    print 'Processing %s' % fn_fastq
    cnt = 0
    cnt1 = 0
    for l, line in enumerate(open(fn_fastq, 'r')):
        if l % 4 != 1:
            continue
        if use_fraction and npr.random() > options.kmer_thresh:
            continue
        cnt += 1
        if not use_fraction and cnt1 > options.kmer_thresh:
            break
        if cnt % 10000 == 0:
            sys.stdout.write('.')
            if cnt % 100000 == 0:
                sys.stdout.write(' processed %i reads - %i (%.2f%%) used for quantification\n' % (cnt, cnt1, cnt1 / float(cnt) * 100))
            sys.stdout.flush()
        sl = line.strip()
        slr = __reverse_complement(sl)
        for s in range(0, len(sl) - options.k + 1, options.step_k):
            try:
                all_kmers1[sl[s:s+options.k]] +=1
                cnt1 += 1
                break
            except KeyError:
                pass
            try:
                all_kmers1[slr[s:s+options.k]] +=1
                cnt1 += 1
                break
            except KeyError:
                pass
            try:
                all_kmers2[sl[s:s+options.k]] +=1
                cnt1 += 1
                break
            except KeyError:
                pass
            try:
                all_kmers2[slr[s:s+options.k]] +=1
                cnt1 += 1
                break
            except KeyError:
                pass

    return sp.array([[sp.sum([all_kmers1[x] for x in kmers1[y]]), sp.sum([all_kmers2[x] for x in kmers2[y]])] for y in range(len(kmers1))], dtype='float').ravel('C') 