import sys
import scipy as sp
import os
import scipy.stats as spst

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASEDIR = '/cluster/work/grlab/projects/m53'
PLOTDIR = os.path.join(BASEDIR, 'plots')
if not os.path.exists(PLOTDIR):
    os.makedirs(PLOTDIR)

### load bam data
data_bam = sp.loadtxt(os.path.join(BASEDIR, 'bam_based.tsv'), dtype='str', delimiter='\t')
data_bam[:, 0] = sp.array([x.split('/')[-1].split('.')[0] for x in data_bam[:, 0]])

Ks = ['0.1', '0.3', '0.5', '10000', '25000', '50000']

for K in Ks:

    print('processing %s' % K)

    data_fastq = sp.loadtxt(os.path.join(BASEDIR, 'fastq_based_%s.tsv' % K), dtype='str', delimiter='\t')
    data_fastq[:, 0] = sp.array([x.split('/')[-1].split('_')[0] for x in data_fastq[:, 0]])

    kidx = sp.in1d(data_fastq[:, 0], data_bam[:, 0])
    assert sp.all(data_fastq[kidx, 0] == data_bam[:, 0])
    data_fastq = data_fastq[kidx, :]

    _pr,_pp = spst.pearsonr(data_bam[:, 1].astype('float'), data_fastq[:, 1].astype('float'))
    _sr,_sp = spst.spearmanr(data_bam[:, 1].astype('float'), data_fastq[:, 1].astype('float'))

    fig = plt.figure(figsize=(10, 10), dpi=200)
    ax = fig.add_subplot(111)
    ax.plot(data_bam[:, 1].astype('float'), data_fastq[:, 1].astype('float'), 'bo')
    ax.set_title('fastq vs bam - %s' % K)
    ax.set_xlabel('bam')
    ax.set_ylabel('fastq (%s)' % K)
    xl = ax.get_xlim()
    yl = ax.get_ylim()
    ax.text(xl[0] + 0.1*(xl[1] - xl[0]), yl[0] + 0.9*(yl[1] - yl[0]), 'Pearson r: %.2f\nSpearman r: %.2f' % (_pr, _sr))

    plt.savefig(os.path.join(PLOTDIR, 'fastq_vs_bam_%s.png' % K), format='png', bbox_inches='tight')
    plt.close(fig)
