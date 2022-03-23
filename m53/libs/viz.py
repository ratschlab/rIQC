import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np


def plotBias(vals, fn_plot, myidx, logScale = False, refname = 'TCGA'):

    iqr    = ( (np.percentile(vals[~myidx],75) - np.percentile(vals[~myidx],25) ) * 1.5)
    iqr2    = ( (np.percentile(vals[myidx],75) - np.percentile(vals[myidx],25) ) * 1.5)

    sidx   = np.argsort(vals)
    vals   = vals[sidx]
    myidx = myidx[sidx]

    fig  = plt.figure(figsize=(12,10))
    ax   = fig.add_subplot(111)
    ax_c = ax.twinx()
    ax.vlines(np.array(np.arange(np.sum(vals.shape[0])))[myidx],[0], vals[myidx], label = '%s Reference'%refname)
    ax.vlines(np.array(np.arange(np.sum(vals.shape[0])))[~myidx],[0], vals[~myidx], color = 'r', label = 'Your Samples')

    ax.plot([0,vals.shape[0]],[3,3], '--', color = 'green')
    ax.plot([0,vals.shape[0]],[5,5] , '--',color = 'green')
    ax.plot([0,vals.shape[0]],[iqr + np.percentile(vals[~myidx], 75),iqr + np.percentile(vals[~myidx], 75)], '--',color = 'green')
    ax.plot([0,vals.shape[0]],[iqr2 + np.percentile(vals[myidx], 75),iqr2 + np.percentile(vals[myidx], 75)], '--',color = 'green')

#    ax.plot([0,vals.shape[0]],[6.25,6.25],'--', color = 'green')
    ax.plot([0,vals.shape[0]],[10,10] , '--',color = 'green')
    ax.set_ylabel('Median 3\'/5\' Bias')
    ax.set_xlim(0,vals.shape[0])
    if logScale:
        ax.set_yscale('log')
        ax_c.set_yscale('log')
    ax_c.set_ylim(ax.get_ylim())

    ### add right side ticks
    if logScale:       
        tick_thresholds = np.array([3,5,iqr+np.percentile(vals[~myidx],75),iqr2 + np.percentile(vals[myidx], 75), 10])
    else:
        tick_thresholds = np.array([3,5,iqr+np.percentile(vals[~myidx],75),iqr2 + np.percentile(vals[myidx], 75), 10])
    tick_idx        = np.argsort(tick_thresholds)
    tick_thresholds = tick_thresholds[tick_idx]
    tick_thresholds = np.around(tick_thresholds, decimals = 2)
    ax_c.set_yticks(tick_thresholds)

    tick_thresholds                = tick_thresholds.astype('|S4')
    tick_thresholds                = tick_thresholds.astype('|S50')
    tick_thresholds[tick_idx == 2] = tick_thresholds[tick_idx == 2][0] + ' (Your Filter)'
#    tick_thresholds[tick_idx == 3] = tick_thresholds[tick_idx == 3][0] + ' (PRAD Filter)'
    tick_thresholds[tick_idx == 3] = tick_thresholds[tick_idx == 3][0] + ' (%s Filter)'%(refname)

    ax_c.set_yticklabels(tick_thresholds)


    ax.grid()
    ax.legend(loc=2)
    plt.tight_layout()
    plt.savefig(fn_plot, dpi = 300)
    plt.clf()

