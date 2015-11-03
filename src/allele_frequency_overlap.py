# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/12/14
content:    Plot allele frequencies in overlap of consecutive fragments.
'''
# Modules
import os
import numpy as np

from hivevo.hivevo.patients import Patient
from filenames import get_figure_folder
from hivwholeseq.controls.check_allele_frequency_overlap import get_allele_frequency_overlap
from util import store_data, load_data, fig_width, fig_fontsize



# Globals
overlaps = ['F1-2', 'F2-3', 'F3-4', 'F4-5', 'F5-6']
cov_min = 1000
qual_min = 30




# Functions
def get_map_overlap(sample, fr1, fr2):
    '''Get a coordinate map of the overlap between the two fragments'''
    import numpy as np
    from seqanpy import align_ladder

    seq1 = sample.get_reference(fr1)
    seq2 = sample.get_reference(fr2)
    (score, ali1, ali2) = align_ladder(seq1, seq2, score_gapopen=-20)
    start2 = len(ali2) - len(ali2.lstrip('-'))
    end1 = len(ali1.rstrip('-'))
    
    mapco = []
    pos1 = start2
    pos2 = 0
    for i in xrange(start2, end1):
        if (ali1[i] != '-') and (ali2[i] != '-'):
            mapco.append((pos1, pos2))

        if ali1[i] != '-':
            pos1 += 1

        if ali2[i] != '-':
            pos2 += 1

    return np.array(mapco, int)


def get_allele_frequency_overlap(sample, overlaps, cov_min=1000,
                                 VERBOSE=0, qual_min=30):
    '''Get allele frequency in the overlaps'''
    data = [] 
    for io, overlap in enumerate(overlaps):
        fr1 = overlap[:2]
        fr2 = 'F'+overlap[-1]

        if VERBOSE >= 1:
            print overlap, samplename

        # FIXME: actually use frequencies
        try:
            ac1 = sample.get_allele_counts(fr1, qual_min=qual_min)
            ac2 = sample.get_allele_counts(fr2, qual_min=qual_min)
        except IOError:
            continue

        coord_map = get_map_overlap(sample, fr1, fr2)

        acjoint = np.zeros((2, ac1.shape[0], coord_map.shape[0]), int)
        acjoint[0] = ac1[:, coord_map[:, 0]]
        acjoint[1] = ac2[:, coord_map[:, 1]]

        # Convert to frequencies
        afjoint = np.ma.array(acjoint)
        cov = acjoint.sum(axis=1)
        mask = (cov < cov_min).any(axis=0)
        mask = np.tile(mask, (afjoint.shape[0], afjoint.shape[1], 1))
        afjoint.mask = mask
        afjoint = (1.0 * afjoint.swapaxes(0, 1) / afjoint.sum(axis=1)).swapaxes(0, 1)

        data.append({'af': afjoint,
                     'samplename': samplename,
                     'overlap': overlap,
                     'io': io,
                     'n_templates': sample.get_n_templates_dilutions(),
                     'coverage': cov})

    return data


def estimate_templates_overlaps(sample, data):
    '''Estimate templates for the overlaps'''
    for datum in data:
        af1, af2 = datum['af']

        # Filter only polymorphic sites
        afmin = 3e-3
        indfm = (af1 >= afmin) & (af1 <= 1 - afmin) & (af2 >= afmin) & (af2 <= 1 - afmin)
        nsites = len(np.unique(indfm.nonzero()[1]))

        # Estimate the template number
        mea = 0.5 * (af1[indfm] + af2[indfm])
        var = ((af1[indfm] - af2[indfm]) / 2)**2

        # In binomial sampling, the variance on k is var(k) = nx (1 - x), so
        # for the frequency var(k/n) = x (1 - x) / n
        n_all = mea * (1 - mea) / var
        
        # NOTE: pseudocounts that come from the F4 dilution estimate, so we
        # only listen to the data if there is enough data points to listen to
        len_pseudo = 1
        n_pseudo = sample.get_n_templates_dilutions()
        n_allp = np.concatenate([n_all, ([n_pseudo] * len_pseudo)])

        if VERBOSE >= 2:
            print datum['overlap'], 'Number of doubly polymorphic sites:', nsites, 'n_pseudo:', n_pseudo

        # NOTE: the estimate of n has a bad distribution because some points are
        # exactly on the diagonal, so we average the inverse (which is well
        # behaved) and also take the medians as alternatives
        n = 1.0 / (1.0 / n_allp).mean()

        datum['n'] = n


def store_data(data, fn):
    '''Store data to file for the plots'''
    import cPickle as pickle
    with open(fn, 'wb') as f:
        pickle.dump(data, f, protocol=-1)


def load_data(fn):
    '''Load the data for the plots'''
    import cPickle as pickle
    with open(fn, 'rb') as f:
        return pickle.load(f)


def plot_allele_frequency_overlap(data, title='', VERBOSE=0, use_logit=False,
                                  fig_filename=None,
                                  separate_axes=False):
    '''Plot allele frequency in the overlap regions'''
    if VERBOSE >= 2:
        print 'Plot allele frequency in overlaps'
    import matplotlib.pyplot as plt
    import seaborn as sns
    from util import add_panel_label

    sns.set_style('darkgrid')
    colors = sns.color_palette('Set1', 5)
    fs = fig_fontsize
    xmin = 1e-3
    
    if not separate_axes:
        fig, ax = plt.subplots(figsize=(fig_width, 0.8 * fig_width))
        axs = [ax]
    else:
        fig = plt.figure(figsize=(1.5 * fig_width, 1.5 * 0.8 * fig_width))
        fun = fig.add_axes
        lpad = 0.07
        hpad = 0.05
        vpad = 0.05
        width = (1.0 - lpad - 3 * hpad) / 3
        height = (1.0 - 4 * vpad) / 2
        axs = [fun([lpad, 0.5 + vpad, width, height]),
               fun([lpad +  hpad + width, 0.5 + vpad, width, height]),
               fun([lpad + 2 * hpad + 2 * width, 0.5 + vpad, width, height]),
               fun([0.5 * (1 - hpad) - width, vpad, width, height]),
               fun([0.5 * (1 + hpad), vpad, width, height]),
              ]

    # NOTE (Fabio): my logit patch has made it to matplotlib master but now there
    # is a weird bug here that crashes the figure. If I manually call
    # plt.yscale('logit')
    # it works, so I don't quite understand. Anyway, this is only aesthetics and
    # does not affect any data.
    if use_logit:
        for iax, ax in enumerate(axs):
            ax.set_xlim(xmin, 1 - xmin)
            ax.set_ylim(xmin, 1 - xmin)
            ax.set_xscale('logit')
            ax.set_yscale('logit')
            ax.xaxis.set_tick_params(labelsize=fs)
            ax.yaxis.set_tick_params(labelsize=fs)

            if iax not in (0, 3):
                ax.set_yticklabels([])

    for ida, datum in enumerate(data):
        if separate_axes:
            ax = axs[ida]

        afjoint = datum['af']
        color = colors[datum['io']]

        x = afjoint[0].ravel()
        y = afjoint[1].ravel()

        ind = ~(x.mask | y.mask)
        x = x[ind]
        y = y[ind]

        ind = (x >= xmin) & (x <= 1 - xmin) & (y >= xmin) & (y <= 1 - xmin)
        x = x[ind]
        y = y[ind]

        ax.scatter(x, y,
                   s=50,
                   color=color,
                   alpha=0.7,
                   edgecolor='none')

        ## Plot stddev in Poisson sampling
        #n = datum['n']
        #x = np.linspace(np.log10(xmin), 0, 1000)
        #x = 1.0 / (1 + 10**(-x))
        #y = x - np.sqrt(x / n)
        #ax.plot(np.concatenate([x, 1 - y[::-1]]), np.concatenate([y, 1 - x[::-1]]),
        #        lw=3, c=color, alpha=0.5)
        #ax.plot(np.concatenate([y, 1 - x[::-1]]), np.concatenate([x, 1 - y[::-1]]),
        #        lw=3, c=color, alpha=0.5,
        #        label=datum['overlap'])

        if separate_axes or (ida == len(data) - 1):
            ax.plot([xmin, 1 - xmin], [xmin, 1 - xmin], lw=2, color='k', alpha=0.5)
            #ax.set_xlabel('SNP frequency leading fragment', fontsize=fs)
            #ax.set_ylabel('SNP frequency trailing fragment', fontsize=fs)


    if not separate_axes:
        add_panel_label(ax, 'C', x_offset=-0.22)

    if title:
        ax.set_title(title)

    if not separate_axes:
        plt.tight_layout()

    if fig_filename is not None:
        for ext in ['.pdf','.svg', '.png']:
            fig.savefig(fig_filename+ext)
            plt.close(fig)

    else:
        plt.ion()
        plt.show()



# Script
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="make figure for SNP correlations")
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    params = parser.parse_args()

    VERBOSE = 2
    pname = 'p11'
    n_time = 4

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'controls')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'allele_frequency_overlap.pickle'

    if not os.path.isfile(fn_data) or params.redo:
        patient = Patient.load(pname)
        samples = patient.samples[n_time]
        data = get_allele_frequency_overlap(sample, overlaps, cov_min=cov_min,
                                            VERBOSE=VERBOSE, qual_min=qual_min)

        estimate_templates_overlaps(sample, data)


        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    filename = foldername+'allele_frequency_overlap'
    plot_allele_frequency_overlap(data, VERBOSE=VERBOSE,
                                  fig_filename=filename,
                                 )
