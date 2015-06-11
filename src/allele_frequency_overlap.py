# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/12/14
content:    Plot allele frequencies in overlap of consecutive fragments.
'''
# Modules
import os
import numpy as np

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.controls.check_allele_frequency_overlap import get_allele_frequency_overlap
from util import store_data, load_data, fig_width, fig_fontsize



# Globals
overlaps = ['F1-2', 'F2-3', 'F3-4', 'F4-5', 'F5-6']
cov_min = 1000
qual_min = 30




# Functions
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


def plot_allele_frequency_overlap(data, title='', VERBOSE=0, use_logit=True,
                                  fig_filename=None):
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
    
    fig, ax = plt.subplots(figsize=(fig_width, 0.8 * fig_width))

    legend = set()
    for ida, datum in enumerate(data):
        afjoint = datum['af']
        color = colors[datum['io']]
        if datum['overlap'] not in legend:
            label = datum['overlap']
            #legend.add(datum['overlap'])
        else:
            label = None
        ax.scatter(afjoint[0].ravel(), afjoint[1].ravel(),
                   s=50,
                   color=color,
                   alpha=0.7,
                   #label=label,
                   edgecolor='none')

        # Plot stddev in Poisson sampling
        n = datum['n']
        x = np.linspace(np.log10(xmin), 0, 1000)
        x = 1.0 / (1 + 10**(-x))
        y = x - np.sqrt(x / n)
        ax.plot(np.concatenate([x, 1 - y[::-1]]), np.concatenate([y, 1 - x[::-1]]),
                lw=3, c=color, alpha=0.5)
        ax.plot(np.concatenate([y, 1 - x[::-1]]), np.concatenate([x, 1 - y[::-1]]),
                lw=3, c=color, alpha=0.5,
                label=datum['overlap'])

    
    ax.plot([xmin, 1 - xmin], [xmin, 1 - xmin], lw=2, color='k', alpha=0.5)

    ax.set_xlabel('SNP frequency leading fragment', fontsize=fs)
    ax.set_ylabel('SNP frequency trailing fragment', fontsize=fs)
    ax.grid(True)

    if use_logit:
        ax.set_xscale('logit')
        ax.set_yscale('logit')
        ax.set_xlim(xmin, 1 - xmin)
        ax.set_ylim(xmin, 1 - xmin)

    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.legend(loc=2, fontsize=fs)
    add_panel_label(ax, 'C', x_offset=-0.22)

    if title:
        ax.set_title(title)

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

    VERBOSE = 2
    pname = 'p11'
    n_time = 4

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'controls')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'allele_frequency_overlap.pickle'

    if not os.path.isfile(fn_data):
        patient = load_patient(pname)
        samples = patient.samples.iloc[[n_time]]
        sample = SamplePat(samples.iloc[0])
        data = get_allele_frequency_overlap(samples, overlaps, cov_min=cov_min,
                                            VERBOSE=VERBOSE, qual_min=qual_min)

        estimate_templates_overlaps(sample, data)


        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    filename = foldername+'allele_frequency_overlap'
    plot_allele_frequency_overlap(data, VERBOSE=VERBOSE,
                                  fig_filename=filename,
                                 )
