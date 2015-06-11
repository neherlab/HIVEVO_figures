# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Plot allele frequency of control/patient sample.
'''
# Modules
import os
import argparse
import numpy as np

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.samples import load_sample_sequenced as lssp
from hivwholeseq.sequencing.samples import load_sample_sequenced as lss
from hivwholeseq.paper_figures.filenames import get_figure_folder
from filenames import get_figure_folder
from util import store_data, load_data, fig_width, fig_fontsize



# Functions
def compress_data(counts, samplename, fragment, data=[]):
    '''Compress data for plots, discarding useless info'''
    def get_minor_freqs(counts, cov_min=1000):
        '''Get minor freqs from counts'''
        import numpy as np
        af = np.zeros(counts.shape[-1])
        for pos, count in enumerate(counts.T):
            cov = count.sum()
            if cov >= cov_min:
                af_pos = 1.0 * count / cov
                af[pos] = np.sort(af_pos)[-2]
        return af

    datum = {'freq_minor': get_minor_freqs(counts),
             'samplename': samplename,
             'fragment': fragment}
    data.append(datum)

    return data


def plot_minor_allele_example(data, title='', VERBOSE=0, fig_filename=None):
    '''Plot minor allele in a typical sample'''
    import matplotlib.pyplot as plt
    import seaborn as sns
    from util import add_panel_label

    plt.ioff()

    if VERBOSE:
        print 'Plot minor alleles of example sample'

    fig_size = (fig_width, 0.8*fig_width)
    fig, axs = plt.subplots(1, 2,
                            figsize=fig_size,
                            sharey=True,
                            gridspec_kw={'width_ratios': [3, 1]})
    sns.set_style('darkgrid')

    labels = ['control', 'patient']
    alphas = [0.6, 1]
    colors = [sns.color_palette()[i] for i in [2, 0]]
    shapes = ['s', 'o']

    for idat, datum in enumerate(data):
        y = datum['freq_minor']
        x = np.arange(len(y))
        #axs[0].plot(x, y, lw=1.5, alpha=0.8)
        axs[0].scatter(x, y,
                       marker=shapes[idat],
                       lw=1.5, edgecolor='none',
                       facecolor=colors[idat],
                       zorder=idat+1)

        h = np.histogram(y, bins=np.logspace(-4, 0, 27))
        axs[1].barh(h[1][:-1], h[0], (h[1][1:] - h[1][:-1]),
                    color=colors[idat],
                    alpha=alphas[idat],
                    zorder=2 - idat)

    axs[0].set_xlabel('Position [bp]', fontsize=fig_fontsize)
    axs[0].set_ylabel('SNP frequency', fontsize=fig_fontsize)
    axs[0].set_yscale('log')
    axs[0].set_ylim(10**(-4), 1)
    axs[0].set_xlim(-20, y.nonzero()[0][-1] + 21)
    axs[0].grid(True)
    axs[0].tick_params(axis='both', labelsize=fig_fontsize)

    axs[1].set_xlabel('Number of positions', fontsize=fig_fontsize)
    axs[1].grid(True)
    axs[1].set_yscale('log')
    axs[1].set_xlim(0.8, 2 * h[0].max())
    axs[1].set_xscale('log')
    axs[1].tick_params(axis='x', labelsize=fig_fontsize)

    add_panel_label(axs[0], 'C', x_offset=-0.4)
    plt.tight_layout(pad=0.1, h_pad=0.001, w_pad=0.001)

    if title:
        fig.suptitle(title)

    if fig_filename is not None:
        for ext in ['.pdf','.svg', '.png']:
            fig.savefig(fig_filename+"_"+measure+ext)
            plt.close(fig)

    else:
        plt.ion()
        plt.show()



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    params = parser.parse_args()

    fragment = 'F1'
    VERBOSE = 2
    username = os.path.split(os.getenv('HOME'))[-1]

    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'minor_alleles_example.pickle'

    if not os.path.isfile(fn_data) or params.redo:
        samplename = 'NL4-3'
        sample = lss(samplename)
        counts = sample.get_allele_counts(fragment, merge_read_types=True)
        data = compress_data(counts, samplename, fragment)

        samplename = '27134'
        sample = lssp(samplename)
        counts = sample.get_allele_counts(fragment, merge_read_types=True)
        data = compress_data(counts, samplename, fragment, data=data)


        store_data(data, fn_data)
    else:
        data = load_data(fn_data)

    plot_minor_allele_example(data,
                              VERBOSE=VERBOSE)
