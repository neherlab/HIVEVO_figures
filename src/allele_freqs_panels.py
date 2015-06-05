# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/03/15
content:    Make three panels with allele frequencies in a short genomic region.
'''
# Modules
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm

from hivevo.hivevo.patients import Patient
from filenames import get_figure_folder
from util import HIVEVO_colormap, store_data, load_data



# Functions
def compress_data(aft, times, pcode, region):
    '''Compress data for plots, discarding useless info'''
    data = []
    datum = {'aft': aft,
             'times': times,
             'pcode': pcode,
             'region': region}
    data.append(datum)

    return data

def plot_allele_freq_example(data, title='', VERBOSE=0, savefig=False):
    '''Plot the frequencies of alleles as trajectories and  
       at 3 representative time points'''
    fig, axs = plt.subplots(2, 3, figsize=(7, 5))
    sns.set_style('darkgrid')
    fs = 18
    datum = data[0]

    # make three panel plot with SNV frequencies at specific time points
    ind = np.arange(len(datum['times']))
    ind = [0, len(ind) // 2, ind[-1]]

    icons0 = datum['aft'][0].argmax(axis=0)

    cmap = HIVEVO_colormap(kind='alternative')
    x = np.arange(datum['aft'].shape[2])
    color = [[float(tmp) for tmp in cmap(p)] for p in np.linspace(0, 1, len(x))]
    for ii, i in enumerate(ind): # loop over times
        ax = axs[0][ii]
        time = datum['times'][i]
        af = datum['aft'][i]

        af_min = [] # plot minor allele frequencies
        for pos, afpos in enumerate(af.T):
            afpos = afpos.copy()
            afpos[icons0[pos]] = 0
            afpos.sort()
            af_min.append(afpos[-1])
        af_min = np.array(af_min)
        ax.scatter(x, af_min, s=100, c=color, edgecolor='none')
        
        ax.set_ylim(1e-2, 1.35)
        ax.set_xlim(-5, len(x) + 5)
        ax.set_xticks(range(0, len(x), 150))
        ax.set_yscale('log')
        ax.grid(True)
        if time>500:
            ax.set_title(str(int(time / 365.25))+' years', fontsize=fs)
        else:
            ax.set_title(str(int(time / 30.5))+' months', fontsize=fs)
        for item in ax.get_xticklabels():
            item.set_fontsize(fs)

        if ii == 0:
            for item in ax.get_yticklabels():
                item.set_fontsize(fs)
        else:
            ax.set_yticklabels([])

    axs[0][1].set_xlabel('Position [bp]', fontsize=fs, labelpad=5)
    fig.text(0.035, 0.5, 'SNP frequency', ha='center', va='center', rotation='vertical',
             fontsize=fs)

    # plot SNV trajectories
    ax = plt.subplot2grid((2, 3), (1, 0), colspan=3)
    tday = datum['times']/365.25
    for pos in xrange(datum['aft'].shape[2]):
        for nuc in xrange(4):
            traj = datum['aft'][:,nuc,pos]
            traj[traj<0.003] = 0.003
            if (traj[0] < 0.5) and (traj.max() > 0.05):
                ax.plot(tday, traj, c=color[pos])

    ax.set_ylim(1e-2, 1.35)
    ax.set_xlim(0, tday[-1] + .1)
    ax.set_xticks([0,2,4,6,8])
    ax.set_yscale('log')
    ax.set_xlabel('ETI [years]', fontsize=fs)
    for item in ax.get_xticklabels() + ax.get_yticklabels():
        item.set_fontsize(fs)

    plt.tight_layout(rect=(0.07, 0.02, 0.98, 0.98), pad=0.05, h_pad=0.5, w_pad=0.4)
    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.ion()
        plt.show()




# Script
if __name__ == '__main__':

    VERBOSE = 2

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'allele_freqs_panels.pickle'

    if not os.path.isfile(fn_data):
        print("Regenerating plot data")
        pcode = 'p3'
        region = 'p17'
        cutoff = 0.01

        patient = Patient.load(pcode)

        # FIXME: the matrix has a compress mask, R. do you know why?
        aft = patient.get_allele_frequency_trajectories(region, error_rate=cutoff)
        times = patient.times()

        data = compress_data(aft, times, pcode, region)
        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)
        
    pcode = data[0]['pcode']
    region = data[0]['region']
    filename = foldername+'_'.join(['allele_freq_example', pcode, region])
    for ext in ['png', 'pdf', 'svg']:
        plot_allele_freq_example(data,
                                    VERBOSE=VERBOSE,
                                    savefig=filename+'.'+ext)

    plot_allele_freq_example(data,
                             VERBOSE=VERBOSE)
