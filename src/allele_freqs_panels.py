# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/03/15
content:    Make three panels with allele frequencies in a short genomic region.
'''
# Modules
import os, argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm

from hivevo.patients import Patient
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


def plot_allele_freq_example(data, title='', VERBOSE=0, fig_filename=None,
                             figtypes=['.png', '.svg', '.pdf']):
    '''Plot the frequencies of alleles as trajectories and  
       at 3 representative time points'''
    fig, axs = plt.subplots(2, 3, figsize=(7, 5))
    sns.set_style('darkgrid')
    fs = 18
    cmap = HIVEVO_colormap(kind='alternative')

    datum = data[0]

    # make three panel plot with SNV frequencies at specific time points
    ind = np.arange(len(datum['times']))
    ind = [0, len(ind) // 2, ind[-1]]

    icons0 = datum['aft'][0].argmax(axis=0)

    # SNP frequencies in panels
    x = np.arange(datum['aft'].shape[2])
    color = [[float(tmp) for tmp in cmap(p)] for p in np.linspace(0, 1, len(x))]
    for ii, i in enumerate(ind): # loop over times
        ax = axs[0][ii]
        time = datum['times'][i]
        af = datum['aft'][i]

        af_min = []
        for pos, afpos in enumerate(af.T):
            afpos = afpos.copy()
            afpos[icons0[pos]] = 0
            afpos.sort()
            af_min.append(afpos[-1])
        af_min = np.array(af_min)
        ax.scatter(x, af_min, s=100, c=color, edgecolor='none')
        
        ax.set_ylim(1e-2, 1.35)
        ax.set_xlim(-5, len(x) + 20)
        ax.set_xticks(range(0, len(x), 150))
        ax.set_yscale('log')
        ax.grid(True)
        if time>500: # label first time point as month, later ones as years
            ax.set_title(str(int(time / 365.25))+' years', fontsize=fs)
        else:
            ax.set_title(str(int(time / 30.5))+' months', fontsize=fs)
        ax.xaxis.set_tick_params(labelsize=fs)

        if ii == 0:
            ax.yaxis.set_tick_params(labelsize=fs)
        else:
            ax.set_yticklabels([])

    axs[0][1].set_xlabel('Position [bp]', fontsize=fs, labelpad=5)

    # plot SNP trajectories
    fig.text(0.035, 0.5, 'SNP frequency', ha='center', va='center',
             rotation='vertical',
             fontsize=fs)

    # SNP trajectories
    ax = plt.subplot2grid((2, 3), (1, 0), colspan=3)
    tyears = datum['times']/365.25
    for pos in xrange(datum['aft'].shape[2]):
        for nuc in xrange(4):
            traj = datum['aft'][:,nuc,pos]
            traj[traj<0.007] = 0.007
            if (traj[0] < 0.5) and (traj.max() > 0.05):
                ax.plot(tyears, traj, c=color[pos], 
                        alpha = max(0.2,1+np.log10(traj.max())/1.0),
                        lw = max(1,3+1*np.log10(traj.max())))

    ax.set_ylim(1e-2, 1.35)
    ax.set_xlim(0, tyears[-1] + .1)
    ax.set_xticks([0,2,4,6,8])
    ax.set_yscale('log')

    ax.set_xlabel('ETI [years]', fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)

    # Final touches
    plt.tight_layout(rect=(0.07, 0.02, 0.98, 0.98), pad=0.05, h_pad=0.5, w_pad=0.4)
    if fig_filename is not None:
        for ext in figtypes:
            fig.savefig(fig_filename+ext)
        #plt.close(fig)
    else:
        plt.ion()
        plt.show()




# Script
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    params = parser.parse_args()

    VERBOSE = 2

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'allele_freqs_panels.pickle'

    if not os.path.isfile(fn_data) or params.redo:
        print("Regenerating plot data")
        pcode = 'p1'
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
    plot_allele_freq_example(data,
                             VERBOSE=VERBOSE,
                             fig_filename=filename)
