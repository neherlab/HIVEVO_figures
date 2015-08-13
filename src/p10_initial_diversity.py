'''
p10 has suspicioulsy high diversity in the first time point, which at the
same time is thought to be really early in infection. This script analyzes
this diversity to establish the cause
'''

import numpy as np
from itertools import izip
from hivevo.patients import Patient
from hivevo.samples import all_fragments
from hivevo.af_tools import LD as LDfunc
from util import store_data, load_data, fig_width, fig_fontsize, HIVEVO_colormap
import os
from filenames import get_figure_folder
import matplotlib.pyplot as plt
import seaborn as sns
plt.ion()

sns.set_style('darkgrid')
cols = HIVEVO_colormap()

if __name__=="__main__":
    p = Patient.load('p10')
    aft = p.get_allele_frequency_trajectories('genomewide')

    af = aft[0]
    consensus_indices = p.get_initial_indices('genomewide')
    minor_af = 1.0 - af.max(axis=0)

    # make a histogram of the minor allele frequencies
    plt.figure()
    plt.hist(minor_af, bins = np.linspace(0,1,51), bottom=0.5)
    plt.yscale('log')
    plt.xlabel('frequency')
    plt.ylabel('number of minor variants')

    # --> there are two clear peaks one around 0.35-0.5 , the other around 0.1-0.15
    variable_pos = minor_af>0.05
    print("number of variable positions:",variable_pos.sum())
    peak1 = minor_af>0.3
    peak1_ii = ((af>0.3)&(af<0.5)).argmax(axis=0)
    peak2 = (minor_af>0.05)&(minor_af<0.2)
    peak2_ii = ((af>0.05)&(af<0.2)).argmax(axis=0)
    print("there are two clear peaks at frequency about 0.15 and 0.4")
    print("peak 1:",peak1.sum())
    print("peak 2:",peak2.sum())
    print("Mutations from both peaks show up in later samples")
    print(" --> hence they are unlikely contaminants")
    # --> these amount to about 100 positions, corresponding to diversity of about 0.5 to 1%
    # --> this seems consistent with donor diversity

    # trajectories of these mutations
    plt.figure()
    plt.title('frequencies trajectories of high peak')
    for pos in np.where(peak1)[0]:
        plt.plot(p.ysi, aft[:,peak1_ii[pos], pos], c=cols(pos*0.0001))
    plt.xlabel('ETI[years]')

    plt.figure()
    plt.title('frequencies trajectories of low peak')
    for pos in np.where(peak2)[0]:
        plt.plot(p.ysi, aft[:,peak2_ii[pos], pos], c=cols(pos*0.0001))
    plt.xlabel('ETI[years]')


    # look at LD between these mutations
    print("LD in the first sample is next to complete. this is consistent with")
    print("  * very early infection with a small number of variants")
    print("  * not yet enough to time to recombine")
    print("  * it is effectively a recombination control with a Patient sample")
    first_sample = p.samples[0]
    fig, axs = plt.subplots(2,3)
    fig.suptitle('Linkage between mutations is strong')
    dists = []
    weights_LD = []
    weights_Dp = []
    cov_min=100
    bins = np.arange(0,401,40)
    binc = (bins[:-1]+bins[1:])*0.5
    for fi, frag in enumerate(all_fragments):
        ax = axs[fi//3][fi%3]
        positions, af2p, cov, af1p = first_sample.get_pair_frequencies(frag, var_min=0.4)
        if positions is not None:
            LD, Dp, p12 =  LDfunc(af2p, af1p, cov, cov_min=100)
            ax.imshow(LD, interpolation='nearest', cmap='jet', vmin=0, vmax=1.0)

        X,Y = np.meshgrid(positions, positions)
        np.fill_diagonal(cov, 0)
        dists.extend(np.abs(X-Y)[cov>=cov_min])
        weights_LD.extend(LD[cov>=cov_min])
        weights_Dp.extend(Dp[cov>=cov_min])

    yn,xn = np.histogram(dists, bins = bins)
    y,x = np.histogram(dists, weights = weights_LD, bins=bins)
    LD_vs_distance = y/(1e-10+yn)
    y,x = np.histogram(dists, weights = weights_Dp, bins=bins)
    Dp_vs_distance=y/(1e-10+yn)

    plt.figure()
    plt.plot(binc, LD_vs_distance, label='r^2')
    plt.plot(binc, Dp_vs_distance, label='D')
    plt.xlabel('distance [bp]')
    plt.ylim([0,1])
    # --> these positions are on strong long range LD, effectively another recombination control. 

