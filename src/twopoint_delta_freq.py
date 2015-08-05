# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       05/08/15
content:    Analyze the dynamics of SNP frequencies that are coupled
'''
from __future__ import division, print_function
import numpy as np
import sys,os
from itertools import izip
from hivevo.patients import Patient
from hivevo.samples import all_fragments
from hivevo.af_tools import LD as LDfunc
from filenames import get_figure_folder
from util import store_data, load_data, fig_width, fig_fontsize
import matplotlib.pyplot as plt
import seaborn as sns
plt.ion()
sns.set_style('darkgrid')

# Script
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    params = parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'LD.pickle'
    patients = ['p' +str(i) for i in xrange(1,12) if i not in [4,7]]

    dmin = 20
    dmin_pad = 100
    var_min = 0.1
    cov_min = 200
    corr_vs_distance = {}
    bins = np.arange(0,401,40)
    binc = (bins[:-1]+bins[1:])*0.5
    all_dists = []
    all_weights = []
    for frag in all_fragments:
        if frag not in ['F'+str(i) for i in xrange(1,7)]:
            continue
        dists = []
        weights = []
        for pcode in patients:
            p = Patient.load(pcode)
            aft = p.get_allele_frequency_trajectories(frag)
            depth = p.get_fragment_depth(pad=False, limit_to_dilution=False)
            depth_pad = p.get_fragment_depth(pad=True, limit_to_dilution=False)
            for si, sample in enumerate(p.samples[:-1]):
                if depth[si][all_fragments.index(frag)]>dmin \
                    or depth_pad[si][all_fragments.index(frag)]>dmin_pad:
                    try:
                        positions, af2p, cov, af1p = sample.get_pair_frequencies(frag, var_min=var_min)
                        majority_nuc = af1p.argmax(axis=0)
                        if positions is None:
                            continue
                        LD, Dp, p12 =  LDfunc(af2p, af1p, cov, cov_min=100)
                        daf = aft[si+1][majority_nuc,positions] - aft[si][majority_nuc,positions]
                        X,Y = np.meshgrid(positions, positions)
                        dp1,dp2 = np.meshgrid(daf, daf)
                        dists.extend(np.abs(X-Y)[cov>cov_min])
                        weights.extend((np.sign(LD)*np.sign(dp1*dp2))[cov>cov_min])
                        print (pcode, si, frag,
                               " # of positions:", len(positions),
                               'depth:', depth[si][all_fragments.index(frag)])
                    except:
                        print('no variable sites')
                else:
                    print (pcode, si, frag, "insufficient depth:",
                           depth[si][all_fragments.index(frag)],
                           depth_pad[si][all_fragments.index(frag)])
        # prune bad sites
        weightsa = np.array(weights)
        distsa = np.array(dists)[~np.isnan(weightsa)]
        weightsa = weightsa[~np.isnan(weightsa)]
        all_weights.extend(weightsa)
        all_dists.extend(distsa)
        yn,xn = np.histogram(distsa, bins=bins)
        y,x = np.histogram(distsa, weights = weightsa, bins=bins)
        corr_vs_distance[frag]=y/(1e-10+yn)

    yn,xn = np.histogram(all_dists, bins=bins)
    y,x = np.histogram(all_dists, weights = all_weights, bins=bins)
    corr_vs_distance['all']=y/(1e-10+yn)

    plt.figure()
    for frag in all_fragments:
        plt.plot(binc, corr_vs_distance[frag], label=frag, lw=2)
    plt.legend()

    plt.figure()
    plt.plot(binc, corr_vs_distance['all'], lw=2)
    plt.ylabel('correlation of SNP frequency change')
    plt.xlabel('distance [bp]')
