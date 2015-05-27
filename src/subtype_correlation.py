import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.HIVreference import HIVreference
from hivevo.hivevo.samples import all_fragments
from util import store_data, load_data, fig_width, fig_fontsize
import os
from filenames import get_figure_folder

def get_quantiles(q, arr):
    from scipy.stats import scoreatpercentile
    thresholds = [scoreatpercentile(arr, 100.0*i/q) for i in range(q+1)]
    return {i: {'range':(thresholds[i],thresholds[i+1]), 
                'ind':((arr>=thresholds[i])*(arr<thresholds[i+1]))}
           for i in range(q)}


def plot_subtype_correlation(data, fig_filename = None, figtypes=['.png', '.svg', '.pdf']):
    ####### plotting ###########
    import seaborn as sns
    from matplotlib import pyplot as plt
    plt.ion()
    sns.set_style('darkgrid')
    figpath = 'figures/'
    fs=fig_fontsize
    fig_size = (fig_width, 0.5*fig_width)
    fig, axs = plt.subplots(1, 2, figsize=fig_size)

    ax=axs[0]
    patients = sorted(data['correlations']['pcode'].unique(), key = lambda x:int(x[1:]))
    colors = {pat:c for pat, c in zip(patients, 
                                      sns.color_palette(n_colors=len(patients)))}

    mean_rho = data['correlations'].groupby(by=['time', 'pcode'], as_index=False).mean().groupby('pcode')
    var_rho = data['correlations'].groupby(by=['time', 'pcode'], as_index=False).var().groupby('pcode')

    for pat in patients:
        ax.errorbar(mean_rho.get_group(pat)['time'],
                    mean_rho.get_group(pat)['rho'],
                    np.sqrt(var_rho.get_group(pat)['rho']),
                    color = colors[pat], ls="none", markersize=8, marker='o', label=pat)

    ax.legend(loc=2, fontsize=fs-3, ncol=2, labelspacing=0.1, columnspacing=0.1)
    ax.set_yticks([0,0.25,0.5, 0.75])
    ax.set_xticks([0,1000,2000,3000])
    ax.set_xlabel('EDI [days]', fontsize=fs)
    ax.set_ylabel("Spearman's r", fontsize=fs)
    for item in ax.get_xticklabels()+ax.get_yticklabels():
        item.set_fontsize(fs)

    ax=axs[1]
    colors = sns.color_palette(n_colors=4)
    time_bins = np.arange(0,4000,500)
    binc = 0.5*(time_bins[1:] + time_bins[:-1])
    div = data['diverse_fraction']
    div.loc[:,'time_bin'] = np.minimum(len(time_bins)-2, np.maximum(0,np.searchsorted(time_bins, div["time"])-1))
    for i in range(4): 
        ent = 'S'+str(i+1)
        div.loc[:,ent] = div.loc[:,ent].astype(float)
        mean_div = div.loc[:,[ent, 'time_bin']].groupby(by=['time_bin'], as_index=False).mean()
        var_div = div.loc[:,[ent, 'time_bin']].groupby(by=['time_bin'], as_index=False).var()
        ax.errorbar(binc, mean_div.loc[:,ent] , np.sqrt(var_div.loc[:,ent]),
                    label='Q'+str(i+1), c=colors[i])

    ax.set_ylim([0,0.3])
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.set_xticks([0,1000,2000,3000])
    ax.set_ylabel('fraction>0.01')
    ax.set_xlabel('EDI [days]', fontsize=fs)
    for item in ax.get_xticklabels()+ax.get_yticklabels():
        item.set_fontsize(fs)
    ax.legend(loc=2, ncol=2,fontsize=fs-3,
              labelspacing=0.1, columnspacing=0.5)

    plt.tight_layout(rect=(0.0, 0.02, 0.98, 0.98), pad=0.05, h_pad=0.5, w_pad=0.4)
    if fig_filename is not None:
        for ext in figtypes:
            fig.savefig(fig_filename+'_sfs'+ext)
    else:
        plt.ion()
        plt.show()

if __name__=="__main__":
    import argparse
    from scipy.stats import spearmanr
    import pandas as pd
    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action = 'store_true', help = 'recalculate data')
    params=parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'subtype_correlation.pickle'

    print(
"""Note: this depends a little on which regions are used. 
stratification is good within regions, but overall entropy 
differences between regions swamp differences within regions.
""")
    if not os.path.isfile(fn_data) or params.redo:
        patients = ['p2', 'p3','p5', 'p8', 'p9', 'p10','p11']
        regions = ['p24', 'p17', 'RT1', 'RT2', 'RT3', 'RT4', 'PR', 
                   'IN1', 'IN2', 'IN3','p15', 'vif', 'nef','gp41','gp1201']
        cov_min = 1000
        af_threshold = 0.01
        hxb2 = HIVreference(refname='HXB2')
        good_pos_in_reference = hxb2.get_ungapped(threshold = 0.05)

        correlations = []
        diverse_fraction = []
        for pi, pcode in enumerate(patients):
            try:
                p = Patient.load(pcode)
            except:
                print "Can't load patient", pcode
            else:
                for region in regions:
                    aft = p.get_allele_frequency_trajectories(region, cov_min=cov_min)
                    if len(aft.mask.shape)<2:
                        aft.mask = np.zeros_like(aft, dtype=bool)

                    # get patient to subtype map and subset entropy vectors
                    patient_to_subtype = p.map_to_external_reference(region, refname = 'HXB2')
                    subtype_entropy = hxb2.get_entropy_in_patient_region(patient_to_subtype)
                    good_ref = good_pos_in_reference[patient_to_subtype[:,0]]
                    # loop over times and calculate the correlation for each value
                    for t, af in izip(p.dsi,aft):
                        patient_entropy = np.maximum(0,-np.sum(af[:5]*np.log(1e-10+af[:5]), axis=0))[patient_to_subtype[:,2]]
                        good_af = (~np.any(af.mask, axis=0)[patient_to_subtype[:,2]]) & good_ref
                        if good_af.sum()>0.5*good_af.shape[0]:
                            rho,pval = spearmanr(patient_entropy[good_af], subtype_entropy[good_af])
                            correlations.append({'pcode':pcode,
                                         'region':region,
                                         'time':t,
                                         'rho':rho,
                                         'pval':pval})
        # determine genome wide fraction of alleles above a threshold
        for pi, pcode in enumerate(patients):
            try:
                p = Patient.load(pcode)
            except:
                print "Can't load patient", pcode
            else:
                for region in regions:
                    aft = p.get_allele_frequency_trajectories(region, cov_min=cov_min)
                    if len(aft.mask.shape)<2:
                        aft.mask = np.zeros_like(aft, dtype=bool)

                    # get patient to subtype map and subset entropy vectors
                    patient_to_subtype = p.map_to_external_reference(region, refname = 'HXB2')
                    subtype_entropy = hxb2.get_entropy_in_patient_region(patient_to_subtype)
                    entropy_quantiles = get_quantiles(4, subtype_entropy)
                    good_ref = good_pos_in_reference[patient_to_subtype[:,0]]
                    # loop over times and calculate the correlation for each value
                    for t, af in izip(p.dsi,aft):
                        good_af = (~np.any(af.mask, axis=0)[patient_to_subtype[:,2]]) & good_ref
                        tmp_af = af[:,patient_to_subtype[:,2]]
                        tmp = {'S'+str(i+1):np.mean(tmp_af[:,Squant['ind']*good_af*good_ref].max(axis=0)\
                                              <tmp_af[:,Squant['ind']*good_af*good_ref].sum(axis=0)-af_threshold)
                                                for i, Squant in entropy_quantiles.iteritems()}
                        tmp.update({'pcode':pcode,'region':region,'time':t})
                        diverse_fraction.append(tmp)

        data={'correlations': pd.DataFrame(correlations),
              'diverse_fraction': pd.DataFrame(diverse_fraction),
              'regions':regions, 'patients':patients, 'cov_min':cov_min, 'threshold':af_threshold}
        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)

plot_subtype_correlation(data, fig_filename=foldername+'entropy_correlation')
