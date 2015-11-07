# Modules
import os, sys
import numpy as np
import pandas as pd
from itertools import izip
from scipy.stats import spearmanr

from hivevo.hivevo.patients import Patient
from hivevo.hivevo.HIVreference import HIVreference, HIVreferenceAminoacid
from hivevo.hivevo.samples import all_fragments

from util import store_data, load_data, fig_width, fig_fontsize, get_quantiles, add_panel_label, patient_colors, patients
from filenames import get_figure_folder



# Functions
def collect_correlations(patients, regions, cov_min=1000, refname='HXB2', min_dsi=1500):
    '''Correlation of entropy between patients'''
    ps = [Patient.load(pcode) for pcode in patients]

    correlations = []
    for region in regions:
        print region
        for pi, p1 in enumerate(ps):
            aft1 = p1.get_allele_frequency_trajectories(region, cov_min=cov_min)
            af1 = aft1[p1.dsi >= min_dsi].mean(axis=0)
            en1 = np.maximum(0,-np.sum(af1[:-1]*np.log(1e-10+af1[:-1]), axis=0))
            ptoref1 = p1.map_to_external_reference(region, refname=refname)
            ptorefd1 = dict(ptoref1[:, ::2])
            seq1 = p1.get_initial_sequence(region)

            for p2 in ps[:pi]:
                aft2 = p2.get_allele_frequency_trajectories(region, cov_min=cov_min)
                af2 = aft2[p2.dsi >= min_dsi].mean(axis=0)
                en2 = np.maximum(0,-np.sum(af2[:-1]*np.log(1e-10+af2[:-1]), axis=0))
                ptoref2 = p2.map_to_external_reference(region, refname=refname)
                ptorefd2 = dict(ptoref2[:, ::2])
                seq2 = p2.get_initial_sequence(region)

                overlap = np.intersect1d(ptoref1[:, 0], ptoref2[:, 0], assume_unique=True)
                af_ov = np.array([(en1[ptorefd1[pos]], en2[ptorefd2[pos]]) for pos in overlap])
                rho, pval = spearmanr(af_ov[:, 0], af_ov[:, 1])

                seq1_ov = np.array([seq1[ptorefd1[pos]] for pos in overlap])
                seq2_ov = np.array([seq2[ptorefd2[pos]] for pos in overlap])
                dist = (seq1_ov != seq2_ov).mean()

                correlations.append({'pcode1': p1.name,
                                     'pcode2': p2.name,
                                     'pcode': p1.name+'-'+p2.name,
                                     'region': region,
                                     'rho': rho,
                                     'distance': dist,
                                     'pval': pval})

    return pd.DataFrame(correlations)


def plot_correlation(data, fig_filename=None, figtypes=['.png', '.svg', '.pdf']):
    '''Plot results'''
    import seaborn as sns
    from matplotlib import pyplot as plt
    plt.ion()

    sns.set_style('darkgrid')
    figpath = 'figures/'
    fs=fig_fontsize
    fig_size = (fig_width, 0.5*fig_width)
    fig, ax = plt.subplots(figsize=fig_size)

    patients = sorted(data['correlations']['pcode'].unique(),
                      key=lambda x:1000 * int(x.split('-')[0][1:]) + int(x.split('-')[1][1:]))
    colorB, colorC, colorAE, colorC_AE = sns.color_palette('deep', n_colors=4) 

    # calculate mean and variance across regions for each time point and patient
    mean_rho = data['correlations'].groupby(by=['pcode'], as_index=True).mean()
    var_rho = data['correlations'].groupby(by=['pcode'], as_index=True).var()

    # loop over patients and plot the mean/std of the previously grouped data 
    for pat in patients:
        pats = pat.split('-')
        if ('p1' in pats) and ('p6' in pats):
            color = colorC_AE
            z = 10
        elif 'p1' in pats:
            color = colorAE
            z = 5
        elif 'p6' in pats:
            color = colorC
            z = 6
        else:
            color = colorB
            z = 4
        ax.errorbar(np.array(mean_rho.loc[pat, 'distance']),
                    np.array(mean_rho.loc[pat, 'rho']),
                    yerr=np.array(np.sqrt(var_rho.loc[pat, 'rho'])),
                    color=color, alpha=0.6,
                    ls="none",
                    zorder=z,
                    markersize=8, marker='o')

    hs = [plt.Line2D((0,0),(0,0), color=colorB, marker='o', ls=''),
          plt.Line2D((0,0),(0,0), color=colorC, marker='o', ls=''),
          plt.Line2D((0,0),(0,0), color=colorAE, marker='o', ls=''),
          plt.Line2D((0,0),(0,0), color=colorC_AE, marker='o', ls='')]
    ls = ['Only subtype B', 'Subtypes B and C', 'Subtypes B and 01_AE', 'Subtypes C and 01_AE']
    ax.legend(hs, ls,
              loc=1, fontsize=fs-3)
    ax.set_yticks([0,0.25,0.5])
    ax.set_xticks([0, 0.1, 0.2])
    ax.set_xlim([0, 0.3])
    ax.set_ylim([-0.1, 0.8])
    ax.set_xlabel('Hamming distance between founders', fontsize=fs)
    ax.set_title(r"Spearman's $\rho$", fontsize=fs)
    for item in ax.get_xticklabels()+ax.get_yticklabels():
        item.set_fontsize(fs)

    # plot output
    plt.tight_layout(rect=(0.0, 0.02, 0.98, 0.98), pad=0.05, h_pad=0.5, w_pad=0.4)
    if fig_filename is not None:
        for ext in figtypes:
            fig.savefig(fig_filename+ext)
    else:
        plt.ion()
        plt.show()



# Script
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="make figure for SNP correlations")
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    # TODO: add choice between nucleotides and amino acids
    parser.add_argument('--type', choices=['nuc', 'aa'], default='nuc',
                        help='Sequence type (nuc or aa)')
    parser.add_argument('--reference', choices=['HXB2', 'NL4-3'], default='HXB2',
                        help='Reference')

    params = parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'interpatient_correlation'
    if params.reference != 'HXB2':
        fn_data = fn_data + '_'+params.reference
    if params.type == 'aa':
        fn_data = fn_data + '_aa'
    fn_data = fn_data + '.pickle'

    if not os.path.isfile(fn_data) or params.redo:
        regions = ['p17', 'p24', 'PR', 'RT', 'p15', 'IN', 'vif', 'gp41', 'gp120', 'nef']
        #regions = ['p24', 'p17', 'RT1', 'RT2', 'RT3', 'RT4', 'PR', 
        #           'IN1', 'IN2', 'IN3','p15', 'vif', 'nef','gp41','gp1201']
        cov_min = 1000

        if params.type == 'nuc':
            # determine correlations between intra patient diversity and subtype diversity
            correlations = collect_correlations(patients, regions, cov_min=cov_min,
                                                refname=params.reference)

        else:
            pass
            #correlations = collect_correlations_aminoacids(patients, regions, cov_min=cov_min,
            #                                    refname=params.reference,
            #                                              )


        data={'correlations': correlations,
              'regions':regions,
              'patients': patients,
              'cov_min': cov_min}
        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)

fig_filename = foldername+'entropy_correlation_interpatient'
if params.type == 'aa':
    fig_filename = fig_filename + '_aa'
plot_correlation(data, fig_filename=fig_filename)
