import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.HIVreference import HIVreference
from hivevo.hivevo.af_tools import divergence
from util import store_data, load_data, fig_width, fig_fontsize, \
                 add_panel_label ,add_binned_column, HIVEVO_colormap
from util import boot_strap_patients, replicate_func
import os
from filenames import get_figure_folder
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
plt.ion()
sns.set_style('darkgrid')


def get_toaway_histograms(subtype, Sc=1):
    '''
    calculate allele frequency histograms for each patient and each time points
    separately for sites that agree or disagree with consensus.
    this can be done for a low and high entropy category with the threshold set by Sc
    '''
    away_histogram = {(pcode, Sbin):{} for Sbin in ['low','high'] for pcode in patients}
    to_histogram = {(pcode, Sbin):{} for Sbin in ['low','high'] for pcode in patients}
    # if subtypes == 'any' meaning comparison to groupM, we can load the reference here
    if subtype=='any':
        hxb2 = HIVreference(refname='HXB2', subtype = subtype)
        good_pos_in_reference = hxb2.get_ungapped(threshold = 0.05)

    # determine divergence and minor variation at sites that agree with consensus or not
    for pi, pcode in enumerate(patients):
        try:
            p = Patient.load(pcode)
        except:
            print "Can't load patient", pcode
        else:
            print('subtype:', subtype, "patient",pcode)
            if subtype == 'patient': # if we take the subtype of the patient, load specific ref alignment here
                hxb2 = HIVreference(refname='HXB2', subtype = p['Subtype'])
                good_pos_in_reference = hxb2.get_ungapped(threshold = 0.05)
            for region in regions:
                aft = p.get_allele_frequency_trajectories(region, cov_min=cov_min)

                # get patient to subtype map and subset entropy vectors, convert to bits
                patient_to_subtype = p.map_to_external_reference(region, refname = 'HXB2')
                subtype_entropy = hxb2.get_entropy_in_patient_region(patient_to_subtype)/np.log(2.0)
                ancestral = p.get_initial_indices(region)[patient_to_subtype[:,2]]
                consensus = hxb2.get_consensus_indices_in_patient_region(patient_to_subtype)
                good_ref = good_pos_in_reference[patient_to_subtype[:,0]]
                away_sites = ancestral==consensus
                aft_HXB2 = aft[:,:,patient_to_subtype[:,2]]

                for H, sites in [(away_histogram, away_sites), (to_histogram, ~away_sites)]:
                    for Sbin in ['low', 'high']:
                        if Sbin=='low':
                            ind = (sites)&(subtype_entropy<Sc)&(good_ref)
                        else:                    
                            ind = (sites)&(subtype_entropy>=Sc)&(good_ref)
                        for ti,t in enumerate(p.dsi):
                            y,x = np.histogram(aft_HXB2[ti,ancestral[ind],np.where(ind)[0]].compressed(), bins=af_bins)
                            H[(pcode, Sbin)][t]=y

    return to_histogram, away_histogram

def bin_time(freq_arrays, time_bins):
    binned_hists = [np.zeros_like(af_binc) for ti in time_bins[1:]]
    for hists in freq_arrays.values():
        for t, y in hists.iteritems():
            ti = np.searchsorted(time_bins, t)
            if ti>0 and ti<len(time_bins):
                binned_hists[ti-1]+=y

    return binned_hists

def plot_spectra():
    def labelfunc(ti, tbins):
        if ti==0: return '<'+str(tbins[ti+1])+'y'
        elif ti==len(tbins)-2: return '>'+str(tbins[-2])+'y'
        else: return str(tbins[ti])+'-'+str(tbins[ti+1])+'y'

    for Sbin in ['low']: #, 'high']:
        fig = plt.figure(figsize = (fig_width, 0.8*fig_width))
        ax = plt.subplot(1,1,1)
        cols = HIVEVO_colormap()
        colors = [cols(x) for x in [0.0, 0.5, 0.99]]
        for toaway, ls, sym, H in [(r'founder $=$ group M', '-','o', away_histogram), (r'founder $\neq$ group M ', '--','d', to_histogram)]:
            for ti,tbin in enumerate(tbins[1:]):
                freqs = H[(tbin,Sbin)]
                counts = freqs.sum()
                freqs = 1.0*freqs/counts
                ax.plot(binc, freqs, '-'+sym, label = toaway+' ETI: '+labelfunc(ti,tbins), c=colors[ti], 
                        ls = ls )

        ax.set_yscale('log')
        ax.set_ylabel('fraction of sites', fontsize = fig_fontsize)
        ax.set_xlabel('frequency of founder allele', fontsize = fig_fontsize)
        ax.legend(loc=9, ncol=1,numpoints=2, labelspacing=0.0, fontsize=fig_fontsize)
        plt.tight_layout()
        plt.savefig('figures/reversions.pdf')

def plot_divergence(time_bins, to_histogram,away_histogram):
    def get_div(afhist, fixed=False):
        '''return the fraction of fixed alleles or the mean divergence'''
        if fixed:
            return afhist[0]/afhist.sum()
        else:
            return np.array(afhist[:-1]*(1-af_binc[:-1])).sum()/afhist.sum()


    from random import choice
    fig = plt.figure(figsize = (fig_width, 0.8*fig_width))
    ax = plt.subplot(1,1,1)
    time_binc = 0.5*(time_bins[1:]+time_bins[:-1])
    sym='o'
    ls='-'
    fs = fig_fontsize
    nreps=100
    for subtype in away_histogram:
        for toaway, H in [(r'founder $=$ '+('group M' if subtype=='any' else 'subtype'),  away_histogram[subtype]), 
                          (r'founder $\neq$ '+('group M' if subtype=='any' else 'subtype'), to_histogram[subtype])]:
            mean_hists = bin_time(H,time_bins)
            div = [get_div(mean_hists[ti]) for ti in range(len(time_bins)-1)]
            # make replicates and calculate bootstrap confidence intervals
            replicates = []
            all_keys = H.keys()
            for ri in xrange(nreps):
                bootstrap_keys = [all_keys[ii] for ii in np.random.randint(len(all_keys), size=len(all_keys))]
                tmp = bin_time({key:H[key] for key in bootstrap_keys}, time_bins)
                replicates.append([get_div(tmp[ti]) for ti in range(len(time_bins)-1)])
            std_dev = np.array(replicates).std(axis=0)
            ax.errorbar(time_binc/365.25, div, std_dev, label = toaway, ls = ls, lw=2)

    plt.xlabel('ETI [years]', fontsize=fs)
    plt.ylabel('divergence from founder sequence', fontsize=fs)
    plt.legend(loc=2, fontsize=fs)
    ax.tick_params(labelsize=fs)
    ax.tick_params(labelsize=fs)
    plt.tight_layout()
    for ext in ['.png', '.svg', '.pdf']:
        fig.savefig('figures/to_away_vs_time'+ext)

if __name__=='__main__':

    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action = 'store_true', help = 'recalculate data')
    params=parser.parse_args()

    Sc=10.5
    cov_min=200
    #patients = ['p1', 'p6'] # other subtypes
    patients = ['p1',  'p2', 'p5', 'p6', 'p8', 'p9','p11'] # all patients with initially homogeneous samples
    regions = ['genomewide'] #'pol', 'gag', 'nef']
    af_bins = np.linspace(0,1,11)
    af_binc = 0.5*(af_bins[:-1]+af_bins[1:])
    time_bins = np.array([-10,500,1000,1500,2000,2500])


    away_histogram = {}; to_histogram={}
    for subtype in ['any', 'patient']:
        to_histogram[subtype], away_histogram[subtype] = get_toaway_histograms(subtype, Sc=10)

    plot_divergence(time_bins, to_histogram, away_histogram)
