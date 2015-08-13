import numpy as np
from itertools import izip
from hivevo.patients import Patient
from hivevo.HIVreference import HIVreference
from hivevo.af_tools import divergence
from util import store_data, load_data, fig_width, fig_fontsize, add_panel_label ,add_binned_column,HIVEVO_colormap
from util import boot_strap_patients, replicate_func
import os
from filenames import get_figure_folder

def collect_to_away(patients, regions, Sbins=[0,0.02, 0.08, 0.25, 2], cov_min=1000, subtype = 'patient'):
    '''
    collect minor variant frequencies, divergences, etc separately for sites that agree or disagree
    with consensus. consensus is either group M consensus (subtype='any') or the subtype of the 
    respective patient (subtype='patient'). In addition, these quantities are stratified by entropy
    '''
    minor_variants = []
    to_away_divergence = []
    to_away_minor = []
    consensus_distance = {}
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
                away_sites = ancestral==consensus
                good_ref = good_pos_in_reference[patient_to_subtype[:,0]]
                consensus_distance[(pcode, region)] = np.mean(~away_sites)
                print pcode, region, "dist:",1-away_sites.mean(), "useful_ref:",good_ref.mean()

                # loop over times and calculate the af in entropy bins
                for t, af in izip(p.dsi,aft):
                    good_af = (((~np.any(af.mask, axis=0))
                                #&(aft[0].max(axis=0)>0.9)
                                &(af.argmax(axis=0)<4))[patient_to_subtype[:,2]]) \
                                & good_ref
                    # make version of all arrays that contain only unmasked sites and are also ungapped
                    clean_af = af[:,patient_to_subtype[:,2]][:5,good_af]
                    clean_away = away_sites[good_af]
                    clean_consensus = consensus[good_af]
                    clean_ancestral = ancestral[good_af]
                    clean_entropy = subtype_entropy[good_af]
                    clean_entropy_bins = [(clean_entropy>=t_lower)&(clean_entropy<t_upper)
                                        for t_lower, t_upper in zip(Sbins[:-1], Sbins[1:])]
                    clean_minor = clean_af.sum(axis=0) - clean_af.max(axis=0)
                    clean_derived = clean_af.sum(axis=0) - clean_af[clean_ancestral,np.arange(clean_ancestral.shape[0])]
                    print pcode, region, t
                    
                    # for each entropy bin, calculate the average divergence and minor variation
                    for sbin, sites in enumerate(clean_entropy_bins):
                        minor_variants.append({'pcode':pcode,'region':region,'time':t,'S_bin':sbin,
                                            'af_away_minor':  np.mean(clean_minor[sites&clean_away]), 
                                            'af_away_derived':np.mean(clean_derived[sites&clean_away]),
                                            'af_to_minor':    np.mean(clean_minor[sites&(~clean_away)]), 
                                            'af_to_derived':  np.mean(clean_derived[sites&(~clean_away)])})

                    # calculate the minor variation at sites were the founder differs from consensus
                    # in different allele frequency bins
                    clean_reversion = clean_af[clean_consensus,np.arange(clean_consensus.shape[0])]*(~clean_away)
                    clean_total_divergence = clean_af.sum(axis=0) - clean_af[clean_ancestral,np.arange(clean_ancestral.shape[0])]
                    to_away_divergence.append({'pcode':pcode,'region':region,'time':t,
                                        'reversion':np.mean(clean_reversion), 
                                        'divergence':np.mean(clean_total_divergence)})

                    af_thres = [0,0.05,0.1, 0.25, 0.5, 0.95, 1.0]
                    rev_tmp = clean_af[clean_consensus,np.arange(clean_consensus.shape[0])][~clean_away]
                    der_tmp = clean_derived[~clean_away] 
                    for ai,(af_lower, af_upper) in enumerate(zip(af_thres[:-1], af_thres[1:])):
                        to_away_minor.append({'pcode':pcode,'region':region,'time':t,'af_bin':ai,
                                        'reversion_spectrum':np.mean(rev_tmp*(rev_tmp>=af_lower)*(rev_tmp<af_upper)),
                                        'minor_reversion_spectrum':np.mean(der_tmp*(der_tmp>=af_lower)*(der_tmp<af_upper))})

    return pd.DataFrame(minor_variants), pd.DataFrame(to_away_divergence),pd.DataFrame(to_away_minor), consensus_distance



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


def plot_to_away(data, fig_filename = None, figtypes=['.png', '.svg', '.pdf']):
    ####### plotting ###########
    import seaborn as sns
    from matplotlib import pyplot as plt
    plt.ion()
    sns.set_style('darkgrid')
    figpath = 'figures/'
    fs=fig_fontsize
    fig_size = (1.4*fig_width, 0.8*fig_width)
    fig, axs = plt.subplots(1, 2, figsize=fig_size)
    nbs=100

    cols = HIVEVO_colormap()
    colors = [cols(x) for x in [0.0, 0.33, 0.66, 0.99]]

    ####################################################################################
    ##### make panel divergence vs entropy
    ####################################################################################
    ax=axs[1]
    Sbins = np.array([0,0.02, 0.08, 0.25, 2])
    Sbinc = 0.5*(Sbins[1:]+Sbins[:-1])
    def get_Sbin_mean(df):
        return df.groupby(by=['S_bin'], as_index=False).mean()
    color_count = 0
    for lblstr, subtype, ls in [('subtype', 'patient', '--'), ('group M', 'any', '-')]:
        mv = data[subtype]['minor_variants']
        # subset to a specific time interval
        mv = mv.loc[(mv.loc[:,'time']>1500)&(mv.loc[:,'time']<2500),:]
        print "average time:", mv.loc[:,'time'].mean()/365
        mv.loc[:,['af_away_minor', 'af_away_derived', 'af_to_minor', 'af_to_derived']] = \
            mv.loc[:,['af_away_minor', 'af_away_derived', 'af_to_minor', 'af_to_derived']].astype(float)
        mean_to_away =get_Sbin_mean(mv)
        bs = boot_strap_patients(mv, eval_func=get_Sbin_mean, n_bootstrap = nbs, 
                             columns=['af_away_minor', 'af_away_derived', 'af_to_minor', 'af_to_derived', 'S_bin'])

        print mean_to_away
        col = 'af_away_derived'
        ax.errorbar(Sbinc, mean_to_away.loc[:,col], 
                    replicate_func(bs, col, np.std, bin_index='S_bin'), ls=ls, 
                    lw = 3, label = 'founder = '+lblstr, c=colors[color_count])
        color_count+=1
        col = 'af_to_derived'
        ax.errorbar(Sbinc, mean_to_away.loc[:,col], 
                    replicate_func(bs, col, np.std, bin_index='S_bin'), ls=ls,
                    lw = 3, label = u'founder \u2260 '+lblstr, c=colors[color_count])
        color_count+=1
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Divergence from founder sequence', fontsize = fig_fontsize)
    ax.set_xlabel('Variability [bits]', fontsize = fig_fontsize)
    add_panel_label(ax, 'B', x_offset=-0.27)
    for item in ax.get_yticklabels()+ax.get_xticklabels():
        item.set_fontsize(fs)
    ax.set_xlim([0.005,2])

    ####################################################################################
    ##### print reversion statistics
    ####################################################################################

    def get_time_bin_means(df):
        return df.loc[:,['divergence', 'reversion','time_bin']].groupby(by=['time_bin'], as_index=False).mean()
    for subtype in ['patient', 'any']:
        to_away = data[subtype]['to_away']
        time_bins = np.array([0,500,1000,1500, 2500, 3500])
        binc = 0.5*(time_bins[1:]+time_bins[:-1])
        add_binned_column(to_away, time_bins, 'time')
        to_away.loc[:,['reversion', 'divergence']] = \
                to_away.loc[:,['reversion', 'divergence']].astype(float)
        rev_div = get_time_bin_means(to_away)
        bs = boot_strap_patients(to_away, get_time_bin_means,  n_bootstrap = nbs, 
                                 columns = ['reversion','divergence','time_bin'])
        reversion_std = replicate_func(bs, 'reversion', np.std, bin_index='time_bin')
        total_div_std = replicate_func(bs, 'divergence', np.std, bin_index='time_bin')
        fraction = rev_div.loc[:,'reversion']/rev_div.loc[:,'divergence']
        print "Comparison:", subtype
        print "Reversions:\n", rev_div.loc[:,'reversion']
        print "Divergence:\n", rev_div.loc[:,'divergence']
        print "Fraction:"
        for frac, total, num_std, denom_std in zip(fraction, rev_div.loc[:,'divergence'],reversion_std, total_div_std):
            print frac, '+/-', np.sqrt(num_std**2/total**2 + denom_std**2*frac**2/total**2)
        #print reversion_std,total_div_std
        print "Consensus!=Founder:",np.mean(data[subtype]['consensus_distance'].values())

    ####################################################################################
    ##### make panel divergence vs time
    ####################################################################################

    to_histogram=data['to_histogram']
    away_histogram=data['away_histogram']
    time_bins=data['time_bins']
    af_bins=data['af_bins']
    af_binc=0.5*(af_bins[1:]+af_bins[:-1])

    def bin_time(freq_arrays, time_bins):
        binned_hists = [np.zeros_like(af_binc) for ti in time_bins[1:]]
        for hists in freq_arrays.values():
            for t, y in hists.iteritems():
                ti = np.searchsorted(time_bins, t)
                if ti>0 and ti<len(time_bins):
                    binned_hists[ti-1]+=y

        return binned_hists

    def get_div(afhist, fixed=False):
        '''return the fraction of fixed alleles or the mean divergence'''
        if fixed:
            return afhist[0]/afhist.sum()
        else:
            return np.array(afhist[:-1]*(1-af_binc[:-1])).sum()/afhist.sum()

    from random import choice
    ax = axs[0]
    time_binc = 0.5*(time_bins[1:]+time_bins[:-1])
    sym='o'
    fs = fig_fontsize
    color_count=0
    for subtype, ls in [('patient', '--'), ('any','-')]:
        for toaway, H in [(u'founder = '+('group M' if subtype=='any' else 'subtype'),  away_histogram[subtype]), 
                          (u'founder \u2260 '+('group M' if subtype=='any' else 'subtype'), to_histogram[subtype])]:
            mean_hists = bin_time(H,time_bins)
            div = [get_div(mean_hists[ti]) for ti in range(len(time_bins)-1)]
            # make replicates and calculate bootstrap confidence intervals
            replicates = []
            all_keys = H.keys()
            for ri in xrange(nbs):
                bootstrap_keys = [all_keys[ii] for ii in np.random.randint(len(all_keys), size=len(all_keys))]
                tmp = bin_time({key:H[key] for key in bootstrap_keys}, time_bins)
                replicates.append([get_div(tmp[ti]) for ti in range(len(time_bins)-1)])
            std_dev = np.array(replicates).std(axis=0)
            ax.errorbar(time_binc/365.25, div, std_dev, ls = ls, lw=3, c=colors[color_count])
            ax.plot(time_binc/365.25, div, label = toaway, ls = ls, lw=3, c=colors[color_count]) # plot again with label to avoid error bars in legend
            color_count+=1
    ax.set_ylim([0,0.15])
    ax.set_xlabel('ETI [years]', fontsize=fs)
    ax.set_ylabel('Divergence from founder sequence', fontsize=fs)
    ax.legend(loc=2, fontsize=fs, labelspacing=0)
    add_panel_label(ax, 'A', x_offset=-0.27)
    ax.tick_params(axis='both', labelsize=fs)
    plt.tight_layout() #rect=(0.0, 0.02, 0.98, 0.98), pad=0.05, h_pad=0.5, w_pad=0.4)
    for ext in figtypes:
        fig.savefig(fig_filename+ext)


if __name__=="__main__":
    import argparse
    import matplotlib.pyplot as plt
    import pandas as pd
    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action = 'store_true', help = 'recalculate data')
    params=parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'to_away.pickle'
   
    #patients = ['p1', 'p6'] # other subtypes
    patients = ['p1', 'p2', 'p3','p5', 'p6', 'p8', 'p9','p11'] # all subtypes
    regions = ['genomewide']
    #regions = ['gag', 'pol', 'nef'] #, 'env']
    #regions = ['p24', 'p17'] #, 'RT1', 'RT2', 'RT3', 'RT4', 'PR', 
    #           'IN1', 'IN2', 'IN3','p15', 'vif', 'nef','gp41','gp1201']
    cov_min = 1000
    Sbins = np.array([0,0.03, 0.08, 0.25, 2])
    Sbinc = 0.5*(Sbins[1:]+Sbins[:-1])

    if not os.path.isfile(fn_data) or params.redo:
        data = {}
        for subtype in ['patient', 'any']:
            minor_variants, to_away_divergence, to_away_minor, consensus_distance = \
                collect_to_away(patients, regions, Sbins=Sbins, cov_min=cov_min, subtype = subtype)

            to_away_minor.loc[:,['reversion_spectrum', 'minor_reversion_spectrum']] = \
                            to_away_minor.loc[:,['reversion_spectrum', 'minor_reversion_spectrum']].astype(float)
            add_binned_column(to_away_minor,  [0,1000,2000,4000], 'time')
            data[subtype] = {'minor_variants':minor_variants, 'to_away':to_away_divergence,'to_away_minor':to_away_minor, 
                    'consensus_distance':consensus_distance, 'Sbins':Sbins, 'Sbinc':Sbinc}

        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)


    af_bins = np.linspace(0,1,11)
    af_binc = 0.5*(af_bins[:-1]+af_bins[1:])
    time_bins = np.array([-10,500,1000,1500,2000,2500])
    data['away_histogram'] = {}; data['to_histogram']={}
    for subtype in ['any', 'patient']:
        data['to_histogram'][subtype], data['away_histogram'][subtype] = get_toaway_histograms(subtype, Sc=10)
    data['time_bins']=time_bins
    data['af_bins']=af_bins
    plot_to_away(data, fig_filename=foldername+'to_away')
