import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.HIVreference import HIVreference
from hivevo.hivevo.af_tools import divergence
from util import store_data, load_data, fig_width, fig_fontsize, add_binned_column
from util import boot_strap_patients, replicate_func
import os
from filenames import get_figure_folder

def collect_to_away(patients, regions, Sbins=[0,0.02, 0.08, 0.25, 2], cov_min=1000):
    hxb2 = HIVreference(refname='HXB2')
    good_pos_in_reference = hxb2.get_ungapped(threshold = 0.05)
    minor_variants = []
    to_away_divergence = []
    consensus_distance = {}

    # determine genome wide fraction of alleles above a threshold
    for pi, pcode in enumerate(patients):
        try:
            p = Patient.load(pcode)
        except:
            print "Can't load patient", pcode
        else:
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
                    
                    for sbin, sites in enumerate(clean_entropy_bins):
                        minor_variants.append({'pcode':pcode,'region':region,'time':t,'S_bin':sbin,
                                            'af_away_minor':  np.mean(clean_minor[sites&clean_away]), 
                                            'af_away_derived':np.mean(clean_derived[sites&clean_away]),
                                            'af_to_minor':    np.mean(clean_minor[sites&(~clean_away)]), 
                                            'af_to_derived':  np.mean(clean_derived[sites&(~clean_away)])})

                    clean_reversion = clean_af[clean_consensus,np.arange(clean_consensus.shape[0])]*(~clean_away)
                    clean_total_divergence = clean_af.sum(axis=0) - clean_af[clean_ancestral,np.arange(clean_ancestral.shape[0])]
                    to_away_divergence.append({'pcode':pcode,'region':region,'time':t,
                                        'reversion':np.mean(clean_reversion), 
                                        'divergence':np.mean(clean_total_divergence)})

    return pd.DataFrame(minor_variants), pd.DataFrame(to_away_divergence), consensus_distance

def plot_to_away(data, fig_filename = None, figtypes=['.png', '.svg', '.pdf']):
    ####### plotting ###########
    import seaborn as sns
    from matplotlib import pyplot as plt
    plt.ion()
    sns.set_style('darkgrid')
    figpath = 'figures/'
    fs=fig_fontsize
    fig_size = (fig_width, 0.8*fig_width)
    fig, axs = plt.subplots(1, 1, figsize=fig_size)

    ax=axs
    Sbins = np.array([0,0.02, 0.08, 0.25, 2])
    Sbinc = 0.5*(Sbins[1:]+Sbins[:-1])
    mv = data['minor_variants']
    mv.loc[:,['af_away_minor', 'af_away_derived', 'af_to_minor', 'af_to_derived']] = \
            mv.loc[:,['af_away_minor', 'af_away_derived', 'af_to_minor', 'af_to_derived']].astype(float)
    def get_Sbin_mean(df):
        return df.groupby(by=['S_bin'], as_index=False).mean()
    mean_to_away =get_Sbin_mean(mv)
    bs = boot_strap_patients(mv, eval_func=get_Sbin_mean, 
                             columns=['af_away_minor', 'af_away_derived', 'af_to_minor', 'af_to_derived', 'S_bin'])

    col = 'af_away_minor'
    ax.errorbar(Sbinc, mean_to_away.loc[:,col], replicate_func(bs, col, np.std), label = 'equal subtype')
    col = 'af_to_minor'
    ax.errorbar(Sbinc, mean_to_away.loc[:,col], replicate_func(bs, col, np.std), label = 'not equal subtype')
    #ax.plot(Sbinc, mean_to_away.loc[mean_to_away.loc[:,'away']==True,'af_derived'] , label = 'away, der')
    #ax.plot(Sbinc, mean_to_away.loc[mean_to_away.loc[:,'away']==False,'af_derived'], label = 'to, der')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('minor SNV frequencies')
    ax.set_xlabel('entropy [bits]')
    ax.set_xlim([0.005,2])
    ax.legend(loc = 'lower right')

    to_away = data['to_away']
    time_bins = np.array([0,500,1000,1500, 2500, 3500])
    binc = 0.5*(time_bins[1:]+time_bins[:-1])
    add_binned_column(to_away, time_bins, 'time')
    to_away.loc[:,['reversion', 'divergence']] = \
            to_away.loc[:,['reversion', 'divergence']].astype(float)

    def get_time_bin_means(df):
        return df.loc[:,['divergence', 'reversion','time_bin']].groupby(by=['time_bin'], as_index=False).mean()
    rev_div = get_time_bin_means(to_away)
    bs = boot_strap_patients(to_away, get_time_bin_means, columns = ['reversion','divergence','time_bin'])
    reversion_std = replicate_func(bs, 'reversion', np.std)
    total_div_std = replicate_func(bs, 'divergence', np.std)
    fraction = rev_div.loc[:,'reversion']/rev_div.loc[:,'divergence']
    print "Reversions:\n", rev_div.loc[:,'reversion']
    print "Divergence:\n", rev_div.loc[:,'divergence']
    print "Fraction:"
    for frac, total, num_std, denom_std in zip(fraction, rev_div.loc[:,'divergence'],reversion_std, total_div_std):
        print frac, '+/-', np.sqrt(num_std**2/total**2 + denom_std**2*frac**2/total**2)
    #print reversion_std,total_div_std
    print "Consensus!=Founder:",np.mean(data['consensus_distance'].values())

    plt.tight_layout(rect=(0.0, 0.02, 0.98, 0.98), pad=0.05, h_pad=0.5, w_pad=0.4)
    if fig_filename is not None:
        for ext in figtypes:
            fig.savefig(fig_filename+'_sfs'+ext)
    else:
        plt.ion()
        plt.show()


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
    
    if not os.path.isfile(fn_data) or params.redo:
        #patients = ['p1', 'p6'] # other subtypes
        patients = ['p2', 'p3','p5','p8', 'p9','p10', 'p11'] # subtype B
        regions = ['genomewide']
        #regions = ['gag', 'pol', 'env']
        #regions = ['p24', 'p17'] #, 'RT1', 'RT2', 'RT3', 'RT4', 'PR', 
        #           'IN1', 'IN2', 'IN3','p15', 'vif', 'nef','gp41','gp1201']
        cov_min = 1000
        Sbins = np.array([0,0.02, 0.08, 0.25, 2])
        Sbinc = 0.5*(Sbins[1:]+Sbins[:-1])

        minor_variants, to_away_divergence, consensus_distance = \
            collect_to_away(patients, regions, Sbins=Sbins, cov_min=cov_min)

        data = {'minor_variants':minor_variants, 'to_away':to_away_divergence, 
                'consensus_distance':consensus_distance, 'Sbins':Sbins, 'Sbinc':Sbinc}
        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)

plot_to_away(data, fig_filename=foldername+'to_away')
#
