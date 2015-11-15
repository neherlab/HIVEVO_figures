'''
script producing the figure comparing evolution at synonymous and nonsynonymous sites
it contains two version of the data collection procedure, which produce identical results
The script ingests data from divergence_diversity_correlation.py
'''
import argparse
import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.samples import all_fragments
from hivevo.hivevo.af_tools import divergence, diversity
from util import store_data, load_data, draw_genome, fig_width, fig_fontsize, add_panel_label, HIVEVO_colormap
from util import boot_strap_patients, replicate_func, add_binned_column
import os
from filenames import get_figure_folder


def collect_data_richard(patients, regions, syn_degeneracy=2):
    '''Collect data for divergence and diversity'''

    syn_divergence = {reg:{p:[] for p in patients} for reg in regions} 
    syn_diversity = {reg:{p:[] for p in patients} for reg in regions}
    nonsyn_divergence = {reg:{p:[] for p in patients} for reg in regions}
    nonsyn_diversity = {reg:{p:[] for p in patients} for reg in regions}
    time_bins = np.array([0, 200, 500, 1000, 1500, 2000, 3000, 5000])

    nbins=10
    sfs_tmin=1000
    sfs = {'syn':np.zeros(nbins, dtype=float), 
           'nonsyn':np.zeros(nbins, dtype='float'),
           'bins':np.linspace(0.01,0.99,nbins+1)}
    time_binc = 0.5*(time_bins[1:]+time_bins[:-1])
    cov_min = 100
    for pi, pcode in enumerate(patients):
        try:
            p = Patient.load(pcode)
        except:
            print "Can't load patient", pcode
        else:
            for region in regions:
                for prot in regions[region]:
                    initial_indices = p.get_initial_indices(prot)
                    aft = p.get_allele_frequency_trajectories(prot, cov_min=cov_min)
                    gaps = p.get_gaps_by_codon(prot, pad=2, threshold=0.05)
                    syn_mask = p.get_syn_mutations(prot)
                    syn_pos = (syn_mask.sum(axis=0)>1)*(gaps==False)
                    nonsyn_pos = (syn_mask.sum(axis=0)<=1)*(p.get_constrained(prot)==False)*(gaps==False)
                    print pcode, prot, syn_pos.sum(), nonsyn_pos.sum()
                    syn_divergence[region][pcode].extend([(t, divergence(af[:,syn_pos], 
                                             initial_indices[syn_pos])) for t,af in zip(p.dsi, aft)])
                    syn_diversity[region][pcode].extend([(t, diversity(af[:,syn_pos])) 
                                                        for t,af in zip(p.dsi, aft)])
                    nonsyn_divergence[region][pcode].extend([(t, divergence(af[:,nonsyn_pos], 
                                             initial_indices[nonsyn_pos])) for t,af in zip(p.dsi, aft)])
                    nonsyn_diversity[region][pcode].extend([(t, diversity(af[:,nonsyn_pos])) 
                                                           for t,af in zip(p.dsi, aft)])

                    syn_derived = syn_mask.copy()
                    syn_derived[initial_indices, np.arange(syn_derived.shape[1])]=False
                    for t,af in izip(p.dsi,aft):
                        if t>sfs_tmin:
                            y,x = np.histogram(af[syn_derived].flatten(), bins=sfs['bins'])
                            sfs['syn']+=y
                    nonsyn_derived = syn_mask==False
                    nonsyn_derived*=(p.get_constrained(prot)==False)*(gaps==False)
                    nonsyn_derived[initial_indices, np.arange(syn_derived.shape[1])]=False
                    for t,af in izip(p.dsi,aft):
                        if t>sfs_tmin:
                            y,x = np.histogram(af[nonsyn_derived], bins=sfs['bins'])
                            sfs['nonsyn']+=y

    for tmp_data in [syn_divergence, syn_diversity, nonsyn_diversity, nonsyn_divergence]:
        for region in regions:
            tmp = np.vstack([np.array(tmp_data[region][p]) for p in patients])
            tmp_clean = tmp[-np.isnan(tmp[:,1]),:]
            y,  x = np.histogram(tmp_clean[:,0],bins = time_bins, weights = tmp_clean[:,1])
            yn, x = np.histogram(tmp_clean[:,0],bins = time_bins)
            tmp_data[region] = {'avg':y/(1e-10+yn), 'bins':time_binc, 'raw':tmp_data[region]}

    data = {'syn_diversity':syn_diversity, 'syn_divergence':syn_divergence,
            'nonsyn_diversity':nonsyn_diversity, 'nonsyn_divergence':nonsyn_divergence,
            'sfs':sfs}

    return data


def collect_data_fabio(patients, regions, cov_min=100, syn_degeneracy=2):
    '''Collect data for divergence and diversity'''
    import pandas as pd
    from itertools import izip

    # Prepare SFS
    nbins=10
    sfs_tmin=1000
    sfs = {'syn': np.zeros(nbins, dtype=float), 
           'nonsyn': np.zeros(nbins, dtype=float),
           'bins': np.linspace(0.01, 0.99, nbins+1),
          }

    # Collect into DataFrame
    data = []
    for pi, pcode in enumerate(patients):
        p = Patient.load(pcode)
        for region, prots in regions.iteritems():
            for prot in prots:
                aft = p.get_allele_frequency_trajectories(prot, cov_min=cov_min)
                initial_indices = p.get_initial_indices(prot)
                gaps = p.get_gaps_by_codon(prot, pad=2, threshold=0.05)

                # Classify syn/nonsyn POSITIONS
                # NOTE: this is not fully correct because some positions (2-fold
                # degenerate) are both syn and nonsyn, but it's close enough
                syn_mask = p.get_syn_mutations(prot)
                syn_sum = syn_mask.sum(axis=0)
                # NOTE: syn_mask == 0 are substitutions, they make up most
                # of the nonsynonymous signal
                pos = {'syn': (syn_sum >= syn_degeneracy) & (-gaps),
                       'nonsyn': (syn_sum <= 1) & (-p.get_constrained(prot)) & (-gaps),
                      }
                
                print pcode, prot, pos['syn'].sum(), pos['nonsyn'].sum()

                # Divergence/diversity
                for t, af in izip(p.dsi, aft):
                    for mutclass, ind in pos.iteritems():
                        data.append({'pcode': pcode,
                                     'time': t,
                                     'region': region,
                                     'protein': prot,
                                     'nsites': ind.sum(),
                                     'mutclass': mutclass,
                                     'divergence': divergence(af[:, ind], initial_indices[ind]),
                                     'diversity': diversity(af[:, ind]),
                                    })


                # Site frequency spectrum
                syn_derived = syn_mask.copy()
                syn_derived[initial_indices, np.arange(syn_derived.shape[1])] = False
                nonsyn_derived = (-syn_mask) & (-p.get_constrained(prot)) & (-gaps)
                nonsyn_derived[initial_indices, np.arange(syn_derived.shape[1])] = False

                for t,af in izip(p.dsi,aft):
                    if t < sfs_tmin:
                        continue

                    sfs['syn'] += np.histogram(af[syn_derived], bins=sfs['bins'])[0]
                    sfs['nonsyn'] += np.histogram(af[nonsyn_derived], bins=sfs['bins'])[0]


    data = pd.DataFrame(data)
    data['divergence'] = data['divergence'].astype(float)
    data['diversity'] = data['diversity'].astype(float)
    return {'divdiv':data, 'sfs':sfs}


def plot_divdiv(data, fig_filename=None, figtypes=['.png', '.svg', '.pdf']):
    '''
    plot divergence and diversity of synonymous and nonsynonymous mutations
    includes:
        - a panel that compares syn diversity/nonsyn divergence
    '''
    n_bootstrap=50
    ####### plotting ###########
    import seaborn as sns
    from matplotlib import pyplot as plt
    plt.ion()
    sns.set_style('darkgrid')
    figpath = 'figures/'
    fs=fig_fontsize
    fig_size = (fig_width, 1.0*fig_width)
    cols = HIVEVO_colormap()

    fig, axs = plt.subplots(2, 2,figsize=fig_size)
    divdiv = data['divdiv']
    # in rough the order in which they most dominantly appear in the plot
    regions = ['envelope', 'accessory', 'structural','enzymes']
    time_bins = np.array([0, 200, 500, 1000, 1500, 2000, 4000])
    time_binc = 0.5*(time_bins[:-1]+time_bins[1:])
    add_binned_column(divdiv, time_bins, 'time')

    # map the regions to rough genomic order to match the genome color map in panel C
    colors = {reg:c for reg, c in zip(regions, [cols(x) for x in [0.66, 0.99, 0.01, 0.33]])}
    def get_time_bin_mean(df):
        return df.loc[:,['time_bin', 'diversity', 'divergence']].groupby(by=['time_bin'], as_index=False).mean()
    def label_func(mutclass, region, divordiv):  # assign labels to Panels A and B separately to make a combinatorial legend (regions vs syn/nonsyn)
        if divordiv=='divergence' and mutclass=='nonsyn':
            return region
        elif divordiv=='diversity' and region=='accessory':
            return mutclass
        else:
            return None

    ########## panel A and B #####################
    csv_out = open(fig_filename+'_AB.tsv', 'w')
    for ax, dtype in izip(axs[0,:], ['divergence', 'diversity']):
        add_panel_label(ax, 'A' if dtype=='divergence' else 'B', x_offset = -0.3)
        for mutclass in ['nonsyn', 'syn']:
            for region in regions:
                ind = (divdiv.loc[:,'region']==region) & (divdiv.loc[:,'mutclass']==mutclass)
                tmp = divdiv.loc[ind,['time_bin', 'diversity', 'divergence', 'pcode']]
                avg_divdiv = get_time_bin_mean(tmp)
                bs = boot_strap_patients(tmp, eval_func = get_time_bin_mean, n_bootstrap=n_bootstrap)
                # plot the same line with and without error bars, labels for legend without
                ax.plot(time_binc/365.25, avg_divdiv.loc[:,dtype], ls='-' if mutclass=='nonsyn' else '--', 
                            c=colors[region], lw=3, label=label_func(mutclass, region, dtype))
                ax.errorbar(time_binc/365.25, avg_divdiv.loc[:,dtype], replicate_func(bs, dtype, np.std, bin_index='time_bin'),
                            ls='-' if mutclass=='nonsyn' else '--', c=colors[region], lw=3)
                csv_out.write('\t'.join(map(str,[dtype, mutclass, region]+list(avg_divdiv.loc[:,dtype])))+'\n')

        ax.legend(loc=2, fontsize=fs-1, numpoints=2, labelspacing = 0)
        ax.set_xticks([0,2,4,6,8])
        if dtype=='divergence':
            ax.set_yticks([0,.02,.04])
            ax.set_ylim([0,.048])
        else:
            ax.set_yticks([0,.01,.02])
            ax.set_ylim([0,.028])            
        ax.set_xlim([0,8.5])
        ax.set_ylabel(dtype)
        ax.tick_params(labelsize=fs-2)
        ax.set_xlabel('ETI [years]', fontsize=fs)
    csv_out.close()

    ########## panel C: anti correlation of syn diversity and nonsyn divergence #############
    csv_out = open(fig_filename+'_C.tsv', 'w')
    (avg_nonsyn_divg, avg_nonsyn_divs, avg_syn_divs) = data['divdiv_corr']
    ax = axs[1,0]
    add_panel_label(ax, 'C', x_offset = -0.3)
    x_data, y_data = avg_nonsyn_divg[::500], avg_syn_divs[::500]
    ax.scatter(x_data, y_data, c=[cols(p) for p in np.linspace(0,1,len(x_data))], s=50)
    csv_out.write('\t'.join(map(str, ["nonsyn_divergence"]+list(x_data)))+'\n')
    csv_out.write('\t'.join(map(str, ["syn_diversity"]+list(y_data)))+'\n')
    csv_out.close()

    ax.set_xlabel('nonsyn divergence', fontsize = fig_fontsize)
    ax.set_ylabel('syn diversity', fontsize = fig_fontsize)
    ax.set_ylim([0,0.028])
    ax.set_xlim([0,0.012])
    ax.set_xticks([0, 0.005,0.01])
    ax.set_yticks([0, 0.01, 0.02])
    ax.tick_params(labelsize=fig_fontsize-2)


    ########## sfs in panel D ##############
    csv_out = open(fig_filename+'_D.tsv', 'w')
    sfs=data['sfs']
    ax = axs[1,1]
    add_panel_label(ax, 'D', x_offset = -0.3)
    colors = sns.color_palette(n_colors=2)
    binc = binc = 0.5*(sfs['bins'][1:]+sfs['bins'][:-1])
    ax.bar(binc-0.045, sfs['syn']/np.sum(sfs['syn']),width = 0.04, label='syn', color=colors[0])
    ax.bar(binc, sfs['nonsyn']/np.sum(sfs['nonsyn']),width = 0.04, label='nonsyn', color=colors[1])
    csv_out.write('\t'.join(map(str, ["bin_centers"]+list(binc)))+'\n')
    csv_out.write('\t'.join(map(str, ["sfs_nonsyn"]+list(sfs['nonsyn']/np.sum(sfs['nonsyn']))))+'\n')
    csv_out.write('\t'.join(map(str, ["sfs_syn"]+list(sfs['syn']/np.sum(sfs['syn']))))+'\n')
    csv_out.close()
    ax.set_ylim([0.005,2.0])
    ax.set_yscale('log')
    ax.set_xlabel('Frequency',fontsize=fs)
    ax.set_ylabel('Fractions of SNPs',fontsize=fs)
    ax.legend(loc=1, fontsize=fs-2)
    ax.tick_params(labelsize=fig_fontsize-2)

    # finalize and save the figure
    plt.tight_layout(rect=(0.0, 0.02, 0.98, 0.98), pad=0.05, h_pad=0.5, w_pad=0.4)
    if fig_filename is not None:
        for ext in figtypes:
            fig.savefig(fig_filename+ext)
    else:
        plt.ion()
        plt.show()

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Make figure for divergence and diversity")
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    params = parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn2_data = fn_data + 'divdiv_correlation.pickle'
    fn_data = fn_data + 'syn_nonsyn_divergence.pickle'

    if not os.path.isfile(fn_data) or params.redo:
        patients = ['p1', 'p2', 'p3','p5', 'p6', 'p8', 'p9', 'p10','p11']
        regions = {'structural':['gag'], #['p17', 'p24'],
                    'enzymes':  ['pol'], #['PR', 'RT', 'p15', 'IN'],
                    'accessory': ['vif', 'nef', 'vpr', 'vpu', 'tat', 'rev'],
                    'envelope': ['env'] #['gp41', 'gp120'],
                    }
        # NOTE: these two give the same result, good
        data = collect_data_fabio(patients, regions)
        #data = collect_data_richard(patients, regions)
        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)
    # this load additional data produced by script divergence_diversity_correlation
    data['divdiv_corr'] = load_data(fn2_data)

    plot_divdiv(data, fig_filename = foldername+'divdiv')
