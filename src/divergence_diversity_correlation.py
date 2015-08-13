'''
script that prepares data for the figure comparing evolution at synonymous and nonsynonymous sites
'''
import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.samples import all_fragments
from util import store_data, load_data, draw_genome, fig_width, fig_fontsize, patients, HIVEVO_colormap
from evolutionary_rates import running_average_masked, weighted_linear_regression
import os
from filenames import get_figure_folder
import matplotlib.pyplot as plt


def get_divergence_trajectory(p, aft=None):
    '''Get divergence in time for a patient'''
    if aft is None:
        aft = p.get_allele_frequency_trajectories('genomewide', cov_min=cov_min)
        aft[aft < 0.002] = 0

    div_traj = []
    ii = p.initial_indices
    ii2 = np.zeros((aft.shape[1], aft.shape[2]), bool)
    ii2[ii, np.arange(len(ii))] = True
    for af in aft:
        # NOTE: af.sum(axis=0) is always one or masked
        d = af.sum(axis=0) - af[ii, np.arange(len(ii))]
        d = np.ma.array(d, shrink=False) 
        div_traj.append(d)

    return np.ma.array(div_traj)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action = 'store_true', help = 'recalculate data')
    params=parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'divdiv_correlation.pickle'

    window_size=1000
    cov_min = 200
    # translated regions that tile the HIV genome to a large extend
    regions = ['gag', 'pol','vif', 'vpu', 'vpr', 'nef', 'env']
    if not os.path.isfile(fn_data) or params.redo:
        print("Regenerating plot data")

        # prepare arrays to accumulate divergence and diversity data
        # all are -1, stuff that stays negative will be eventually masked
        HXB2_syn_divs = -np.ones((len(patients), 10000), dtype=float)
        HXB2_nonsyn_divs = -np.ones((len(patients), 10000), dtype=float)
        HXB2_nonsyn_divg = -np.ones((len(patients), 10000), dtype=float)
        for pi, pcode in enumerate(patients):
            print("patient:",pcode)
            p = Patient.load(pcode)
            for region in regions:
                # map each regional alignment to HXB2, exclude regions gapped in the global alignmnt
                toHXB2 = p.map_to_external_reference(region)
                aft = p.get_allele_frequency_trajectories(region, cov_min=cov_min)
                initial_indices = p.get_initial_indices(region)
                diversity = (aft*(1-aft)).sum(axis=1)[p.ysi>4].mean(axis=0)
                divergence = (1-aft[:,initial_indices,np.arange(len(initial_indices))])[-1]/p.ysi[-1]
                gaps = p.get_gaps_by_codon(region, pad=2, threshold=0.05)
                syn_mask = p.get_syn_mutations(region)
                syn_pos = (syn_mask.sum(axis=0)>1)*(gaps==False)*(~divergence.mask)
                #nonsyn_pos = (syn_mask.sum(axis=0)<=1)*(p.get_constrained(region)==False)*(gaps==False)
                nonsyn_pos = (syn_mask.sum(axis=0)<=1)*(gaps==False)*(~divergence.mask)

                divs_syn = -np.ones_like(diversity)
                divs_syn[syn_pos] = diversity[syn_pos]
                divs_nonsyn = -np.ones_like(diversity)
                divs_nonsyn[nonsyn_pos] = diversity[nonsyn_pos]
                divg_nonsyn = -np.ones_like(divergence)
                divg_nonsyn[nonsyn_pos] = divergence[nonsyn_pos]

                HXB2_syn_divs[pi,toHXB2[:,0]] = divs_syn[toHXB2[:,2]]
                HXB2_nonsyn_divs[pi,toHXB2[:,0]] = divs_nonsyn[toHXB2[:,2]]
                HXB2_nonsyn_divg[pi,toHXB2[:,0]] = divg_nonsyn[toHXB2[:,2]]
                
        # all HXB2 arrays now contain data where appropriate
        # negative values are masked (i.e. positions that are never syn)
        # and we take the average over patients
        HXB2_syn_divs = np.ma.array(HXB2_syn_divs)
        HXB2_syn_divs.mask = HXB2_syn_divs<0
        avg_HXB2_syn_divs = HXB2_syn_divs.mean(axis=0)

        HXB2_nonsyn_divs = np.ma.array(HXB2_nonsyn_divs)
        HXB2_nonsyn_divs.mask = HXB2_nonsyn_divs<0
        avg_HXB2_nonsyn_divs = HXB2_nonsyn_divs.mean(axis=0)

        HXB2_nonsyn_divg = np.ma.array(HXB2_nonsyn_divg)
        HXB2_nonsyn_divg.mask = HXB2_nonsyn_divg<0
        avg_HXB2_nonsyn_divg = HXB2_nonsyn_divg.mean(axis=0)

        # determine the running average over positions
        avg_syn_divs = running_average_masked(avg_HXB2_syn_divs, window_size, 0.15)
        avg_nonsyn_divs = running_average_masked(avg_HXB2_nonsyn_divs, window_size, 0.3)
        avg_nonsyn_divg = running_average_masked(avg_HXB2_nonsyn_divg, window_size, 0.3)

        store_data((avg_nonsyn_divg, avg_nonsyn_divs, avg_syn_divs), fn_data)
    else:
        (avg_nonsyn_divg, avg_nonsyn_divs, avg_syn_divs) = load_data(fn_data)


    plt.ion()
    plt.plot(avg_nonsyn_divg)
    plt.plot(avg_syn_divs)
    plt.plot(avg_nonsyn_divs)

    fig,axs = plt.subplots(1,2,sharey=True, figsize = (fig_width, 0.6*fig_width))
    cols = HIVEVO_colormap()
    ax = axs[0]
    ax.scatter(avg_nonsyn_divg[::(window_size/2)], avg_nonsyn_divs[::(window_size/2)],  label='nonsynonymous',
                   c=[cols(p) for p in np.linspace(0,1,len(avg_nonsyn_divg[::(window_size/2)]))], s=50)
    ax.set_ylabel('diversity', fontsize = fig_fontsize)
    ax.set_xlabel('nonsyn divergence', fontsize=fig_fontsize)
    ax.set_xlim([0,0.012])
    ax.set_xticks([0, 0.005, 0.01])
    ax.tick_params(labelsize=fig_fontsize)
    #ax.legend(loc=2)

    ax = axs[1]
    ax.scatter(avg_nonsyn_divg[::(window_size/2)], avg_syn_divs[::(window_size/2)], 
                   c=[cols(p) for p in np.linspace(0,1,len(avg_nonsyn_divg[::(window_size/2)]))], s=50, label='synonymous')
    ax.set_xlabel('nonsyn divergence', fontsize = fig_fontsize)
    ax.set_xlim([0,0.012])
    ax.set_xticks([0, 0.005,0.01])
    ax.tick_params(labelsize=fig_fontsize)
    ax.set_ylim([0,0.03])
    plt.tight_layout()
