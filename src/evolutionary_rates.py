# Modules
import os, sys
import numpy as np
from itertools import izip
from matplotlib import pyplot as plt

from hivevo.hivevo.patients import Patient
from hivevo.hivevo.samples import all_fragments
from hivevo.hivevo.HIVreference import HIVreference, HIVreferenceAminoacid

from util import store_data, load_data, draw_genome, fig_width, fig_fontsize
from filenames import get_figure_folder



# Functions
def running_average_masked(obs, ws, min_valid_fraction=0.95):
    '''
    calculates a running average via convolution, fixing the edges
    obs     --  observations (a masked array)
    ws      --  window size (number of points to average)
    '''
    #tmp_vals = np.convolve(np.ones(ws, dtype=float), obs*(1-obs.mask), mode='same')
    tmp_vals = np.convolve(np.ones(ws, dtype=float), obs.filled(0), mode='same')

     # if the array is not masked, edges needs to be explictly fixed due to smaller counts
    if len(obs.mask.shape) == 0:
        tmp_valid = ws*np.ones_like(tmp_vals)
        # fix the edges. using mode='same' assumes zeros outside the range
        if ws%2==0:
            tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2,ws)
            if ws//2>1:
                tmp_vals[-ws//2+1:]*=float(ws)/np.arange(ws-1,ws//2,-1.0)
        else:
            tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2+1,ws)
            tmp_vals[-ws//2:]*=float(ws)/np.arange(ws,ws//2,-1.0)

    # if the array is masked, then we get the normalizer from counting the unmasked values
    else:
        tmp_valid = np.convolve(np.ones(ws, dtype=float), (1-obs.mask), mode='same')

    run_avg = np.ma.array(tmp_vals / tmp_valid)
    run_avg.mask = tmp_valid < ws * min_valid_fraction

    return run_avg


def weighted_linear_regression(x, y):
    data = np.array([(tmpx, tmpy) for tmpx, tmpy, m in zip(x,y,y.mask) if not m])
    if len(data)>2:
        weights = data[:,1] + 3e-3  #shot noise + sequencing error
        slope = np.sum(data[:,0]*data[:,1]/weights)/np.sum(data[:,0]**2/weights)
        gof = np.corrcoef(x,y)[0,1]
        return slope, gof
    else:
        return np.nan, np.nan


def get_divergence_trajectory(p, cov_min=100, sequence_type='nuc',
                              only_substitutions=False):
    '''Get divergence in time for a patient'''
    import numpy as np

    div_traj = np.ma.masked_all((len(p.samples),
                                 len(p.get_initial_sequence('genomewide'))),
                                float)

    if sequence_type == 'nuc':
        regions = ['genomewide']
    else:
        regions = ['p17', 'p24', 'p6', 'p7', 'p2', 'p1',
                   'PR', 'RT', 'p15', 'IN',
                   'vif', 'vpu', 'vpr',
                   'tat',
                   'gp41', 'gp120', 'nef']

    for region in regions:
        print region
        if region == 'genomewide':
            reg_coo = np.arange(div_traj.shape[-1])
        else:
            reg_coo = list(p.annotation[region])[::3]

        aft = p.get_allele_frequency_trajectories(region, cov_min=cov_min,
                                                  type=sequence_type)
        aft[aft < 0.002] = 0

        ii = p.get_initial_indices(region, type=sequence_type)
        ii2 = np.zeros((aft.shape[1], aft.shape[2]), bool)
        ii2[ii, np.arange(len(ii))] = True
        for it, af in enumerate(aft):
            if not only_substitutions:
                # NOTE: af.sum(axis=0) is always one or masked
                d = af.sum(axis=0) - af[ii, np.arange(len(ii))]

            else:
                # NOTE: this criterion for substitutions is quite arbitrary
                pos_subst = (af > 0.9) & (~ii2)
                d = np.zeros_like(ii2, float)
                d[pos_subst] = af[pos_subst]
                d = d.sum(axis=0)
                d = np.ma.array(d, mask=af.mask.any(axis=0))

            d = np.ma.array(d, shrink=False)
            div_traj[it, reg_coo] = d

    return div_traj


def plot_evo_rates(data, fig_filename=None, figtypes=['.png', '.svg', '.pdf'],
                   include_substitutions=False,
                   refname='HXB2'):
    '''Plot evolutionary rate in a sliding window'''
    print 'Plot evolutionary rates'
    import seaborn as sns
    from matplotlib import pyplot as plt

    plt.ioff()
    sns.set_style('darkgrid')
    figpath = 'figures/'
    fs=fig_fontsize

    fig, axs = plt.subplots(2, 1,
                            sharex=True,
                            figsize=(fig_width, 0.8*fig_width),
                            gridspec_kw={'height_ratios':[8, 1]})

    # plot the divergence rates for each of the patients
    ax=axs[0]
    ref_masked = np.ma.array(data['rates'])
    ref_masked.mask = data['rates']<0
    for pi,pcode in enumerate(data['patients']):
        ax.plot(np.arange(ref_masked.shape[1])[~ref_masked.mask[pi]],
                 ref_masked[pi][~ref_masked.mask[pi]], alpha = 0.5, label = pcode)

    # plot the average
    ax.plot(np.arange(ref_masked.shape[1]),
            np.exp(np.log(ref_masked).mean(axis=0)),
            c='k', lw=3,
            label='average')

    # plot the average of substitutions only
    if include_substitutions:
        ref_subst = np.ma.array(data['rates_substitutions'])
        ref_subst.mask = data['rates_substitutions']<0
        ax.plot(np.arange(ref_subst.shape[1]),
                np.exp(np.log(ref_subst).mean(axis=0)),
                c=[0.4] * 3,
                lw=1.5,
                label='substitutions')

    ax.yaxis.set_tick_params(labelsize=fs)
    ax.set_ylabel('Divergence rate [1/site/year]', fontsize=fs)

    ax.legend(loc='upper left', ncol=3, fontsize=fs-3 ,title='Patients')
    ax.set_ylim([2e-4, 8e-2])
    ax.set_yscale('log')

    # add genome annotation
    ax=axs[1]
    sns.set_style('white')
    from hivevo.hivevo.HIVreference import HIVreference
    refseq = HIVreference(refname)
    genome_annotations = {name:refseq.annotation[name]
                          for name in ["LTR5'", 'gag', 'pol', 'vif','vpr','vpu',
                                       'gp120', 'RRE', 'gp41', 'nef', "LTR3'"]}
    draw_genome(ax, genome_annotations,fs=7)
    ax.set_yticks([])
    ax.set_xlabel('Position [bp]', fontsize=fs)

    ax.xaxis.set_tick_params(labelsize=fs)
    # Final touches
    plt.tight_layout(rect=(0.0, 0.02, 0.98, 0.98), pad=0.1, h_pad=0.5, w_pad=0.4)

    if fig_filename is not None:
        for ext in figtypes:
            fig.savefig(fig_filename+ext)
        #plt.close(fig)
    else:
        plt.ion()
        plt.show()

    # output statistics
    print "genome wide variation:", np.std(np.log2(ref_masked).mean(axis=0))
    print "position wide variation:", np.mean(np.log2(ref_masked+.0001).std(axis=0))



# Script
if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    parser.add_argument('--type', choices=['nuc', 'aa'], default='nuc',
                        help='Sequence type (nuc or aa)')
    parser.add_argument('--reference', choices=['HXB2', 'NL4-3'], default='HXB2',
                        help='Reference')
    params = parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'evolutionary_rates'
    if params.reference != 'HXB2':
        fn_data = fn_data + '_'+params.reference
    if params.type == 'aa':
        fn_data = fn_data + '_aa'
    fn_data = fn_data + '.pickle'

    patients = ['p1', 'p2', 'p3','p5', 'p6', 'p8', 'p9', 'p11']

    rate_or_gof = 0
    window_size = 300
    cov_min = 200

    if not os.path.isfile(fn_data) or params.redo:
        print("Regenerating plot data")
        cats = [{'name': 'total', 'only_substitutions': False},
                {'name': 'substitutions', 'only_substitutions': True},
               ]
        ref = {key: -np.ones((len(patients), 10000), dtype=float) for key in ['total', 'substitutions']}
        evo_rates = {key: {} for key in ref}
        for pi, pcode in enumerate(patients):
            p = Patient.load(pcode)
            to_ref = p.map_to_external_reference('genomewide')

            for cat in cats:
                div_traj = get_divergence_trajectory(p, cov_min=cov_min,
                                                     sequence_type=params.type,
                                                     only_substitutions=cat['only_substitutions'])

                print (pcode, cat['name']+' divergence',
                       zip(np.round(p.ysi),
                           [[np.round(x[x<th].sum()) for th in [.1, .5, 0.95, 1.0]] for x in div_traj]))

                if params.type == 'nuc':
                    min_valid_fraction = 0.95
                else:
                    # Two out of three are masked by design
                    min_valid_fraction = 0.30
                smoothed_divergence = np.ma.array([running_average_masked(div, window_size,
                                                                          min_valid_fraction=min_valid_fraction)
                                                   for div in div_traj])

                evo_rates[cat['name']][pcode] = \
                        np.array([weighted_linear_regression(p.ysi, smoothed_divergence[:,i])[rate_or_gof]
                                  for i in xrange(smoothed_divergence.shape[1])])

                ref[cat['name']][pi, to_ref[:,0]] = evo_rates[cat['name']][pcode][to_ref[:,1]]

        data = {'rates': ref['total'],
                'rates_substitutions': ref['substitutions'],
                'patients': patients,
               }

        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)

    fig_filename = foldername+'evolutionary_rates_withsubst'
    if params.reference != 'HXB2':
        fig_filename = fig_filename + '_' + params.reference
    if params.type == 'aa':
        fig_filename = fig_filename + '_aa'
    plot_evo_rates(data, fig_filename=fig_filename)

    # calculate the overall correlation with diversity
    ref = HIVreference(refname=params.reference, subtype='any')
    ref_masked = np.ma.array(data['rates'])
    ref_masked.mask = data['rates']<0
    avg_log_div = np.log(ref_masked).mean(axis=0)
    from scipy.stats import spearmanr
    print ("Correlation between groupM entropy and divergence,",
           spearmanr(avg_log_div[:len(ref.entropy)], running_average_masked(np.ma.array(ref.entropy),300)))

