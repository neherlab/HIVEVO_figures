import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.samples import all_fragments
from util import store_data, load_data, draw_genome, fig_width, fig_fontsize
import os
from filenames import get_figure_folder

def running_average_masked(obs, ws):
    '''
    calculates a running average
    obs     --  observations
    ws      --  winodw size (number of points to average)
    '''
    try:
        tmp_vals = np.convolve(np.ones(ws, dtype=float), obs*(1-obs.mask), mode='same')
        if len(obs.mask.shape)==0:
            tmp_valid = ws*np.ones_like(tmp_vals)
            # fix the edges. using mode='same' assumes zeros outside the range
            if ws%2==0:
                tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2,ws)
                if ws//2>1:
                    tmp_vals[-ws//2+1:]*=float(ws)/np.arange(ws-1,ws//2,-1.0)
            else:
                tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2+1,ws)
                tmp_vals[-ws//2:]*=float(ws)/np.arange(ws,ws//2,-1.0)
        else:
            tmp_valid = np.convolve(np.ones(ws, dtype=float), (1-obs.mask), mode='same')

        run_avg = np.ma.array(tmp_vals/tmp_valid)
        run_avg.mask = tmp_valid<ws*0.95
    except:
        import pdb; pdb.set_trace()

    return run_avg


def weighted_linear_regression(x,y):
    data = np.array([(tmpx, tmpy) for tmpx, tmpy, m in zip(x,y,y.mask) if not m])
    if len(data)>2:
        weights = data[:,1]+3e-3  #shot noise + sequencing error
        slope = np.sum(data[:,0]*data[:,1]/weights)/np.sum(data[:,0]**2/weights) 
        gof = np.corrcoef(x,y)[0,1]
        return slope, gof
    else:
        return np.nan, np.nan


def plot_evo_rates(data, fig_filename=None, figtypes=['.png', '.svg', '.pdf']):
    '''Plot evolutionary rate in a sliding window'''
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

    # Evolutionary rate
    ax=axs[0]
    HXB2_masked = np.ma.array(data['rates'])
    HXB2_masked.mask = data['rates']<0
    for pi,pcode in enumerate(data['patients']):
        ax.plot(np.arange(HXB2_masked.shape[1])[-HXB2_masked.mask[pi]], 
                 HXB2_masked[pi][-HXB2_masked.mask[pi]], alpha = 0.5, label = pcode)

    ax.plot(np.arange(HXB2_masked.shape[1]),
            np.exp(np.log(HXB2_masked).mean(axis=0)),
            c='k', lw=3,
            label='average')

    ax.yaxis.set_tick_params(labelsize=fs)
    ax.set_ylabel('Evolutionary rate [1/site/year]', fontsize=fs)
    ax.legend(loc='upper left', ncol=3, fontsize=fs-3 ,title='Patients')
    ax.set_ylim([2e-4, 4e-2])
    ax.set_yscale('log')

    # Genome annotations
    ax=axs[1]
    sns.set_style('white')
    from hivevo.hivevo.HIVreference import HIVreference
    refseq = HIVreference('HXB2')
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
        plt.close(fig)
    else:
        plt.ion()
        plt.show()

    print "genome wide variation:", np.std(np.log2(HXB2_masked).mean(axis=0))
    print "position wide variation:", np.mean(np.log2(HXB2_masked+.0001).std(axis=0))



if __name__=="__main__":
    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'evolutionary_rates.pickle'

    patients = ['p1', 'p2', 'p3','p5', 'p6', 'p8', 'p9', 'p11']
    if not os.path.isfile(fn_data):
        print("Regerating plot data")
        rate_or_gof = 0
        window_size=300
        cov_min = 200
        HXB2 = -np.ones((len(patients), 10000), dtype=float)
        evo_rates = {}
        for pi, pcode in enumerate(patients):
            try:
                p = Patient.load(pcode)
            except:
                print "Can't load patient", pcode
            else:
                toHXB2 = p.map_to_external_reference('genomewide')
                aft = p.get_allele_frequency_trajectories('genomewide', cov_min=cov_min)
                aft[aft<0.002]=0
                div_traj = [np.ma.array(af.sum(axis=0) - af[p.initial_indices, np.arange(len(p.initial_indices))], shrink=False) 
                            for af in aft]
                print 'total divergence', zip(p.ysi, [x.sum() for x in div_traj])
                smoothed_divergence = np.ma.array([ running_average_masked(div, window_size) for div in div_traj])
                evo_rates[pcode] =  np.array([weighted_linear_regression(p.ysi, smoothed_divergence[:,i])[rate_or_gof]
                                    for i in xrange(smoothed_divergence.shape[1])])
                HXB2[pi,toHXB2[:,0]] = evo_rates[pcode][toHXB2[:,1]]
        data = {'rates': HXB2, 'patients': patients}
        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)

    plot_evo_rates(data, fig_filename=foldername+'evolutionary_rates')
