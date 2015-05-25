import numpy as np
import sys,os
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.samples import all_fragments
from hivevo.hivevo.af_tools import LD as LDfunc
from hivwholeseq.filenames import root_data_folder
from filenames import get_figure_folder
from util import store_data, load_data


def control_LD(PCR='PCR1', fragment='F3', var_min = 0.2):
    control_fn = (root_data_folder+'specific/PCR_recombination/'+
              'RNA_mix'+PCR+'_cocounts_'+fragment+'.pickle')
    def load_cocounts(PCR, fragment):
        '''Load cocounts from file'''
        import cPickle as pickle
        with open(control_fn, 'rb') as f:
            return pickle.load(f)

    # Globals
    samplenames = {'PCR1': 'MIX1_new_PCR1',
                   'PCR2': 'MIX1_new_PCR2'}
    refnames = ['LAI-III', '38304']
    # NOTE: F1 is not very good (some shorter PCR products)
    # F2 and F3 have no indels between the references, which makes life easier
    # F4 and F5 have indels

    data = load_cocounts(PCR, fragment)
    acc = data['cocounts']
    af1p = np.array(acc[:,:,np.arange(acc.shape[2]), np.arange(acc.shape[3])].sum(axis=1),dtype=float)
    af1p = af1p/(1e-10+af1p.sum(axis=0))
    variable_sites = np.sum(af1p**2, axis=0)<1.-var_min
    reduced_af1p = af1p[:,variable_sites]
    positions = np.where(variable_sites)[0]
    n_variable_sites = reduced_af1p.shape[-1]
    reduced_af2p = np.zeros((acc.shape[0], acc.shape[1], n_variable_sites, n_variable_sites), dtype = float)
    for di,dsite in enumerate(positions):
        reduced_af2p[:,:,di,:] = acc[:,:,dsite,variable_sites]

    reduced_acc_cov = np.array(reduced_af2p.sum(axis=1).sum(axis=0), dtype=int)
    reduced_af2p /= (1e-10+reduced_acc_cov) 
    return (positions, reduced_af2p, reduced_acc_cov, af1p[:,variable_sites])

def plot_LD(data, fig_filename = None):
    import seaborn as sns
    from matplotlib import pyplot as plt
    plt.ion()
    sns.set_style('darkgrid')
    figpath = 'figures/'
    fs=16
    fig_size = (5.5, 4.3)

    binc = data['binc']

    for LD, label, measure in [(data['Dp'], "linkage disequilibrium D'", 'Dp'),
                       (data['LDrsq'], 'linkage disequilibrium r^2', 'rsq')]:
        fig = plt.figure(figsize=fig_size)
        ax = plt.subplot(111)
        for frag in all_fragments:
            y = LD[frag]
            plt.plot(binc, y, label=frag, lw=3)
        plt.plot(binc, LD['PCR1'], label="control", lw=2, ls='--', c='k')
    
        ax.set_xticks(range(0,401,100))
        ax.set_yticks(np.arange(0.,1.01,0.2))
        for item in ax.get_xticklabels() + ax.get_yticklabels():
            item.set_fontsize(fs)
        plt.ylabel(label,fontsize = fs)
        plt.xlabel("distance [bp]",fontsize = fs)
        plt.legend(loc=(0.5,0.45), ncol=2, fontsize = fs-3, columnspacing=1)
        plt.tight_layout(rect=(0, 0, 0.98, 1))
        if fig_filename is not None:
            for ext in ['.pdf','.svg', '.png']:
                plt.savefig(fig_filename+"_"+measure+ext)

if __name__=="__main__":
    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'LD.pickle'
    patients = ['p' +str(i) for i in xrange(1,12) if i not in [4,7]]

    if not os.path.isfile(fn_data):
        print("Regerating plot data")
        dmin = 40
        dmin_pad = 200
        var_min = 0.2
        cov_min = 200
        LD_vs_distance = {}
        Dp_vs_distance = {}
        bins = np.arange(0,401,40)
        binc = (bins[:-1]+bins[1:])*0.5
        for frag in all_fragments:
            if frag not in ['F'+str(i) for i in xrange(1,7)]:
                continue
            dists = []
            weights_LD = []
            weights_Dp = []
            for pcode in patients:
                try:
                    p = Patient.load(pcode)
                    depth = p.get_fragment_depth(pad=False, limit_to_dilution=False)
                    depth_pad = p.get_fragment_depth(pad=True, limit_to_dilution=False)
                except:
                    print "Can't load patient", pcode
                else:
                    for si, sample in enumerate(p.samples):
                        if depth[si][all_fragments.index(frag)]>dmin \
                            or depth_pad[si][all_fragments.index(frag)]>dmin_pad:
                            positions, af2p, cov, af1p = sample.get_pair_frequencies(frag, var_min=var_min)
                            if positions is None:
                                continue
                            LD, Dp, p12 =  LDfunc(af2p, af1p, cov, cov_min=100)
                            X,Y = np.meshgrid(positions, positions)
                            np.fill_diagonal(cov, 0)
                            dists.extend(np.abs(X-Y)[cov>=cov_min])
                            weights_LD.extend(LD[cov>=cov_min])
                            weights_Dp.extend(Dp[cov>=cov_min])
                            print pcode, si, frag, " # of positions:", len(positions), 'depth:', depth[si][all_fragments.index(frag)]
                        else:
                            print pcode, si, frag, "insufficient depth:", depth[si][all_fragments.index(frag)], depth_pad[si][all_fragments.index(frag)]

            yn,xn = np.histogram(dists, bins = bins)
            y,x = np.histogram(dists, weights = weights_LD, bins=bins)
            LD_vs_distance[frag] = y/(1e-10+yn)
            y,x = np.histogram(dists, weights = weights_Dp, bins=bins)
            Dp_vs_distance[frag]=y/(1e-10+yn)

        for pcr in ['PCR1', 'PCR2']:
            positions, af2p, cov, af1p = control_LD(pcr, var_min=var_min)
            LD, Dp, p12 =  LDfunc(af2p, af1p, cov, cov_min=100)

            X,Y = np.meshgrid(positions, positions)
            np.fill_diagonal(cov, 0)
            dists = np.abs(X-Y)[cov>=cov_min].flatten()
            weights_LD = LD[cov>=cov_min].flatten()
            weights_Dp = Dp[cov>=cov_min].flatten()

            yn,xn = np.histogram(dists, bins = bins)
            y,x = np.histogram(dists, weights = weights_LD, bins=bins)
            LD_vs_distance[pcr] = y/(1e-10+yn)
            y,x = np.histogram(dists, weights = weights_Dp, bins=bins)
            Dp_vs_distance[pcr]=y/(1e-10+yn)
        data = {'Dp':Dp_vs_distance, 'LDrsq':LD_vs_distance, "bins":bins, "binc":binc,
                'var_min':var_min, 'cov_min':cov_min, 'dmin':dmin, 'dmin_pad':200, 'patients':patients}
        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)

    plot_LD(data, fig_filename = foldername+'LD')
