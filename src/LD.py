import numpy as np
import sys
from itertools import izip
from hivevo.patients import Patient
from hivevo.samples import all_fragments
from hivevo.af_tools import LD as LDfunc
from hivwholeseq.filenames import root_data_folder


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


patients = ['p' +str(i) for i in xrange(1,12) if i!=4]


if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description="produce linkage disequilibrium figure")
    parser.add_argument('--fragment', default = 'F1', nargs='+',
                        help = 'fragmet')
    parser.add_argument('--plot', action = 'store_true', help = 'plot LD figures')
    parser.add_argument('--controls', action = 'store_true', help = 'include control experiment in plot')
    params=parser.parse_args()

    dmin = 40
    dmin_pad = 200
    var_min = 0.2
    cov_min = 200
    LD_vs_distance = {}
    Dp_vs_distance = {}
    bins = np.arange(0,401,40)
    binc = (bins[:-1]+bins[1:])*0.5
    for frag in params.fragment:
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

    if params.plot:
        import seaborn as sns
        from matplotlib import pyplot as plt
        plt.ion()
        sns.set_style('darkgrid')
        figpath = 'figures/'
        fs=16
        fig_size = (5.5, 4.3)

        ########### figure with r^2 #################
        fig = plt.figure(1, figsize=fig_size)
        ax = plt.subplot(111)
        for frag in params.fragment:
            y = LD_vs_distance[frag]
            plt.plot(binc, y, label=frag, lw=3)
        if params.controls:
            plt.plot(binc, LD_vs_distance['PCR1'], label="control", lw=2, ls='--', c='k')

        ax.set_xticks(range(0,401,100))
        ax.set_yticks(np.arange(0.,1.01,0.2))
        for item in ax.get_xticklabels() + ax.get_yticklabels():
            item.set_fontsize(fs)
        plt.ylabel('linkage disequilibrium r^2',fontsize = fs)
        plt.xlabel("distance [bp]",fontsize = fs)
        plt.legend(loc='center right', ncol=2)
        plt.tight_layout(rect=(0, 0, 0.98, 1))
        for ext in ['pdf','svg']:
            plt.savefig(figpath+'LD_rsq.'+ext)

        ########### figure with Dp #################
        fig = plt.figure(2, figsize=fig_size)
        ax = plt.subplot(111)
        for frag in params.fragment:
            y = Dp_vs_distance[frag]
            plt.plot(binc, y, label=frag, lw=3)

        if params.controls:
            plt.plot(binc, Dp_vs_distance['PCR1'], label='control', lw=2, ls='--', c='k')

        ax.set_xticks(range(0,401,100))
        ax.set_yticks(np.arange(0.,1.01,0.2))
        for item in ax.get_xticklabels() + ax.get_yticklabels():
            item.set_fontsize(fs)
        plt.ylabel("linkage disequilibrium D'",fontsize = fs)
        plt.xlabel("distance [bp]",fontsize = fs)
        plt.legend(loc='center right', ncol=2)
        plt.tight_layout(rect=(0, 0, 0.98, 1))
        for ext in ['pdf','svg', 'png']:
            plt.savefig(figpath+'LD_Dp.'+ext)


        ########### figure with Controls only #################
        if params.controls:
            fig = plt.figure(3, figsize=fig_size)
            ax = plt.subplot(111)
    
            plt.plot(binc, Dp_vs_distance['PCR1'], label='PCR1', lw=2, ls='-')
            plt.plot(binc, Dp_vs_distance['PCR2'], label='PCR2', lw=2, ls='-')
    
            ax.set_xticks(range(0,401,100))
            ax.set_yticks(np.arange(0.,1.01,0.2))
            for item in ax.get_xticklabels() + ax.get_yticklabels():
                item.set_fontsize(fs)
            plt.ylabel("linkage disequilibrium D'",fontsize = fs)
            plt.xlabel("distance [bp]",fontsize = fs)
            plt.legend(loc=3)
            plt.tight_layout(rect=(0, 0, 0.98, 1))
            for ext in ['pdf','svg', 'png']:
                plt.savefig(figpath+'LD_Dp_control.'+ext)
    