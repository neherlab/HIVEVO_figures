import sys
sys.path.append('/ebio/ag-neher/share/users/rneher/HIVEVO_access')
import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.samples import all_fragments
from util import store_data, load_data, draw_genome, fig_width, fig_fontsize, add_panel_label
from scipy.stats import spearmanr

if __name__=="__main__":
    patients = ['p'+str(i) for i in range(1,12) if i not in [4,7]]
    import glob
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.ion()
    sns.set_style('darkgrid')
    fs=fig_fontsize
    fig1 = plt.figure(1, figsize=(fig_width, 0.8*fig_width))
    ax1 = plt.subplot(111)
    add_panel_label(ax1, 'A', x_offset=-0.15)
    fig2 = plt.figure(2, figsize=(fig_width, 0.8*fig_width))
    ax2 = plt.subplot(111)
    add_panel_label(ax2, 'B', x_offset=-0.15)
    pat_colors = sns.color_palette(n_colors=len(patients))
    frag_colors = sns.color_palette(n_colors=6)
    
    depth_estimates = []
    total_viral_load_dilutions_list = []
    overlap_dilution_list = {i:[] for i in range(6)}
    for pi, pcode in enumerate(patients):
        try:
            p = Patient.load(pcode)
            frag_depth = p.get_fragment_depth() #pad = True, limit_to_dilution = True)
            depth_estimates_pat = []
            for si, s in enumerate(p.samples):
                depth_estimates_pat.append(np.ma.concatenate([[p.n_templates_viral_load[si], p.n_templates_dilutions[si]], frag_depth[si]]))
            de_pat = np.array(depth_estimates_pat)
            plt.figure(1)
            plt.plot(de_pat[:,0], de_pat[:,1], 'o', c=pat_colors[pi], label=pcode)
            total_viral_load_dilutions_list.extend(zip(de_pat[:,0], de_pat[:,1]))
            plt.figure(2)
            for frag in range(6):
                good_vals= (de_pat[:,2+frag]>1)&(~np.isnan(de_pat[:,1]))
                plt.plot(de_pat[good_vals,1], de_pat[good_vals,2+frag], 'o', c=pat_colors[frag],
                         label='frag '+str(frag+1) if pcode=='p1' else None)
                overlap_dilution_list[frag].extend(zip(de_pat[:,1], de_pat[:,2+frag]))


            depth_estimates.extend(depth_estimates_pat)
        except (RuntimeError, IndexError, UnboundLocalError) as e:
            print '############\n',pcode,'\n############\n', 'ERROR',e

    total_viral_load_dilutions_list =  np.array(total_viral_load_dilutions_list)
    for f, val in overlap_dilution_list.iteritems():
        tmp=np.array(val)
        overlap_dilution_list[f]=tmp
    good_vals= (total_viral_load_dilutions_list[:,1]>1)&(~np.isnan(total_viral_load_dilutions_list[:,0]))
    print "viral load dilution rank correlation", spearmanr(total_viral_load_dilutions_list[good_vals,0], 
                                                            total_viral_load_dilutions_list[good_vals,1])
    for f, val in overlap_dilution_list.iteritems():
        good_vals= (val[:,1]>1)&(~np.isnan(val[:,0]))
        #print val[good_vals,:]
        print "dilution overlap, fragment",str(f+1), spearmanr(val[:,0], val[:,1])
    tmp = np.vstack(overlap_dilution_list.values())
    print "dilution overlap, all fragments", spearmanr(tmp[:,0], tmp[:,1])

    plt.figure(1)
    plt.plot([1,1e5], [1,1e5], lw=3, c='grey')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1,1e5)
    plt.xlim(1,1e5)
    for item in ax1.get_xticklabels() + ax1.get_yticklabels():
        item.set_fontsize(fs)
    plt.ylabel('Estimate from dilution', fontsize=fs)
    plt.xlabel('Estimate from viral load', fontsize=fs)
    plt.legend(loc=2, fontsize=fs*0.75, ncol=2)
    plt.tight_layout(rect=(0, 0, 0.98, 1))
    for fmt in ['pdf', 'svg', 'png']:
        plt.savefig('figures/dilutions_vs_viral_load.'+fmt)


    plt.figure(2)
    plt.plot([1,1e5], [1,1e5], lw=3, c='grey')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1,1e5)
    plt.xlim(1,1e5)
    for item in ax2.get_xticklabels() + ax2.get_yticklabels():
        item.set_fontsize(fs)
    plt.xlabel('Estimate from dilutions', fontsize=fs)
    plt.ylabel('Estimate from overlaps', fontsize=fs)
    plt.legend(loc=2, fontsize=fs*0.75, ncol=2)
    plt.tight_layout(rect=(0, 0, 0.98, 1))
    for fmt in ['pdf', 'svg', 'png']:
        plt.savefig('figures/fragments_vs_dilutions.'+fmt)
