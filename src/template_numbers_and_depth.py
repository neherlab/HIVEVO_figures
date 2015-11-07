# Modules
import os, sys
import glob
import numpy as np
from itertools import izip
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt

from hivevo.hivevo.patients import Patient
from hivevo.hivevo.samples import all_fragments

from util import store_data, load_data, draw_genome, fig_width, fig_fontsize, add_panel_label



# Script
if __name__=="__main__":

    patients = ['p'+str(i) for i in range(1,12) if i not in [4,7]]
    
    # make two figures, each showing one method of template quantification
    sns.set_style('darkgrid')
    fs=fig_fontsize
    fig1, ax1 = plt.subplots(figsize=(fig_width, 0.8*fig_width))
    fig2, ax2 = plt.subplots(figsize=(fig_width, 0.8*fig_width))
    add_panel_label(ax1, 'A', x_offset=-0.15)
    add_panel_label(ax2, 'B', x_offset=-0.15)

    # define colors for patients and fragments
    pat_colors = sns.color_palette(sns.color_palette(['#a6cee3', '#1f78b4',
                                                      '#b2df8a', '#33a02c',
                                                      '#fb9a99', '#e31a1c',
                                                      '#fdbf6f', '#ff7f00',
                                                      '#cab2d6'],
                                                     n_colors=len(patients)))
    frag_colors = sns.color_palette(n_colors=6)
    
    depth_estimates = []
    total_viral_load_dilutions_list = []
    overlap_dilution_list = {i:[] for i in range(6)}
    for pi, pcode in enumerate(patients):
        p = Patient.load(pcode)

        try:
            frag_depth = p.get_fragment_depth() #pad = True, limit_to_dilution = True)
            depth_estimates_pat = []
            for si, s in enumerate(p.samples):
                depth_estimates_pat.append(np.ma.concatenate([[p.n_templates_viral_load[si],
                                                               p.n_templates_dilutions[si]],
                                                              frag_depth[si]]))
            de_pat = np.array(depth_estimates_pat)
            # plot the virus load against template quantification by limiting dilution
            # one color and set of points per patient
            ax1.plot(de_pat[:,0], de_pat[:,1], 'o',
                     c=pat_colors[pi], alpha=0.5,
                     label=pcode)
            total_viral_load_dilutions_list.extend(zip(de_pat[:,0], de_pat[:,1]))

            # plot the fragment specific estimates, per fragment for time points of this patient
            for frag in range(6):
                good_vals= (de_pat[:,2+frag]>1)&(~np.isnan(de_pat[:,1]))
                ax2.plot(de_pat[good_vals,1], de_pat[good_vals,2+frag], 'o',
                         c=pat_colors[frag], alpha=0.5,
                         label='frag '+str(frag+1) if pcode=='p1' else None)
                overlap_dilution_list[frag].extend(zip(de_pat[:,1], de_pat[:,2+frag]))


            depth_estimates.extend(depth_estimates_pat)
        except (RuntimeError, IndexError, UnboundLocalError) as e:
            print '############\n',pcode,'\n############\n', 'ERROR', e


    # Calculate statistics on dilutions
    total_viral_load_dilutions_list =  np.array(total_viral_load_dilutions_list)
    dils = total_viral_load_dilutions_list[:, 1]
    print 'Estimated templates from dilutions:',
    print 'median:', np.median(dils),
    print 'quartiles:', np.percentile(dils, [25, 75])

    # calculate statistics of estimate concordance
    # mask nan and zero values
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

    ax1.plot([1,1e5], [1,1e5], lw=3, c='grey')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim(1,1e5)
    ax1.set_xlim(1,1e5)
    for item in ax1.get_xticklabels() + ax1.get_yticklabels():
        item.set_fontsize(fs)
    ax1.set_ylabel('Estimate from dilution', fontsize=fs)
    ax1.set_xlabel('Estimate from viral load', fontsize=fs)
    ax1.legend(loc=2, fontsize=fs*0.75, ncol=2)
    plt.tight_layout(rect=(0, 0, 0.98, 1))
    #for fmt in ['pdf', 'svg', 'png']:
    #    plt.savefig('figures/dilutions_vs_viral_load.'+fmt)


    ax2.plot([1,1e5], [1,1e5], lw=3, c='grey')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(1,1e5)
    ax2.set_xlim(1,1e5)
    for item in ax2.get_xticklabels() + ax2.get_yticklabels():
        item.set_fontsize(fs)
    ax2.set_xlabel('Estimate from dilutions', fontsize=fs)
    ax2.set_ylabel('Estimate from overlaps', fontsize=fs)
    ax2.legend(loc=2, fontsize=fs*0.75, ncol=2)
    plt.tight_layout(rect=(0, 0, 0.98, 1))
    #for fmt in ['pdf', 'svg', 'png']:
    #    plt.savefig('figures/fragments_vs_dilutions.'+fmt)
