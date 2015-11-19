import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.samples import all_fragments
from util import store_data, load_data, fig_width, fig_fontsize, patients, patient_colors
import os
from filenames import get_figure_folder
import matplotlib.pyplot as plt
import seaborn as sns
plt.ion()

sns.set_style('darkgrid')


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action = 'store_true', help = 'recalculate data')
    params=parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'genomewide_divdiv.pickle'

    csv_out = open(foldername+'/genomewide_divdiv.tsv','w')
    if not os.path.isfile(fn_data) or params.redo:
        #####
        ## recalculate diversity and divergence from allele frequencies
        #####
        diversity = {}
        divergence = {}
        for pcode in patients:
            p = Patient.load(pcode)
            for frag in all_fragments:
                diversity[(pcode,frag)] = (p.ysi, p.get_diversity(frag))
                divergence[(pcode,frag)] = (p.ysi, p.get_divergence(frag))
        store_data((diversity, divergence), fn_data)
    else:
        print("Loading data from file")
        diversity, divergence = load_data(fn_data)

#####
## plot diversity 
#####
    fig, axs = plt.subplots(2,3, sharey=True, sharex=True)
    for pcode in patients:
        for fi, frag in enumerate(all_fragments):
            ax = axs[fi//3][fi%3]
            ax.plot(diversity[(pcode,frag)][0], diversity[(pcode,frag)][1], 
                    '-o', label=pcode, c=patient_colors[pcode])
            csv_out.write('\t'.join(map(str, ['diversity_'+pcode+'_'+frag]+list(diversity[(pcode,frag)][1])))+'\n')
            csv_out.write('\t'.join(map(str, ['time_'+pcode+'_'+frag]+list(diversity[(pcode,frag)][0])))+'\n')
    for ax in axs[:,0]:
        ax.set_ylabel('diversity')
        ax.locator_params(nbins=5)
        ax.tick_params(axis='both', labelsize = fig_fontsize)
    for ax in axs[-1,:]:
        ax.set_xlabel('Years since EDI]')
        ax.locator_params(nbins=5)
        ax.tick_params(axis='both', labelsize = fig_fontsize)
    axs[0][0].set_ylim([0,0.017])
    axs[0][0].set_xlim([-0.1,4])
    for fi, frag in enumerate(all_fragments):
        ax = axs[fi//3][fi%3]
        ax.text(0.02,0.9, frag, fontsize=1.5*fig_fontsize, transform=ax.transAxes)
    axs[0][0].legend(loc=1, ncol=2)
    plt.tight_layout()
    for fmt in ['.pdf', '.svg', '.png']:
        plt.savefig(foldername+'genomewide_diversity'+fmt)

#####
## plot divergence
#####
    fig, axs = plt.subplots(2,3, sharey=True, sharex=True)
    for pcode in patients:
        for fi, frag in enumerate(all_fragments):
            ax = axs[fi//3][fi%3]
            ax.plot(divergence[(pcode, frag)][0], divergence[(pcode, frag)][1], 
                    '-o', label=pcode, c=patient_colors[pcode])
            csv_out.write('\t'.join(map(str, ['divergence_'+pcode+'_'+frag]+list(divergence[(pcode,frag)][1])))+'\n')
            csv_out.write('\t'.join(map(str, ['time_'+pcode+'_'+frag]+list(divergence[(pcode,frag)][0])))+'\n')
    for ax in axs[:,0]:
        ax.set_ylabel('diversity')
        ax.locator_params(nbins=5)
        ax.tick_params(axis='both', labelsize = fig_fontsize)
    for ax in axs[-1,:]:
        ax.set_xlabel('Years since EDI]')
        ax.locator_params(nbins=5)
        ax.tick_params(axis='both', labelsize = fig_fontsize)
    axs[0][0].set_ylim([0,0.05])
    for fi, frag in enumerate(all_fragments):
        ax = axs[fi//3][fi%3]
        ax.text(0.8,0.02, frag, fontsize=1.5*fig_fontsize, transform=ax.transAxes)
    axs[0][0].legend(loc=2, ncol=2)
    plt.tight_layout()
    for fmt in ['.pdf', '.svg', '.png']:
        plt.savefig(foldername+'genomewide_divergence'+fmt)
    csv_out.close()
