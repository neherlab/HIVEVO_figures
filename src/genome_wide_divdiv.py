import numpy as np
from itertools import izip
from hivevo.patients import Patient
from hivevo.samples import all_fragments
from util import store_data, load_data, fig_width, fig_fontsize
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

    patients = ['p1', 'p2', 'p3','p5', 'p6', 'p8', 'p9', 'p11']
    if not os.path.isfile(fn_data) or params.redo:
        diversity = {}
        divergence = {}
        for pcode in patients:
            p = Patient.load(pcode)
            for frag in all_fragments:
                diversity[(pcode,frag)] = (p.ysi, p.get_diversity(frag))
                divergence[(pcode,frag)] = (p.ysi, p.get_divergence(frag))
    else:
        print("Loading data from file")
        data = load_data(fn_data)

    fig, axs = plt.subplots(2,3, sharey=True, sharex=True)
    for pcode in patients:
        for fi, frag in enumerate(all_fragments):
            ax = axs[fi//3][fi%3]
            ax.plot(diversity[(pcode,frag)][0], diversity[(pcode,frag)][1], 
                    '-o', label=pcode)
    axs[0][0].legend(loc=2, ncol=2)

    fig, axs = plt.subplots(2,3, sharey=True, sharex=True)
    for pcode in patients:
        for fi, frag in enumerate(all_fragments):
            ax = axs[fi//3][fi%3]
            ax.plot(divergence[(pcode, frag)][0], divergence[(pcode, frag)][1], 
                    '-o', label=pcode)
    axs[0][0].legend(loc=2, ncol=2)
