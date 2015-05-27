import numpy as np
from itertools import izip
from hivevo.hivevo.patients import Patient
from hivevo.hivevo.HIVreference import HIVreference
from hivevo.hivevo.samples import all_fragments
from util import store_data, load_data, fig_width, fig_fontsize
import os
from filenames import get_figure_folder

def plot_subtype_correlation(data, fig_filename = None, figtypes=['.png', '.svg', '.pdf']):
    ####### plotting ###########
    import seaborn as sns
    from matplotlib import pyplot as plt
    plt.ion()
    sns.set_style('darkgrid')
    figpath = 'figures/'
    fs=fig_fontsize
    fig_size = (fig_width, 0.5*fig_width)
    fig, axs = plt.subplots(1, 2, sharey=True,figsize=fig_size)
    ax=axs[0]
    colors = {reg:c for reg, c in zip(data['spearmanr'].keys(), 
                                      sns.color_palette(n_colors=7))}

    for p,corrs in data['spearmanr'].iteritems():
        ax.plot(corrs['genomewide'][:,0], corrs['genomewide'][:,1], c=colors[p], label=p)
    if fig_filename is not None:
        for ext in figtypes:
            fig.savefig(fig_filename+'_sfs'+ext)
    else:
        plt.ion()
        plt.show()

if __name__=="__main__":
    import argparse
    from scipy.stats import spearmanr
    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action = 'store_true', help = 'recalculate data')
    params=parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'subtype_correlation.pickle'

    if not os.path.isfile(fn_data) or params.redo:
        patients = ['p2', 'p3','p5', 'p8', 'p9', 'p10','p11']
        time_bins = np.array([0, 200, 500, 1000, 1500, 2000, 3000, 5000])
        cov_min = 1000
        hxb2 = HIVreference(refname='HXB2')
        good_pos_in_reference = hxb2.get_ungapped(threshold = 0.05)
        subtype_correlation = {}
        for pi, pcode in enumerate(patients):
            subtype_correlation[pcode] = {}
            try:
                p = Patient.load(pcode)
            except:
                print "Can't load patient", pcode
            else:
                for region in ['genomewide', 'gag', 'pol', 'env']:
                    aft = p.get_allele_frequency_trajectories(region, cov_min=cov_min)
                    if len(aft.mask.shape)<2:
                        aft.mask = np.zeros_like(aft, dtype=bool)
                    corrs = []
                    patient_to_subtype = p.map_to_external_reference(region, refname = 'HXB2')
                    subtype_entropy = hxb2.get_entropy_in_patient_region(patient_to_subtype)
                    good_ref = good_pos_in_reference[patient_to_subtype[:,0]]
                    for t, af in izip(p.dsi,aft):
                        patient_entropy = np.maximum(0,-np.sum(af[:5]*np.log(1e-10+af[:5]), axis=0))[patient_to_subtype[:,2]]
                        good_af = (~np.any(af.mask, axis=0)[patient_to_subtype[:,2]]) & good_ref
                        if good_af.sum()>0.5*good_af.shape[0]:
                            corrs.append((t, spearmanr(subtype_entropy[good_af], patient_entropy[good_af])[0]))
                    subtype_correlation[pcode][region]=np.array(corrs)

        data = {'spearmanr':subtype_correlation}
        store_data(data, fn_data)
    else:
        print("Loading data from file")
        data = load_data(fn_data)

plot_subtype_correlation(data, fig_filename=None)
