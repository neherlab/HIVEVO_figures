# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/05/15
content:    Make figure for the substitutions and CTL epitopes.
'''
# Modules
import os
import argparse
from itertools import izip
import numpy as np
import pandas as pd

from hivevo.hivevo.sequence import alpha
from hivevo.hivevo.patients import Patient
from util import store_data, load_data, draw_genome, fig_width, fig_fontsize
from filenames import get_figure_folder



# Functions
def collect_ctl_data(patients, regions, ctl_kind='mhci=80'):
    data_ctl = []

    for pi, pcode in enumerate(patients):
        p = Patient.load(pcode)

        # Add predicted epitopes
        ctl_table = p.get_ctl_epitopes(kind=ctl_kind, regions=regions)
        ctl_table['pcode'] = p.name
        data_ctl.append(ctl_table)
    data_ctl = pd.concat(data_ctl)
    return data_ctl


def collect_substitution_data(patients, regions, cov_min=100):
    from Bio.Seq import translate
    data = []
    for pi, pcode in enumerate(patients):
        p = Patient.load(pcode)

        for region in regions:
            print p.name, region

            initial_indices = p.get_initial_indices(region)
            aft = p.get_allele_frequency_trajectories(region, cov_min=cov_min)
            if np.isscalar(aft.mask):
                aft.mask = np.zeros_like(aft, bool)

            coomap = p.map_to_external_reference(region)[:, ::2]
            coomapd = {'pat_to_subtype': dict(coomap[:, ::-1]),
                       'subtype_to_pat': dict(coomap)}

            for posdna in xrange(aft.shape[-1]):
                # Get the position in reference coordinates
                if posdna not in coomapd['pat_to_subtype']:
                    pos_sub=-1 #continue
                else:
                    pos_sub = coomapd['pat_to_subtype'][posdna]

                # Get allele frequency trajectory
                aftpos = aft[:, :, posdna]
                ind = -aftpos[:, 0].mask
                if ind.sum() == 0:
                    continue
                aftpos = aftpos[ind]
                timespos = p.dsi[ind]

                # Ancestral allele
                ianc = initial_indices[posdna]
                anc = alpha[ianc]

                # Ignore indels
                if ianc >= 4:
                    continue

                # Check for fixation
                if (aftpos[0, ianc] < 0.7) or np.min(aftpos[:, ianc]) > 0.2:
                    continue

                # Get codon
                ci = posdna // 3
                rf = posdna % 3
                cod_anc = ''.join(alpha[initial_indices[ci * 3: (ci + 1) * 3]])
                if '-' in cod_anc:
                    continue
                aa_anc = translate(cod_anc)

                # Check which allele (if any) is fixing
                for inuc, nuc in enumerate(alpha[:4]):
                    if nuc == anc:
                        continue
                    
                    if aftpos[-1, inuc] < 0.95:
                        continue

                    # NOTE: OK, it's a substitution (max 1 per site)
                    break
                else:
                    continue

                # Assign a time to the substitution
                ist = (aftpos[:, inuc] > 0.5).nonzero()[0][0]
                tsubst = 0.5 * (timespos[ist - 1] + timespos[ist])

                nuc = alpha[inuc]
                mut = anc+'->'+nuc

                # Define transition/transversion
                if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                    trclass = 'ts'
                else:
                    trclass = 'tv'

                # Check syn/nonsyn
                cod_nuc = cod_anc[:rf] + nuc + cod_anc[rf+1:]
                aa_nuc = translate(cod_nuc)
                is_syn = aa_nuc == aa_anc

                datum = {'pcode': p.name,
                         'region': region,
                         'pos_patient': posdna,
                         'pos_ref': pos_sub,
                         'mut': mut,
                         'trclass': trclass,
                         'syn': is_syn,
                         'time': tsubst,
                        }

                data.append(datum)


    data = pd.DataFrame(data)
    return data


def correlate_epitope_substitution(ds, dctl):
    '''Correlate presence of a substitution with epitope'''
    from hivwholeseq.sequencing.primer_info import primers_coordinates_HXB2_outer
    start_F1 = primers_coordinates_HXB2_outer['F1'][0][1]
    end_F6 = primers_coordinates_HXB2_outer['F6'][1][0]

    ds = ds.copy()

    dg = []
    for pcode, datum in dctl.groupby('pcode'):
        a = np.arange(start_F1, end_F6)
        b = np.zeros(len(a), bool)
        for _, epi in datum.iterrows():
            b[(a >= epi['start_HXB2']) & (a < epi['end_HXB2'])] = True
        c = np.zeros(len(a), bool)
        datum = ds.loc[ds['pcode'] == pcode]
        # Keep only nonsyn substitutions
        datum = datum.loc[datum['syn'] == False]
        c[datum['pos_ref'] - a[0]] = True
        dat = {'pos': a,
               'epitope': b,
               'substitution': c,
               }
        dat = pd.DataFrame(dat)
        dat['pcode'] = pcode
        dg.append(dat)
    dg = pd.concat(dg)

    # Exclude env because it has antibody-related substitutions
    from hivwholeseq.reference import load_custom_reference
    from hivwholeseq.utils.sequence import find_annotation
    ref = load_custom_reference('HXB2', 'gb')
    start_env = find_annotation(ref, 'gp41').location.nofuzzy_start 
    end_env = find_annotation(ref, 'gp41').location.nofuzzy_end - 450
    dg = dg.loc[(dg['pos'] < start_env) | (dg['pos'] >= end_env)]

    M = dg.groupby(['epitope', 'substitution']).size().unstack()
    Ma = np.array(M)
    xp = 1.0*Ma[1,0]/Ma[0,0]*Ma[0,1]
    xs = Ma[1,1] - xp
    print M
    from scipy.stats import fisher_exact
    print 'Fisher\'s exact enrichment:', fisher_exact(Ma)[0]
    print 'Fisher\'s exact P value:', fisher_exact(Ma)[1]
    print 'expected:', xp
    print 'excess:',  xs, 'per patient:',xs/9.0

    pos_epi = dg.loc[dg['epitope'] == True]['pos'].unique()
    dg2 = dg.loc[dg['pos'].isin(pos_epi)].copy()
    M2 = dg2.groupby(['epitope', 'substitution']).size().unstack()
    M2a = np.array(M2)
    xp = 1.0*M2a[1,0]/M2a[0,0]*M2a[0,1]
    xs = M2a[1,1] - xp
    print M2
    print '\nFisher\'s exact enrichment:', fisher_exact(M2a)[0]
    print 'Fisher\'s exact P value:', fisher_exact(M2a)[1]
    print 'expected:', xp
    print 'excess:',  xs, 'per patient:',xs/9.0
    return {'dg': dg,
            'dg2': dg2,
           }


def number_substitutions_per_patient(data):
    '''Calculate the median and quartiles of substitutions per patient'''
    print data.groupby('pcode').size()


def plot_ctl_epitopes(data, ax=None, yoffset=0, colormap=None, fs=None):
    '''Plot the predicted CTL epitopes'''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import seaborn as sns

    sns.set_style('darkgrid')

    if fs is None:
        from util import fig_fontsize
        fs = fig_fontsize

    if colormap is None:
        colormap = cm.jet
    elif isinstance(colormap, basestring):
        from util import HIVEVO_colormap
        colormap = HIVEVO_colormap(kind=colormap)

    pnames = sorted(data['pcode'].unique(), key=lambda x: int(x[1:]))
    if ax is None:
        from util import fig_width
        fig, ax = plt.subplots(figsize=(2 * fig_width, 1 + 0.2 * len(pnames)))
    for pcode, datum in (data
                         .loc[:, ['pcode', 'start_HXB2', 'end_HXB2']]
                         .groupby('pcode')):

        y = pnames.index(pcode) + yoffset
        for _, d in datum.iterrows():
            ax.plot([d['start_HXB2'], d['end_HXB2']], [y] * 2, lw=2,
                   color=colormap(1.0 * y / len(pnames))
                   )

    ax.set_xlabel('Position in reference [bp]', fontsize=fs)
    ax.set_yticks(range(len(pnames)))
    ax.set_yticklabels(pnames, fontsize=fs)
    ax.set_ylim(len(pnames) - 0.5 + yoffset, -0.5)
    ax.grid(axis='y')
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)

    plt.tight_layout()



def plot_substitutions(data, fig_filename=None, figtypes=['.png', '.svg', '.pdf']):
    '''Plot substitution locations and CTL epitopes'''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import seaborn as sns

    sns.set_style('darkgrid')

    from util import HIVEVO_colormap, fig_width, fig_fontsize
    fs = fig_fontsize
    colormap = HIVEVO_colormap(kind='website')

    pnames = sorted(data['ctl']['pcode'].unique(), key=lambda x: int(x[1:]))
    fig, ax = plt.subplots(figsize=(2 * fig_width, 1 + 0.2 * len(pnames)))

    for pcode, datum in (data['substitutions']
                      .loc[:, ['pcode', 'pos_ref']]
                      .groupby('pcode')):

        x = np.array(datum['pos_ref'])
        y = np.repeat(pnames.index(pcode), len(x))

        ax.scatter(x, y,
                   s=20,
                   color=colormap(1.0 * pnames.index(pcode) / len(pnames)),
                  ) 

    plot_ctl_epitopes(data['ctl'], ax=ax, yoffset=0.4, colormap=colormap, fs=fs)
    if fig_filename is not None:
        for ext in figtypes:
            fig.savefig(fig_filename+ext)
        plt.close(fig)
    else:
        plt.ion()
        plt.show()





if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Make figure for substitutions and CTL epitopes")
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    params = parser.parse_args()

    VERBOSE = 2
    window_size = 300
    ctl_kind = 'mhci=80'

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    fn_data = fn_data + 'substitutions_CTL.pickle'

    if not os.path.isfile(fn_data) or params.redo:
        patients = ['p1', 'p2', 'p3', 'p5', 'p6', 'p8', 'p9', 'p10', 'p11']
        # FIXME: add more regions
        regions = ['gag', 'pol', 'gp120_noVloops', 'gp41', 'vif', 'vpu', 'vpr', 'nef']

        ds = collect_substitution_data(patients, regions)

        dctl = collect_ctl_data(patients, regions, ctl_kind=ctl_kind)
        correlate_epitope_substitution(ds, dctl)

        data = {'substitutions': ds,
                'ctl': dctl,
               }

        store_data(data, fn_data)
    else:
        data = load_data(fn_data)


    correlate_epitope_substitution(data['substitutions'], data['ctl'])

    number_substitutions_per_patient(data['substitutions'])

    #plot_ctl_epitopes(data['ctl'])
    #plot_substitutions(data, fig_filename=foldername+'substitutions')
