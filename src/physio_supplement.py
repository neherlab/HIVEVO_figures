# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/11/15
content:    Plot viral load and CD4+ counts as a supplementary, 9-panel figure.
'''
# Modules
import os, sys
import numpy as np
import pandas as pd
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import seaborn as sns

from hivevo.hivevo.patients import Patient

from util import fig_width, fig_fontsize


# Globals
pnumbers = [1, 2, 3, 5, 6, 8, 9, 10, 11]
colors = {'CD4': 'darkred',
          'VL': 'steelblue'}


# Functions
def get_table_filename():
    return 'data/More CD4 and VL 151105.xlsx'


def get_physio(sheet, pnumbers=pnumbers):
    fn = get_table_filename()
    df = (pd.read_excel(fn, sheet, parse_cols=[0, 3, 4])
            .rename(columns={'Pnr': 'patient number',
                             'DecVal': sheet}))
    df_list = {'p'+str(pn): d[['ETI', sheet]]
               for pn, d in df.groupby('patient number')
               if pn in pnumbers}
    print df_list.keys()
    return df_list


def get_CD4():
    return get_physio('CD4')


def get_VL():
    return get_physio('VL')


def plot_physio(data):
    pass


def collect_data():
    data = {'CD4': get_CD4(),
            'VL': get_VL(),
            'deep sequencing': {}}

    for pn in pnumbers:
        pcode = 'p'+str(pn)
        p = Patient.load(pcode)
        data['deep sequencing'][pcode] = p.dsi

    return data



# Script
if __name__ == '__main__':
    data = collect_data()

    fig_size = (2 * fig_width, 2 * fig_width)
    fig, axs = plt.subplots(3, 3, figsize=fig_size)
    axs = axs.ravel()
    sns.set_style('dark')

    for i, ax in enumerate(axs):
        pn = pnumbers[i]
        pcode = 'p'+str(pn)
        
        for iax, dtype in enumerate(['VL', 'CD4']):
            datum = data[dtype][pcode]

            x = np.array(datum['ETI'], int)
            y = np.array(datum[dtype], float)
            
            ax.plot(x, y, '-o',
                    lw=2,
                    zorder=100*(iax+1),
                    color=colors[dtype])

            ax.set_xlim(0, 3100)
            
            if (iax == 0) and (i == 3):
                ax.set_ylabel('Viral load [counts/ml]')
                fig.patches.append(Circle((0.025, 0.415), radius=0.008,
                                       edgecolor='none',
                                       facecolor=colors[dtype],
                                       transform=fig.transFigure))
                

            if (iax == 1) and (i == 5):
                ax.set_ylabel('CD4+ counts [cells/ul]', rotation=270,
                              labelpad=17)
                fig.patches.append(Circle((0.974, 0.61), radius=0.008,
                                       edgecolor='none',
                                       facecolor=colors[dtype],
                                       transform=fig.transFigure))

            if (iax == 0) and (i == 7):
                ax.set_xlabel('ETI [days]')

            if iax == 0:
                ax.set_title(pcode)
                ax.set_ylim(1e1, 1e6)
                ax.set_yscale('log')

                if i not in [0, 3, 6]:
                    ax.set_yticklabels([])


                if i not in [6, 7, 8]:
                    ax.set_xticklabels([])

                # Add sequencing arrows
                ds = data['deep sequencing'][pcode]
                for dstime in ds:
                    ax.annotate('',
                                xy=(dstime, 3e1),
                                xytext=(dstime, 1e1),
                                arrowprops=dict(arrowstyle='->'))

                ax.grid(True, axis='both')

                ax = ax.twinx()
            else:
                ax.set_ylim(0, 1500)
                ax.set_yticks([0, 300, 600, 900, 1200, 1500])

                if i not in [2, 5, 8]:
                    ax.set_yticklabels([])

    plt.tight_layout()

    fig.savefig('figures/physiological.svg')
    fig.savefig('figures/physiological.pdf')
    fig.savefig('figures/physiological.png')

    plt.ion()
    plt.show()
