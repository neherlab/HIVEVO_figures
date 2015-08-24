# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Plot the quality along the read pair.
'''
# Modules
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from filenames import get_figure_folder
from util import HIVEVO_colormap, tree_from_json
from util import fig_width, fig_fontsize, draw_tree


# Globals
regions = ['p17', 'V3']
cutoff = 0.04

web_colormap = HIVEVO_colormap()



# Functions
def collect_data(regions):
    '''Get the trees in JSON'''
    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'

    data = []
    for region in regions:
        fn = fn_data+'haplotype_tree_p1_'+region+'.json'
        data.append({'region': region,
                     'tree': tree_from_json(fn)})

    return data


def plot_haplotype_trees(data, fig_filename=None):
    '''Plot tree of minor haplotypes in a typical patient'''
    import seaborn as sns
    from matplotlib import pyplot as plt

    def assign_color(tree, cmap='jet', attrname='DSI'):
        '''Assign color to leaves based on a colormap'''
        if isinstance(cmap, basestring):
            from matplotlib import cm
            cmap = getattr(cm, cmap)

        attr_max = max(getattr(leaf, attrname) for leaf in tree.get_terminals())
        def get_color(node):
            return map(int, np.array(cmap(getattr(node, attrname)/attr_max)[:-1]) * 255)
    
        for node in tree.get_terminals():
            node.color = get_color(node)

        # For internal nodes, set the attribute (e.g. age) as the arithmetic mean of
        # the children clades
        for node in tree.get_nonterminals(order='postorder'):
            setattr(node, attrname, np.mean([getattr(c, attrname) for c in node.clades]))
            node.color = get_color(node)


    plt.ioff()

    if VERBOSE:
        print 'Plot haplotype tree of typical patient'

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    sns.set_style('white')
    ax.grid(False)

    x_offset = 0
    y_offset = 15
    y_padding = 15

    for datum in data:
        tree = datum['tree']
        tree.root.branch_length = 0.01

        times = sorted(set([leaf.DSI for leaf in tree.get_terminals()]))

        assign_color(tree, attrname='DSI', cmap = web_colormap)
        labels = ['{:>8s}'.format('{:>3.0%}'.format(leaf.frequency)+' '+
                   '{:>2d}'.format(int(float(leaf.DSI) / 30.5))+
                  'm')
                  for leaf in tree.get_terminals()]
        depths = tree.depths()
        maxdepth = max(depths.itervalues())
        mindepth = min(depths.itervalues())

        # Collect data for circle plot
        rmin = 5
        rmax = 150
        rfun = lambda hf: rmin + (rmax - rmin) * (hf - 0.5 * mindepth)**(0.5)
        data_circles = []
        for il, leaf in enumerate(tree.get_terminals(), 1):
            it = times.index(leaf.DSI)
            hf = leaf.frequency
            r = rfun(hf)
            y = il + y_offset
            x = depths[leaf] + x_offset
            c = [tone / 255.0 for tone in leaf.color.to_rgb()]
            cs = [tone / 255.0 * 0.7 for tone in leaf.color.to_rgb()]
            data_circles.append((x, y, 2 * r, c, cs))

        # Draw the tree
        draw_tree(tree, show_confidence=False,
                   label_func=lambda x: '',
                   axes=ax,
                  x_offset=x_offset,
                  y_offset=y_offset,
                   do_show=False)

        # Add circles to the leaves
        (x, y, s, c,cs) = zip(*data_circles)
        ax.scatter(x, y, s=s, c=c, edgecolor=cs, zorder=2)
        ax.set_xlim(-0.04 * maxdepth, 1.04 * maxdepth)

        # Add region name
        ax.text(tree.root.branch_length // 2, y_offset + 10,
                datum['region'], fontsize=16,
                va='top')

        y_offset += tree.count_terminals() + y_padding

    ax.set_ylim((y_offset + y_padding, 0))
    ax.set_ylabel('')
    ax.set_yticklabels([])
    ax.set_axis_off()
    ax.xaxis.set_tick_params(labelsize=16)
    ax.set_xlabel('Genetic distance [changes / site]',
                  fontsize=16,
                  labelpad=10)


    # Draw a "legend" for sizes
    datal = [{'hf': 0.05, 'label': '5%'},
             {'hf': 0.20, 'label': '20%'},
             {'hf': 1.00, 'label': '100%'}]
    ax.text(0.98 * maxdepth, 0.03 * ax.get_ylim()[0],
            'Haplotype frequency:', fontsize=16, ha='right')
    for idl, datuml in enumerate(datal):
        r = rfun(datuml['hf'])
        y = (0.07 + 0.07 * idl) * ax.get_ylim()[0]
        ax.scatter(0.85 * maxdepth, y, s=r,
                   facecolor='k',
                   edgecolor='none')
        ax.text(0.98 * maxdepth, y + 0.02 * ax.get_ylim()[0],
                datuml['label'], ha='right',
                fontsize=14)

    # Draw legend for times
    tone_fun = lambda t: np.array(web_colormap(1.0 * t / max(times))[:3])
    datal = [{'time': t,
              'color': tone_fun(t).tolist(),
              'colorstroke': (0.7 * tone_fun(t)).tolist(),
             }
             for t in [times[0], times[len(times) // 2], times[-1]]]
    xtext = 0
    ytext = (0.97 - 0.06 * min(8, len(datal))) * ax.get_ylim()[0]
    ax.text(0.01 * maxdepth, ytext, 'Time:', fontsize=16)
    for ico, datuml in enumerate(datal):
        ytext += 0.06 * (ax.get_ylim()[0] - ax.get_ylim()[1])

        if ico == 8:
            ytext -= 4 * 0.06 * (ax.get_ylim()[0] - ax.get_ylim()[1])
            xtext += 0.28 * maxdepth

        ax.scatter(xtext, ytext, s=rfun(0.5),
                   facecolor=datuml['color'],
                   edgecolor=datuml['colorstroke'])
        t_text = int(datuml['time'] / 30.5)
        if t_text == 1:
            t_text = str(t_text)+ 'month'
        else:
            t_text = str(t_text)+ 'months'
        ax.text(xtext + 0.21 * maxdepth,
                ytext + 0.02 * ax.get_ylim()[0],
                t_text,
                ha='right',
                fontsize=14,
               )

    # Draw scale bar
    xbar = (0.3 + 0.3 * (len(datal) >= 9)) * maxdepth
    ybar = 0.90 * ax.get_ylim()[0]
    lbar = 0.05 * maxdepth
    lbar_label = '{:.1G}'.format(lbar)
    lbar = float(lbar_label)
    ax.plot([xbar, xbar + lbar], [ybar, ybar], lw=4, c='k')
    ax.text(xbar + 0.5 * lbar, ybar + 0.08 * ax.get_ylim()[0],
            lbar_label, fontsize=14,
            ha='center')


    plt.tight_layout(rect=(0, 0, 0.98, 1))

    if fig_filename is not None:
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout(rect=(0, 0, 0.98, 1))
        plt.ion()
        plt.show()




# Script
if __name__ == '__main__':

    VERBOSE = 2
    data = collect_data(regions)
            
    #filename = get_tree_figure_filename(patient.code,
    #                                    region,
    #                                    format='svg')

    plot_haplotype_trees(data,
                         fig_filename=None)
