# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/05/15
content:    Generic utils for the paper figures.
'''
# Modules
import numpy as np

fig_width = 5  
fig_fontsize = 12  

# Functions
def HIVEVO_colormap(kind='website'):
    from scipy.interpolate import interp1d
    maps = {'website': ["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52",
                        "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"],
            'alternative': ["#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                            "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"],
           }
    colors = maps[kind]
    rgb_colors = []
    for c in colors:
        rgb_colors.append([int(c[i:i+2],16) for i in [1,3,5]]+[255])
    tmp =interp1d(np.linspace(0,1,len(colors)), np.array(rgb_colors, dtype = float).T/255.0)
    cmap = lambda x: [c for c in tmp(x)]
    return cmap

def add_binned_column(df, bins, to_bin):
    df.loc[:,to_bin+'_bin'] = np.minimum(len(bins)-2, np.maximum(0,np.searchsorted(bins, df.loc[:,to_bin])-1))

def store_data(data, fn):
    '''Store data to file for the plots'''
    import cPickle as pickle
    with open(fn, 'wb') as f:
        pickle.dump(data, f, protocol=-1)


def load_data(fn):
    '''Load the data for the plots'''
    import cPickle as pickle
    with open(fn, 'rb') as f:
        return pickle.load(f)

def draw_genome(ax, annotations,rows=3, readingframe=True,fs=9):
    from matplotlib.patches import Rectangle
    y1 ,height, pad = 0, 1, 0.2
    ax.set_ylim([-pad,rows*(height+pad)])
    anno_elements = []
    for name, feature in annotations.iteritems():
        x = [feature.location.nofuzzy_start, feature.location.nofuzzy_end]
        anno_elements.append({'name': name,
                     'x1': x[0], 'x2': x[1], 'width': x[1] - x[0]})
    anno_elements.sort(key = lambda x:x['x1'])
    for ai, anno in enumerate(anno_elements):
        if readingframe:
            anno['y1'] = y1 + (height + pad) * (2 - (anno['x1'])%3)
        else:
            anno['y1'] = y1 + (height + pad) * (ai%rows)
        anno['y2'] = anno['y1'] + height

    for anno in anno_elements:
        r = Rectangle((anno['x1'], anno['y1']),
                      anno['width'],
                      height,
                      facecolor=[0.8] * 3,
                      edgecolor='k',
                      label=anno['name'])
        
        xt = anno['x1'] + 0.5 * anno['width']
        yt = anno['y1'] + 0.2 * height + height*(anno['width']<500)

        ax.add_patch(r)
        ax.text(xt, yt,
                anno['name'],
                color='k', 
                fontsize=fs,
                ha='center')