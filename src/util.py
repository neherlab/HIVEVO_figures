# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/05/15
content:    Generic utils for the paper figures.
'''
# Modules
import numpy as np


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

