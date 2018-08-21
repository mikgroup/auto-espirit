#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from matplotlib import colors, ticker, cm
from scipy.io import loadmat
from math import ceil, floor, log10

##########################################################################################

dct = loadmat('../dem_code/res/map_dem.mat')

maps  = dct['maps']
fix_k = dct['fix_k'][0]
lst_w = dct['lst_w'][0]
lst_c = dct['lst_c'][0]

concat = np.concatenate

# maps are indexed as [lst_c, lst_w, imSize, nc, nc]

fig = plt.figure()
for cdx in range(len(lst_c)):
    for wdx in range(len(lst_w)):
        c = lst_c[cdx]
        w = lst_w[wdx]
        plt.subplot(len(lst_c), len(lst_w), len(lst_w) * cdx + wdx + 1)
        mp = np.squeeze(maps[cdx, wdx, :, :, 0:4, 7])
        im = concat((concat((mp[:, :, 0], mp[:, :, 1]), axis=1), concat((mp[:, :, 2], mp[:, :, 3]), axis=1)), axis = 0)
        plt.imshow(np.abs(im), cmap = 'gray')
        plt.xticks([])
        plt.yticks([])
        if cdx == 0:
            plt.title('w = %0.1f' % (w))
        if wdx == 0:
            plt.ylabel('c = %0.2f' % (c))

plt.tight_layout(pad=0.5, w_pad=-17.5, h_pad=0.5)
plt.savefig('pdf/variability.pdf'); plt.close(fig);
