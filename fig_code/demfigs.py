#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import matplotlib

from matplotlib import colors, ticker, cm
from scipy.io import loadmat
from math import ceil, floor, log10
from cfl import readcfl

##########################################################################################

plt.rc('font',   size=10)       # controls default text sizes
plt.rc('axes',   titlesize=10)  # fontsize of the axes title
plt.rc('axes',   labelsize=10)  # fontsize of the x and y labels
plt.rc('xtick',  labelsize=10)  # fontsize of the tick labels
plt.rc('ytick',  labelsize=10)  # fontsize of the tick labels
plt.rc('legend', fontsize=10)   # legend fontsize
plt.rc('figure', titlesize=10)  # fontsize of the figure title

concat = np.concatenate
cmap = matplotlib.cm.jet
cmap.set_bad(color='pink')
cmap = matplotlib.cm.gray
cmap.set_bad(color='pink')
##########################################################################################

f = lambda x: np.divide(1, x, out=np.zeros_like(x), where=(x!=0))
g = lambda x: np.min(x[np.abs(x) > 0])
h = lambda x: f(np.mean(x[np.abs(x) > 0]))

dct  = loadmat('../dem_code/res/map_dem.mat')
otl  = np.abs(readcfl('../data/brain_outline'))
msk  = np.abs(readcfl('../data/brain_fov'))
mskt = concat([concat([msk, msk], axis=0), concat([msk, msk], axis=0)], axis=1);

def ic(x):
  img = x.copy()
  img[(abs(msk) > 0) & (abs(x) == 0)] = np.inf
  return img

maps      = dct['maps']
gfactor1  = dct['gfactor_rx_1_ry_2']
gfactor2  = dct['gfactor_rx_2_ry_2']
igfactor1 = f(gfactor1)
igfactor2 = f(gfactor2)
fix_k     = dct['fix_k'][0]
lst_w     = dct['lst_w'][0]
lst_c     = dct['lst_c'][0][::2]

fig, axes = plt.subplots(nrows=len(lst_c), ncols=len(lst_w))
ctr = 0;
for ax in axes.flat:
  wdx = ctr % len(lst_w)
  cdx = ctr // len(lst_w)
  c = lst_c[cdx]
  w = lst_w[wdx]
  mp = np.squeeze(maps[cdx, wdx, :, :, 0:4, 7])
  im = concat((concat((ic(mp[:, :, 0]), ic(mp[:, :, 1])), axis=1), \
       concat((ic(mp[:, :, 2]), ic(mp[:, :, 3])), axis=1)), axis = 0)
  im = ax.imshow(np.abs(im), vmin=0, vmax=1, cmap='gray')
  ax.set_xticks([])
  ax.set_yticks([])
  ll, bb, ww, hh = ax.get_position().bounds
  ax.set_position([ll - wdx * ww * 0.05, bb + cdx * hh * 0.3, ww, hh])
  if cdx == 0:
    ax.set_title('w = %0.1f' % (w))
  if wdx == 0:
    ax.set_ylabel('c = %0.2f' % (c))
  ctr = ctr + 1

cbar = fig.colorbar(im, ax=axes.ravel().tolist())

plt.savefig('res/variability.pdf'); plt.close(fig);

fig, axes = plt.subplots(nrows=len(lst_c), ncols=len(lst_w))
ctr = 0;
for ax in axes.flat:
  wdx = ctr % len(lst_w)
  cdx = ctr // len(lst_w)
  c = lst_c[cdx]
  w = lst_w[wdx]
  img = ic(msk * np.squeeze(igfactor2[cdx, wdx, :, :]))
  ax.imshow(np.abs(msk), cmap = 'gray')
  im  = ax.imshow(np.abs(img), cmap = 'jet', vmin=0, vmax=1)
  ax.set_xticks([])
  ax.set_yticks([])
  ll, bb, ww, hh = ax.get_position().bounds
  ax.set_position([ll - wdx * ww * 0.05, bb + cdx * hh * 0.3, ww, hh])
  txt = ax.text(90, 205, r'$g_{max}^{-1}=%0.2f$' % g(img), \
          bbox=dict(facecolor='black', alpha=0.5, pad=1), ha='center', va='center', \
          color='white', fontsize=6, fontweight='bold')
  txt = ax.text(90, 30, r'$g_{avg}^{-1}=%0.2f$' % h(gfactor2[cdx, wdx, :, :]), \
          bbox=dict(facecolor='black', alpha=0.5, pad=1), ha='center', va='center', \
          color='white', fontsize=6, fontweight='bold')
  if cdx == 0:
    ax.set_title('w = %0.1f' % (w))
  if wdx == 0:
    ax.set_ylabel('c = %0.2f' % (c))
  ctr = ctr + 1

cbar = fig.colorbar(im, ax=axes.ravel().tolist())

plt.savefig('res/variability_gfactor_rx_2_ry_2.pdf'); plt.close(fig);
