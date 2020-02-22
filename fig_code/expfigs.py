#!/usr/bin/env python3

import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from matplotlib import colors, ticker, cm
from scipy.io import loadmat
from math import ceil, floor, log10
from os.path import isfile

DPI = 100
matplotlib.rcParams.update({'font.size': 32})

concat = np.concatenate

# Experiment 1

dct = loadmat('../exp_code/res/experiment_1_results.mat')

lst_k = np.squeeze(dct['lst_k'])
lst_c = np.squeeze(dct['lst_c'])
lst_w = np.squeeze(dct['lst_w'])

trueMSE     = np.squeeze(dct['trueMSE'])
sureFullMSE = np.squeeze(dct['sureFullMSE'])
sureCalMSE  = np.squeeze(dct['sureCalMSE'])

minVal      = min(trueMSE.min(), sureFullMSE.min());
sureCalMSE  = sureCalMSE - sureCalMSE.min() + minVal;
maxVal      = max(trueMSE.max(), sureCalMSE.max(), sureFullMSE.max())

levels  = np.arange(np.floor(np.log10(minVal) * 10)/10, np.ceil(np.log10(maxVal) * 10)/10, 0.001)

if not isfile('res/experiment_1_3.png'):
  fig = plt.figure(figsize=(15, 15))
  cp = plt.contourf(lst_w, lst_c, np.log10(trueMSE), levels=levels, cmap='jet', vmin=levels[0], vmax=levels[-1])
  plt.title('True Squared Error'); plt.xlabel('Subspace size (WNSVN)'); plt.ylabel('Crop threshold(c)');
  cb = plt.colorbar()
  cb.ax.set_yticklabels([(r'$10^{%0.2f}$' % float(elm.get_text())) for elm in cb.ax.get_yticklabels()])
  plt.savefig('res/experiment_1_1.png', dpi=DPI); plt.close(fig)

  fig = plt.figure(figsize=(15, 15))
  cp = plt.contourf(lst_w, lst_c, np.log10(sureFullMSE), levels=levels, cmap='jet', vmin=levels[0], vmax=levels[-1])
  plt.title('SURE (Full Data)'); plt.xlabel('Subspace size (WNSVN)'); plt.ylabel('Crop threshold(c)');
  cb = plt.colorbar()
  cb.ax.set_yticklabels([(r'$10^{%0.2f}$' % float(elm.get_text())) for elm in cb.ax.get_yticklabels()])
  plt.savefig('res/experiment_1_2.png', dpi=DPI); plt.close(fig)

  fig = plt.figure(figsize=(15, 15))
  cp = plt.contourf(lst_w, lst_c, np.log10(sureCalMSE), levels=levels, cmap='jet', vmin=levels[0], vmax=levels[-1])
  plt.title('Normalized SURE (Calib. Data)'); plt.xlabel('Subspace size (WNSVN)'); plt.ylabel('Crop threshold(c)');
  cb = plt.colorbar()
  cb.ax.set_yticklabels([(r'$10^{%0.2f}$' % float(elm.get_text())) for elm in cb.ax.get_yticklabels()])
  plt.savefig('res/experiment_1_3.png', dpi=DPI); plt.close(fig)

# Experiment 2

dct = loadmat('../exp_code/res/experiment_2_results.mat')

lst_k = np.squeeze(dct['lst_k'])
lst_c = np.squeeze(dct['lst_c'])
lst_w = np.squeeze(dct['lst_w'])

trueMSE     = np.squeeze(dct['trueMSE'])
sureFullMSE = np.squeeze(dct['sureFullMSE'])
sureCalMSE  = np.squeeze(dct['sureCalMSE'])

minVal      = min(trueMSE.min(), sureFullMSE.min());
sureCalMSE  = sureCalMSE - sureCalMSE.min() + minVal;
maxVal      = max(trueMSE.max(), sureCalMSE.max(), sureFullMSE.max())

if True:
  fig = plt.figure(figsize=(15, 15))
  ax  = plt.gca();
  ax.plot(lst_c, np.log10(trueMSE),     'r-',  linewidth=5);
  ax.plot(lst_c, np.log10(sureFullMSE), 'b--', linewidth=5, alpha=0.75);
  ax.plot(lst_c, np.log10(sureCalMSE),  'g--', linewidth=5, alpha=0.75);
  plt.savefig('res/experiment_2.png', dpi=10); # This is to force text to be printed.
  ax.set_yticklabels([(r'$10^{%0.2f}$' % float(elm.get_text())) for elm in ax.get_yticklabels()])
  plt.grid(True, linestyle=':')
  plt.xlabel('Crop threshold (c)'); plt.ylabel('Squared Error');
  plt.legend(['True Squared Error', 'SURE (Full Data)', 'Normalized SURE (Calib. Data)']); 
  plt.savefig('res/experiment_2.png', dpi=DPI); plt.close(fig);

# Experiment 3.1

dct = loadmat('../exp_code/res/experiment_3_results_noweight.mat')

lst_k = np.squeeze(dct['lst_k'])

trueMSE     = np.squeeze(dct['trueMSE'])
sureFullMSE = np.squeeze(dct['sureFullMSE'])
sureCalMSE  = np.squeeze(dct['sureCalMSE'])

trueMSE     = trueMSE.min(axis=(1, 2))
sureFullMSE = sureFullMSE.min(axis=(1, 2))
sureCalMSE  = sureCalMSE.min(axis=(1, 2))

minVal      = min(trueMSE.min(), sureFullMSE.min());
sureCalMSE  = sureCalMSE - sureCalMSE.min() + minVal;
maxVal      = max(trueMSE.max(), sureCalMSE.max(), sureFullMSE.max())

if not isfile('res/experiment_3_1.png'):
  fig = plt.figure(figsize=(15, 15))
  ax  = plt.gca()
  ax.plot(lst_k, np.log10(trueMSE),     'r-',  linewidth=5);
  ax.plot(lst_k, np.log10(sureFullMSE), 'b--', linewidth=5, alpha=0.75);
  ax.plot(lst_k, np.log10(sureCalMSE),  'g--', linewidth=5, alpha=0.75);
  plt.savefig('res/experiment_3_1.png', dpi=10); # This is to force text to be printed.
  ax.set_yticklabels([(r'$10^{%0.2f}$' % float(elm.get_text())) for elm in ax.get_yticklabels()])
  plt.legend(['True Squared Error', 'SURE (Full Data)', 'Normalized SURE (Calib. Data)']); 
  plt.grid(True, linestyle=':')
  plt.xlabel('Kernel size (k)'); plt.ylabel('Squared Error');
  plt.savefig('res/experiment_3_1.png', dpi=DPI); plt.close(fig);

# Experiment 3.2

dct = loadmat('../exp_code/res/experiment_3_results_weight.mat')

lst_k = np.squeeze(dct['lst_k'])

trueMSE     = np.squeeze(dct['trueMSE'])
sureFullMSE = np.squeeze(dct['sureFullMSE'])
sureCalMSE  = np.squeeze(dct['sureCalMSE'])

trueMSE     = trueMSE.min(axis=(1))
sureFullMSE = sureFullMSE.min(axis=(1))
sureCalMSE  = sureCalMSE.min(axis=(1))

minVal      = min(trueMSE.min(), sureFullMSE.min());
sureCalMSE  = sureCalMSE - sureCalMSE.min() + minVal;
maxVal      = max(trueMSE.max(), sureCalMSE.max(), sureFullMSE.max())

if not isfile('res/experiment_3_2.png'):
  fig = plt.figure(figsize=(15, 15))
  ax  = plt.gca()
  plt.plot(lst_k, np.log10(trueMSE),     'r-',  linewidth=5)
  plt.plot(lst_k, np.log10(sureFullMSE), 'b--', linewidth=5, alpha=0.75)
  plt.plot(lst_k, np.log10(sureCalMSE),  'g--', linewidth=5, alpha=0.75)
  plt.savefig('res/experiment_3_2.png', dpi=10); # This is to force text to be printed.
  ax.set_yticklabels([(r'$10^{%0.2f}$' % float(elm.get_text())) for elm in ax.get_yticklabels()])
  plt.legend(['True Squared Error', 'SURE (Full Data)', 'Normalized SURE (Calib. Data)']); 
  plt.grid(True, linestyle=':')
  plt.xlabel('Kernel size (k)'); plt.ylabel('Squared Error');
  plt.savefig('res/experiment_3_2.png', dpi=DPI); plt.close(fig);

# Experiment 4
dct = loadmat('../exp_code/res/experiment_4_results.mat')
nrm = lambda x: (x - np.min(x))/np.max(x - np.min(x))

lst_k      = np.squeeze(dct['lst_k'])
lst_c      = np.squeeze(dct['lst_c'])
gmax       = np.squeeze(dct['gmax'])
gavg       = np.squeeze(dct['gavg'])
proj       = np.squeeze(dct['proj'])
trueMSE    = nrm(np.squeeze(dct['trueMSE']))
sureFulMSE = nrm(np.squeeze(dct['sureFullMSE']))
sureCalMSE = nrm(np.squeeze(dct['sureCalMSE']))
revmax     = gmax * 0
revmax[np.abs(gmax) > 0] = 1/gmax[np.abs(gmax) > 0];
revavg     = gavg * 0
revavg[np.abs(gavg) > 0] = 1/gavg[np.abs(gavg) > 0];

print("Max 1/g_max:    %f" % np.max(revmax));
print("Location (c):   %f" % lst_c[np.argmax(revmax)]);
print("Min sureCalMSE: %f" % np.min(sureCalMSE));
print("Location (c):   %f" % lst_c[np.argmin(sureCalMSE)]);
print("Min trueMSE:    %f" % np.min(trueMSE));
print("Location (c):   %f" % lst_c[np.argmin(trueMSE)]);

if not isfile('res/experiment_4_2.png'):
  fig, ax1 = plt.subplots(figsize=(13, 10))

  ax1.set_xlabel('Crop threshold (c)')
  ax1.set_ylabel(r'$1/g$', color='black')
  line1 = ax1.plot(lst_c, revmax, color='tab:red', linewidth=5)
  line2 = ax1.plot(lst_c, revavg, color='lime',    linewidth=5)
  ax1.tick_params(axis='y', labelcolor='black')
  ax1.set_xlim([-0.1, 1.1])
  ax1.set_ylim([-0.1, 1.1])

  red  = mpatches.Patch(color='red',   label=r'$1/g_{max}$')
  lime = mpatches.Patch(color='lime', label=r'$1/g_{avg}$')

  ax1.legend(handles=[red, lime])

  #ax1.set_xticks(np.array([0, 0.25, 0.5, 0.75, 1]) * np.max(revmax))
  #ax1.set_xticklabels(np.array([0, 0.25, 0.5, 0.75, 1]) * np.max(revmax))

  ax2 = ax1.twinx()

  color = 'tab:blue'
  ax2.set_ylabel('Normalized Squared Error', color=color)
  ax2.plot(lst_c, trueMSE, color=color, linewidth=5)
  ax2.tick_params(axis='y', labelcolor=color)
  ax2.axvline(x=lst_c[np.argmin(trueMSE)], color='g', linestyle='--', alpha=0.75, linewidth=3)
  ax2.axvline(x=lst_c[np.argmax(revmax)], color='goldenrod', linestyle='--', alpha=0.75, linewidth=3)
  ax2.axvline(x=0.7, color='m', linestyle='--', alpha=0.75, linewidth=3)
  #ax2.set_xticks(np.array([0, 0.25, 0.5, 0.75, 1]) * np.max(trueMSE))
  #ax2.set_xticklabels(np.array([0, 0.25, 0.5, 0.75, 1]) * np.max(trueMSE))
  ax2.set_xlim([-0.1, 1.1])
  ax2.set_ylim([-0.1, 1.1])

  ax2.plot(lst_c[np.argmin(trueMSE)], trueMSE[np.argmin(trueMSE)], 'gx', mew=2, ms=15, linewidth=3)
  ax1.plot(lst_c[np.argmin(trueMSE)],  revmax[np.argmin(trueMSE)], 'gx', mew=2, ms=15, linewidth=3)
  ax1.plot(lst_c[np.argmin(trueMSE)],  revavg[np.argmin(trueMSE)], 'gx', mew=2, ms=15, linewidth=3)

  ax2.plot(lst_c[np.argmax(revmax)], trueMSE[np.argmax(revmax)], color='orange', marker='x', mew=2, ms=15, linewidth=3)
  ax1.plot(lst_c[np.argmax(revmax)],  revmax[np.argmax(revmax)], color='orange', marker='x', mew=2, ms=15, linewidth=3)
  ax1.plot(lst_c[np.argmax(revmax)],  revavg[np.argmax(revmax)], color='orange', marker='x', mew=2, ms=15, linewidth=3)

  ax2.plot(0.7, trueMSE[lst_c == 0.7], color='m', marker='x', mew=2, ms=15, linewidth=3)
  ax1.plot(0.7,  revmax[lst_c == 0.7], color='m', marker='x', mew=2, ms=15, linewidth=3)
  ax1.plot(0.7,  revavg[lst_c == 0.7], color='m', marker='x', mew=2, ms=15, linewidth=3)

  ax1.xaxis.grid(True, linestyle='--')
  ax1.yaxis.grid(True, linestyle='--')

  plt.savefig('res/experiment_4_1.png', dpi=DPI); plt.close(fig);

  fig, ax = plt.subplots(figsize=(10, 9))

  proj[:, proj.shape[1]//2:, :] = 4 * proj[:, proj.shape[1]//2:, :]
  proj = concat((proj[:,:,0], proj[:,:,1], proj[:,:,2]),axis=1);

  ax.imshow(np.abs(proj), cmap='gray', vmax=0.75)
  ax.axis('off')

  plt.savefig('res/experiment_4_2.png', dpi=DPI); plt.close(fig);

# Experiment 5

dct = loadmat('../exp_code/res/experiment_5_results.mat')

lst_r   = np.squeeze(dct['lst_r'])
trueMSE = np.squeeze(dct['trueMSE'])

if not isfile('res/experiment_5.png'):
  fig, ax = plt.subplots(figsize=(12, 10))
  ax.plot(lst_r, trueMSE, 'b-', linewidth=5)
  ax.set_xlabel('Calibration size (r)')
  ax.set_ylabel('True Squared Error');
  ax.ticklabel_format(axis='y', scilimits=(0, 0), useMathText=True)
  ax.grid()
  plt.savefig('res/experiment_5.png', dpi=DPI); plt.close(fig);

# Experiment 6

dct = loadmat('../exp_code/res/experiment_6_results.mat')

lst_k     = np.squeeze(dct['lst_k'])
revmax2x1 = np.squeeze(dct['revmax2x1'])
revmax2x2 = np.squeeze(dct['revmax2x2'])
revavg2x1 = np.squeeze(dct['revavg2x1'])
revavg2x2 = np.squeeze(dct['revavg2x2'])

if not isfile('res/experiment_6.png'):
  fig, ax = plt.subplots(figsize=(15, 12))
  ax.plot(lst_k, revavg2x1, 'r-', linewidth=5)
  ax.plot(lst_k, revavg2x2, 'b-', linewidth=5)
  ax.plot(lst_k, revmax2x1, 'r:', linewidth=5)
  ax.plot(lst_k, revmax2x2, 'b:', linewidth=5)

  ax.set_xlabel('Kernel Size (k)')
  ax.set_ylabel(r'$1/g$');
  ax.set_ylim([0, 1])
  ax.ticklabel_format(axis='y', scilimits=(0, 0), useMathText=True)
  ax.grid()
  ax.legend([r'$R = (2 \times 1)\;\rightarrow\;1/g_{avg}$', \
             r'$R = (2 \times 2)\;\rightarrow\;1/g_{avg}$', \
             r'$R = (2 \times 1)\;\rightarrow\;1/g_{max}$', \
             r'$R = (2 \times 2)\;\rightarrow\;1/g_{max}$'])
  plt.savefig('res/experiment_6.png', dpi=DPI); plt.close(fig);
