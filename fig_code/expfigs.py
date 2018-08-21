#!/usr/bin/env python3

import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import colors, ticker, cm
from scipy.io import loadmat
from math import ceil, floor, log10

matplotlib.rcParams.update({'font.size': 18})

## Helper code.
#
# http://stackoverflow.com/a/25983372

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
fmt = ticker.FuncFormatter(fmt)

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

step   = 10 ** (floor(log10(minVal))-1);
levels = list(range(floor(minVal), ceil(maxVal + step), step))

fig = plt.figure(figsize=(10, 8))
cp = plt.contourf(lst_w, lst_c, trueMSE, levels = levels, cmap = plt.get_cmap("jet"))
plt.title('True Squared Error'); plt.xlabel('Subspace size (WNSVN)'); plt.ylabel('Crop threshold(c)');
#cb = plt.colorbar(format=fmt)
cb = plt.colorbar()
plt.savefig('pdf/experiment_1_1.pdf'); plt.close(fig)

fig = plt.figure(figsize=(10, 8))
cp = plt.contourf(lst_w, lst_c, sureFullMSE, levels = levels, cmap = plt.get_cmap("jet"))
plt.title('SURE (Full Data)'); plt.xlabel('Subspace size (WNSVN)'); plt.ylabel('Crop threshold(c)');
#cb = plt.colorbar(format=fmt)
cb = plt.colorbar()
plt.savefig('pdf/experiment_1_2.pdf'); plt.close(fig)

fig = plt.figure(figsize=(10, 8))
cp = plt.contourf(lst_w, lst_c, sureCalMSE, levels = levels, cmap = plt.get_cmap("jet"))
plt.title('Normalized SURE (Calib. Data)'); plt.xlabel('Subspace size (WNSVN)'); plt.ylabel('Crop threshold(c)');
#cb = plt.colorbar(format=fmt)
cb = plt.colorbar()
plt.savefig('pdf/experiment_1_3.pdf'); plt.close(fig)

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

fig = plt.figure(figsize=(10, 9))
#plt.semilogy(lst_c, trueMSE, 'r-'); plt.semilogy(lst_c, sureFullMSE, 'g+'); plt.semilogy(lst_c, sureCalMSE, 'bo');
plt.plot(lst_c, trueMSE, 'r-'); plt.plot(lst_c, sureFullMSE, 'g+'); plt.plot(lst_c, sureCalMSE, 'bo');
#plt.legend(['True Squared Error', 'SURE (Full Data)', 'Normalized SURE (Calib. Data)']); 
plt.xlabel('Crop threshold (c)'); plt.ylabel('Squared Error');
plt.savefig('pdf/experiment_2.pdf'); plt.close(fig);

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

fig = plt.figure(figsize=(10, 9))
#plt.semilogy(lst_k, trueMSE, 'r-'); plt.semilogy(lst_k, sureFullMSE, 'g+'); plt.semilogy(lst_k, sureCalMSE, 'bo');
plt.plot(lst_k, trueMSE, 'r-'); plt.plot(lst_k, sureFullMSE, 'g+'); plt.plot(lst_k, sureCalMSE, 'bo');
plt.ylim((0, 1E7));
#plt.legend(['True Squared Error', 'SURE (Full Data)', 'Normalized SURE (Calib. Data)']); 
plt.xlabel('Kernel size (k)'); plt.ylabel('Squared Error');
plt.savefig('pdf/experiment_3_1.pdf'); plt.close(fig);

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

fig = plt.figure(figsize=(10, 9))
#plt.semilogy(lst_k, trueMSE, 'r-'); plt.semilogy(lst_k, sureFullMSE, 'g+'); plt.semilogy(lst_k, sureCalMSE, 'bo');
plt.plot(lst_k, trueMSE, 'r-'); plt.plot(lst_k, sureFullMSE, 'g+'); plt.plot(lst_k, sureCalMSE, 'bo');
plt.ylim((0, 1E7));
plt.legend(['True Squared Error', 'SURE (Full Data)', 'Normalized SURE (Calib. Data)']); 
plt.xlabel('Kernel size (k)'); plt.ylabel('Squared Error');
plt.savefig('pdf/experiment_3_2.pdf'); plt.close(fig);
