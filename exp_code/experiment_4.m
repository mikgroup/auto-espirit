%% EXPERIMENT 4
%
% Weighting: 
%   yes 
%
% Constant parameter:
%   kernel size
%
% Varied parameter:
%   crop threshold
%
% This experiments explores g-factor as a function of parameters.

close all; clear all; clc;
addpath('../utils/', '../');
nrm = @(x) x/max(abs(x(:)));

%% Common variables

r     = 24;
stdev =  4;

m  = readcfl('../data/maps');
X  = readcfl('../data/brain_clean');
Y  = readcfl('../data/brain_noise');
nc = size(Y, 3);
x  = nrm(sum(bsxfun(@times, F_inv(X, [1,2]), conj(m)), 3));
x  = x(:,:,:,1);
y  = F_inv(Y, [1, 2]);

fix_k = 8;
fix_c = 0.95;
fix_w = 1;

% Uncomment to test figure generation
%lst_k = 8;
%lst_c = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99];

% Uncomment to run full experiment
lst_k = 8;
lst_c = 0:0.01:0.99;

% Mask.
mask = Y(:,:,1,1) * 0;
for xdx=1:2:size(Y, 1)
  for ydx=1:2:size(Y, 2)
    mask(xdx, ydx) = 1;
  end
end
mask(floor((size(Y, 1) - r)/2) + [1:r], floor((size(Y, 2) - r)/2) + [1:r]) = 1;

% Undersampled data.
Z = bsxfun(@times, Y, mask);
z = F_inv(Z, [1, 2]);

fileID = fopen('log/log_experiment_4.txt', 'w');
[trueMSE, trueOptParam, sureFullMSE, sureFullOptParam, sureCalMSE, sureCalOptParam, gmax] = ...
  compareMSE(fileID, X, Y, r, true, stdev, fix_k, lst_c, true, true);
fclose(fileID);

rev = gmax * 0;
rev(abs(gmax(:)) > 0) = 1./gmax(abs(gmax(:)) > 0);
[~, idx] = max(rev);
[~, jdx] = min(trueMSE);
cvals = sort([0.7, lst_c(idx), lst_c(jdx)]);
C = extractCalreg(Y,     r);
A = calreg2calmat(C, fix_k); 
[U, S, V] = svd(A); s = diag(S);
weights   = svWeightsSURE(s, stdev, size(A));
if (length(weights) < size(V, 2)) 
  weights = [weights; zeros([size(V, 2) - length(weights), 1])];
end 
Vt     = V * diag(weights); 
kernel = reshape(Vt, fix_k, fix_k, nc, size(Vt, 2));
[M, W] = kernelEig(kernel, [size(Y, 1), size(Y, 2)]);

maps = [];
proj = [];
for cdx=1:numel(cvals)
  m = cropEigvec(M, W, cvals(cdx));
  p = espiritlsqr(Z(:), mask, m);
  p = nrm(p(:, :, end, end));
  maps = cat(5, maps, m);
  proj = cat(5, proj, p);
end

trueImage = x;

save('res/experiment_4_results.mat', 'lst_k', 'lst_c', 'trueMSE', 'trueOptParam', 'sureFullMSE', 'sureFullOptParam', 'sureCalMSE', 'sureCalOptParam', 'gmax', 'maps', 'proj', 'trueImage');
