%% Map variability
%
% Illustrates the variability of estimated maps as a function of parameters.

close all; clear all; clc;
addpath('../utils/', '../');

r      = 24;
stdev  = 10;

X = readcfl('../data/brain_clean');
Y = readcfl('../data/brain_noise');
nc = size(X, 3); 
 
fix_k = 3;
lst_w = [0.50, 1.00, 2.00, 4.00, 8.00];
lst_c = [0.75, 0.8, 0.85, 0.9, 0.95];

maps = zeros([length(lst_c), length(lst_w), size(X), nc]);

C = extractCalreg(X, r);
A = calreg2calmat(C, fix_k); [U, S, V] = svd(A); s = diag(S);

wdx = 0;
for w = lst_w
    wdx = wdx + 1;

    n = floor(w * fix_k * fix_k); 
    vecWeight = ones([n, 1]);
    vecWeight = [vecWeight; zeros([size(V, 2) - n, 1])];
    wV = V * diag(vecWeight);
    kernel = reshape(wV, fix_k, fix_k, nc, size(wV, 2));
    [M, W] = kernelEig(kernel, [size(X, 1), size(X, 2)]);

    cdx = 0;
    for c = lst_c
        cdx = cdx + 1;
        maps(cdx, wdx, :, :, :, :) = cropEigvec(M, W, c); 
    end

end

save('res/map_dem.mat', 'fix_k', 'lst_c', 'lst_w', 'maps');
