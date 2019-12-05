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
 
fix_k = 6;
lst_w = [0.50, 1.00, 2.00, 4.00, 8.00];
lst_c = [0.75, 0.80, 0.85, 0.90, 0.95];

maps = zeros([length(lst_c), length(lst_w), size(X), nc]);
gfactor_rx_1_ry_2 = zeros([length(lst_c), length(lst_w), size(X, 1), size(X, 2)]);
gfactor_rx_2_ry_2 = zeros([length(lst_c), length(lst_w), size(X, 1), size(X, 2)]);

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
    m = cropEigvec(M, W, c); 
    gfactor_rx_1_ry_2(cdx, wdx, :, :) = calcGFactor(squeeze(m(:,:,:,end)), 1, 2);
    gfactor_rx_2_ry_2(cdx, wdx, :, :) = calcGFactor(squeeze(m(:,:,:,end)), 2, 2);
    maps(cdx, wdx, :, :, :, :) = m;
  end
end

save('res/map_dem.mat', 'fix_k', 'lst_c', 'lst_w', 'maps', 'gfactor_rx_1_ry_2', 'gfactor_rx_2_ry_2');

% Reference: /https://web.stanford.edu/class/ee369c/restricted/Solutions/assignment_4_solns.pdf
function gmap = calcGFactor(m, rx, ry)
  gmap = zeros(size(m, 1), size(m, 2));
  assert(mod(size(m, 1), rx) == 0);
  assert(mod(size(m, 2), ry) == 0);
  shift_x = size(m, 1)/rx;
  shift_y = size(m, 2)/ry;
  f = @(x) x(:);
  mat = Inf;
  for x=1:1:size(m, 1)
    for y=1:1:size(m, 2)
      for xx=0:(rx-1)
        for yy=0:(ry-1)
          if (xx == 0 && yy == 0)
            mat = f(m(x, y, :));
          else
            xdx = mod((x - 1) + xx * shift_x, size(m, 1)) + 1;
            ydx = mod((y - 1) + yy * shift_y, size(m, 2)) + 1;
            mat = [mat, f(m(xdx, ydx, :))];
          end
        end
      end
      scs  = mat' * mat;
      scsi = pinv(scs);
      cnd  = sqrt(scsi(1, 1) * scs(1, 1));
      gmap(x, y) = abs(cnd);
    end
  end
end
