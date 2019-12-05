%% EXPERIMENT 6
%
% Weighting: 
%   yes 
%
% Varied parameter:
%   kernel size
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
M  = readcfl('../data/brain_FOV');
nc = size(Y, 3);
x  = nrm(sum(bsxfun(@times, F_inv(X, [1,2]), conj(m)), 3));
x  = x(:,:,:,1);
y  = F_inv(Y, [1, 2]);

lst_k   = 4:1:12;
gmax2x1 = 0 * lst_k;
gmax2x2 = 0 * lst_k;
gavg2x1 = 0 * lst_k;
gavg2x2 = 0 * lst_k;

Y = reshape(Y, size(Y, 1), size(Y, 2), 1, size(Y, 3));
ctr = 0;
for k=lst_k
  ctr = ctr + 1;
  maps = squeeze(bart(sprintf('ecalib -v %f -k %d -c %d -a', stdev^2, k, r), Y)); % Using BART as it is faster.
  g1 = calcGFactor(maps, 2, 1);
  g2 = calcGFactor(maps, 2, 2);
  g1 = g1(M(:) > 0.5);
  g2 = g2(M(:) > 0.5);
  gmax2x1(ctr) = max(g1(:));
  gmax2x2(ctr) = max(g2(:));
  gavg2x1(ctr) = mean(g1(:));
  gavg2x2(ctr) = mean(g2(:));
end

revmax2x1 = 1./gmax2x1;
revmax2x2 = 1./gmax2x2;
revavg2x1 = 1./gavg2x1;
revavg2x2 = 1./gavg2x2;

save('res/experiment_6_results.mat', 'lst_k', 'r',      ...
     'gmax2x1',   'gmax2x2',   'gavg2x1',   'gavg2x2',  ... 
     'revmax2x1', 'revmax2x2', 'revavg2x1', 'revavg2x2');
