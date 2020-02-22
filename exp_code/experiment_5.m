%% EXPERIMENT 5
%
% Weighting: 
%   yes 
%
% Constant parameter:
%   kernel size
%
% Varied parameter:
%   calibration size, crop threshold

close all; clear all; clc;
addpath('../utils/', '../');
nrm = @(x) x/max(abs(x(:)));

%% Common variables

stdev =  4;

m  = readcfl('../data/maps');
X  = readcfl('../data/brain_clean');
Y  = readcfl('../data/brain_noise');
nc = size(Y, 3);
x  = nrm(sum(bsxfun(@times, F_inv(X, [1,2]), conj(m)), 3));
x  = x(:,:,:,1);
y  = F_inv(Y, [1, 2]);

lst_k = 8;
lst_c = 0.70:0.01:0.99;
lst_r = 16:2:32;

trueMSE = lst_r * 0;

for rdx=1:numel(lst_r)
  fprintf('On %02d of %02d.\n', rdx, numel(lst_r));
  fileID = fopen(sprintf('log/log_experiment_5_%d.txt', lst_r(rdx)), 'w');
  [trueMSEr, trueOptParamr, ~, ~, ~, ~, ~, ~] = ...
    compareMSE(fileID, X, Y, lst_r(rdx), true, stdev, lst_k, lst_c);
  fclose(fileID);
  trueMSE(rdx) = min(trueMSEr(:));
end

save('res/experiment_5_results.mat', 'lst_k', 'lst_c', 'lst_r', 'trueMSE');
