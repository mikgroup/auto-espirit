%% EXPERIMENT 3
%
% Experiment 1 and Experiment 2 repeated over different kernel sizes. The minimum MSE over other parameters is taken to plot
% MSE as a function of kernel size.

close all; clear all; clc;
addpath('../utils/', '../');

%% Common variables

r     = 24;
stdev =  4;

X  = readcfl('../data/brain_clean');
Y  = readcfl('../data/brain_noise');
nc = size(X, 3);

fix_k = 6;
fix_c = 0.95;
fix_w = 1;

% Uncomment to test figure generation
%lst_k = [4, 6];
%lst_c = 0.95:0.02:0.99;
%lst_w = [0.5, 1];

% Uncomment to run full experiment
lst_k = 4:1:10;
lst_c = 0:0.05:0.95;
lst_w = 0:0.25:nc;

fileID = fopen('log/log_experiment_3_noweight.txt', 'w');
[trueMSE, trueOptParam, sureFullMSE, sureFullOptParam, sureCalMSE, sureCalOptParam, ~, ~] = ...
  compareMSE(fileID, X, Y, r, false, stdev, lst_k, lst_c, lst_w);
fclose(fileID);

trueMSE_k = zeros([length(lst_k), 1]);
sureFullMSE_k = zeros([length(lst_k), 1]);
sureCalMSE_k = zeros([length(lst_k), 1]);

for idx=1:1:length(lst_k)
  tmp = trueMSE(idx, :, :);
  trueMSE_k(idx) = min(tmp(:));

  tmp = sureFullMSE(idx, :, :);
  sureFullMSE_k(idx) = min(tmp(:));

  tmp = sureCalMSE(idx, :, :);
  sureCalMSE_k(idx) = min(tmp(:));
end

save('res/experiment_3_results_noweight.mat', 'lst_k', 'lst_c', 'lst_w', 'trueMSE', 'trueOptParam', 'sureFullMSE', 'sureFullOptParam', 'sureCalMSE', 'sureCalOptParam');

fileID = fopen('log/log_experiment_3_weight.txt', 'w');
[trueMSE, trueOptParam, sureFullMSE, sureFullOptParam, sureCalMSE, sureCalOptParam, ~, ~] = ...
  compareMSE(fileID, X, Y, r, true, stdev, lst_k, lst_c);
fclose(fileID);

trueMSE_k = zeros([length(lst_k), 1]);
sureFullMSE_k = zeros([length(lst_k), 1]);
sureCalMSE_k = zeros([length(lst_k), 1]);

for idx=1:1:length(lst_k)
  tmp = trueMSE(idx, :, :);
  trueMSE_k(idx) = min(tmp(:));

  tmp = sureFullMSE(idx, :, :);
  sureFullMSE_k(idx) = min(tmp(:));

  tmp = sureCalMSE(idx, :, :);
  sureCalMSE_k(idx) = min(tmp(:));
end

save('res/experiment_3_results_weight.mat', 'lst_k', 'lst_c', 'lst_w', 'trueMSE', 'trueOptParam', 'sureFullMSE', 'sureFullOptParam', 'sureCalMSE', 'sureCalOptParam');
