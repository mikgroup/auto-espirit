%% EXPERIMENT 1
%
% Weighting: 
%   no 
%
% Constant parameters:
%   kernel size
%
% Varied parameter:
%   subspace size, crop threshold

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
lst_k = 3:1:10;
lst_c = 0:0.01:0.99;
lst_w = 0:1/(fix_k^2):nc;

fileID = fopen('log/log_experiment_1.txt', 'w');
[trueMSE, trueOptParam, sureFullMSE, sureFullOptParam, sureCalMSE, sureCalOptParam, ~] = ...
  compareMSE(fileID, X, Y, r, false, stdev, fix_k, lst_c, lst_w);
fclose(fileID);

save('res/experiment_1_results.mat', 'lst_k', 'lst_c', 'lst_w', 'trueMSE', 'trueOptParam', 'sureFullMSE', 'sureFullOptParam', 'sureCalMSE', 'sureCalOptParam');
