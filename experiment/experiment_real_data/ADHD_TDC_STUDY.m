%% Real data experiment
% The goal of this experiment is to study the effective connectivity
% differences between ADHD and TDC by using proposed formulations.
% Setting
% - use cvx-DGN, cvx-FGN to estimate concatenation of K=18 subjects in each of ADHD and TDC data sets
% - use cvx-CGN to estimate two models from K=18 ADHD and TDC respectively
% the common part of GC matrix in each model will be considered as group
% level GC.

%% Data concatenation
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0); % load data in format (n,T,K)
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
y_total = cat(3,y_TDC,y_ADHD_C);
y_total = y_total-mean(y_total,2); % detrend in time axis
%% cvx-DGN estimation
p = 1; % VAR order
GridSize = 30; % resolution of regularization grid
data_concat = 1; % set to 1 for time-series concatenation without dependency 
toggle = 'adaptive_L'; % set weight, toggle = 'static' is to set weight to unity
M = test_cvxformulation_D(y_total,p,GridSize,toggle,data_concat);% data with dimension (n,T*K,2), K is # subjects in each TDC, ADHD
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected','M')
%% cvx-FGN estimation
p = 1;
GridSize = 30;
data_concat = 1;
toggle = 'adaptive_D';
M = test_cvxformulation_S(y_total,p,GridSize,toggle,data_concat);% data with dimension (n,T*K,2), K is # subjects in each TDC, ADHD
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected','M')
%% cvx-CGN estimation
clear M
p = 1; % VAR order
GridSize = 30;
data_concat = 1;
toggle = 'adaptive_L';
M.TDC = test_cvxformulation_C(y_TDC,p,GridSize,toggle); % data with dimension (n,p,K)
M.ADHD_C = test_cvxformulation_C(y_ADHD_C,p,GridSize,toggle);
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered','M')
%% Perform model selection correction by setting noise correlation structure to be diagonal: FGN
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected','M')
p = 1;
GridSize = 30;
data_concat = 1;
toggle = 'adaptive_D';
M = correction_S(y_total,M,'LLH_hetero',p,GridSize,toggle,data_concat);
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected_LLHcorrected','M')
%% Perform model selection correction by setting noise correlation structure to be diagonal: DGN
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected','M')
p = 1;
GridSize = 30;
data_concat = 1;
toggle = 'adaptive_L';
M = correction_S(y_total,M,'LLH_hetero',p,GridSize,toggle,data_concat);
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected_LLHcorrected','M')