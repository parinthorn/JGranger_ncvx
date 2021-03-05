%% This script investigate the differences of old version of code.
% load data
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0); % load data in format (n,T,K)
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
T = size(y_TDC,2);
y_total = cat(3,y_TDC,y_ADHD_C);
y_total = y_total-mean(y_total,2); % detrend in time axis
%% Perform model selection correction by setting noise correlation structure to be diagonal: FGN
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected','M')
p = 1;
GridSize = 30;
data_concat = 1;
toggle = 'adaptive_D';
M_old = augment_score(M,T*K-p*K,'LLH_hetero');
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected_LLHcorrected','M')
% M_new = correction_S(y_total,M,'LLH_hetero',p,GridSize,toggle,data_concat);
M_new = M;
GC.new = M_new.model(M_new.index.eBIC).GC(:);
GC.old = M_old.model(M_old.index.eBIC).GC(:);
extra = length(setdiff(find(GC.new),find(GC.old)));
missing = length(setdiff(find(GC.old),find(GC.new)));
fprintf('#extra link:%d, #missing link:%d\n',extra,missing)

%% Perform model selection correction by setting noise correlation structure to be diagonal: DGN
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected','M')
p = 1;
GridSize = 30;
data_concat = 1;
toggle = 'adaptive_L';
M_old = augment_score(M,T*K-p*K,'LLH_hetero');
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected_LLHcorrected','M')
% M_new = correction_S(y_total,M,'LLH_hetero',p,GridSize,toggle,data_concat);
M_new = M;
GC.new = M_new.model(M_new.index.eBIC).GC(:);
GC.old = M_old.model(M_old.index.eBIC).GC(:);
extra = length(setdiff(find(GC.new),find(GC.old)));
missing = length(setdiff(find(GC.old),find(GC.new)));
fprintf('#extra link:%d, #missing link:%d\n',extra,missing)
%% Conclusion

% This script confirmed that the version of xxx_timecorrected_LLHcorrected.mat is equivalent to  xxx_timecorrected.mat with LLH_hetero



