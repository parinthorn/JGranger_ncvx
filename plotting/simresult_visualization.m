%% This script is for printing result
clear
clc
clf
close all
performance_path = './experiment/result_to_plot/';
load([performance_path,'summaryK5'])

experiment_name = {'magda','C_cvx','C'; ...
                   'skripD','D_cvx','D'; ...
                   'skripS','S_cvx','S'};
for exp=1:3
    for type=1:3
        
        
    end
end

% COMPARISON WITH OTHER WORK
% Formulation C
% 1. Magda
% 2. Proposed Convex
% 3. Proposed Non-convex



% Formulation D
% 1. Skripnikov
% 2. Proposed Convex
% 3. Proposed Non-convex


% Formulation S
% 1. Skripnikov
% 2. Proposed Convex
% 3. Proposed Non-convex


%%
% Formulation D K=5
% 1. Skripnikov
% 2. Proposed Convex
% 3. Proposed Non-convex
% Formulation D K=50
% 1. Skripnikov
% 2. Proposed Convex
% 3. Proposed Non-convex
%%
% SAMPLE GC VISUALIZATION(Maximum performance in K=5, Circular graph, black=common, red=differential)
% Formulation C
% 1. Magda
% 2. Proposed Convex
% 3. Proposed Non-convex
% Formulation D
% 1. Skripnikov
% 2. Proposed Convex
% 3. Proposed Non-convex
% Formulation S
% 1. Skripnikov
% 2. Proposed Convex
% 3. Proposed Non-convex



%%
% CONVEX NON-CONVEX COMPARISON
% Formulation C

% Formulation D

% Formulation S


%% reformatting result
% clear
% clc
% main_dir = './experiment/result_to_plot/';
% toggle = {'total','common','differential'};
% 
% f_list = {'formulation_C','formulation_D', ...
%     'formulation_S','formulation_C_cvx', ...
%     'formulation_D_cvx','formulation_S_cvx'};
% 
% for f=1:length(f_list)
%     fname = f_list{f};
%     allname = ['adaptive_',fname,'_ALL_RESULT_K5'];
%     indexname= ['adaptive_',fname,'_result_K5'];
%     acc_list = {'TPR','FPR','ACC','F1','MCC'};
%     load([main_dir,allname])
%     load([main_dir,indexname])
%     
%     for ii =1:2
%         for jj=1:size(R.index,2)
%             for t=1:3
%                 for kk=1:length(acc_list)
%                     tmp.(toggle{t}).(acc_list{kk})(ii,jj) = ALL_RESULT(ii,jj).model_acc(R.index(ii,jj).eBIC).(toggle{t}).(acc_list{kk});
%                 end
%                 
%             end
%             tmp.bias(ii,jj) = ALL_RESULT(ii,jj).model_acc(R.index(ii,jj).eBIC).bias;
%         end
%     end
%     clear R
%     R = tmp;
%     save([main_dir,fname,'_eBICresult'],'R')
% end
%% skripnikov D cleanup K=5
% clear
% clc
% 
% main_dir = './experiment/result_to_plot/';
% fname = 'skrip_formulationS_accuracy_K5.mat';
% load([main_dir,fname])
% toggle = {'total','common','differential'};
% acc_list = {'TPR','FPR','ACC','F1','MCC'};
% for ii=1:2
%         for t=1:length(toggle)
%             for s=1:length(acc_list)
%         R.(toggle{t}).(acc_list{s})(ii,:) = score(ii).(toggle{t}).(acc_list{s});
%             end
%         end
% end
% save([main_dir,'skripS_result'],'R')
%% all_result summary file

% clear
% clc
% % COMPARISON WITH OTHER WORK
% % Formulation C
% % 1. Magda
% % 2. Proposed Convex
% % 3. Proposed Non-convex
% 
% load('./experiment/result_to_plot/magda_result.mat')
% result.magda = R;
% load('./experiment/result_to_plot/formulation_C_cvx_eBICresult.mat')
% result.C_cvx = R.common;
% load('./experiment/result_to_plot/formulation_C_eBICresult.mat')
% result.C = R.common;
% 
% 
% % Formulation D
% % 1. Skripnikov
% % 2. Proposed Convex
% % 3. Proposed Non-convex
% 
% load('./experiment/result_to_plot/skripD_result.mat')
% result.skripD = R;
% load('./experiment/result_to_plot/formulation_D_cvx_eBICresult.mat')
% result.D_cvx = R;
% load('./experiment/result_to_plot/formulation_D_eBICresult.mat')
% result.D = R;
% 
% 
% % Formulation S
% % 1. Skripnikov
% % 2. Proposed Convex
% % 3. Proposed Non-convex
% load('./experiment/result_to_plot/skripS_result.mat')
% result.skripS = R;
% load('./experiment/result_to_plot/formulation_S_cvx_eBICresult.mat')
% result.S_cvx = R;
% load('./experiment/result_to_plot/formulation_S_eBICresult.mat')
% result.S = R;
% 
% save('./experiment/result_to_plot/summaryK5.mat','result')