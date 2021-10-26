clear
clc
inpath = './experiment/model_parameters/';
resultpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx\';
performance_path = './results2plot/';
mname = {'1','5'};
dd = length(mname);
% dd=1;
% dd=2;
K=5;
p_true = 3;
T = 150;

load([inpath,'compare_convex_model_K',int2str(K),'_p',int2str(p_true)]) % struct E
name_list = {'bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
ii=2;
%%
realization = 100;
for jj=1:realization
    fprintf('(%d,%d)\n',ii,jj)
    GTmodel = E{2,jj};
    fname = [resultpath,'estim_CGN_T150_',mname{ii},'percent_lag3_K',int2str(K),'_',int2str(jj)];
    load(fname)
%     M = augment_score_old(M,T);
    model_acc = performance_eval(M,GTmodel);
    toggle_list = {'total','common','differential'};
    %         M.index.bic=best_index(jj);
    ALL_RESULT(ii,jj).model_acc = model_acc;
    for kk=1:length(name_list)
        R.index(ii,jj).(name_list{kk}) = M.index.(name_list{kk});
    end
    for tt = 1:length(toggle_list)
        toggle = toggle_list{tt};
        R.(toggle).F1(ii,jj) =model_acc(M.index.eBIC).(toggle).F1;
        R.(toggle).MCC(ii,jj) =model_acc(M.index.eBIC).(toggle).MCC;
        R.(toggle).TPR(ii,jj) =model_acc(M.index.eBIC).(toggle).TPR;
        R.(toggle).FPR(ii,jj) =model_acc(M.index.eBIC).(toggle).FPR;
        R.(toggle).ACC(ii,jj) =model_acc(M.index.eBIC).(toggle).ACC;
    end
    R.bias(ii,jj) =model_acc(M.index.eBIC).bias;
    fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n', ...
        mean(R.common.F1(ii,1:jj)),mean(R.common.MCC(ii,1:jj)),mean(R.common.ACC(ii,1:jj)),mean(R.common.FPR(ii,1:jj)),mean(R.common.TPR(ii,1:jj)))
    
end
save([performance_path,'T150_CGN_result_K',int2str(K)],'R')
save([performance_path,'T150_CGN_ALL_RESULT_K',int2str(K)],'ALL_RESULT')
%%
clear R
clear ALL_RESULT
realization = 100;
for jj=1:realization
    fprintf('(%d,%d)\n',ii,jj)
    GTmodel = E{2,jj};
    fname = [resultpath,'estim_CGN_cvx_T150_',mname{ii},'percent_lag3_K',int2str(K),'_',int2str(jj)];
    load(fname)
%     M = augment_score_old(M,T);
    model_acc = performance_eval(M,GTmodel);
    toggle_list = {'total','common','differential'};
    %         M.index.bic=best_index(jj);
    ALL_RESULT(ii,jj).model_acc = model_acc;
    for kk=1:length(name_list)
        R.index(ii,jj).(name_list{kk}) = M.index.(name_list{kk});
    end
    for tt = 1:length(toggle_list)
        toggle = toggle_list{tt};
        R.(toggle).F1(ii,jj) =model_acc(M.index.eBIC).(toggle).F1;
        R.(toggle).MCC(ii,jj) =model_acc(M.index.eBIC).(toggle).MCC;
        R.(toggle).TPR(ii,jj) =model_acc(M.index.eBIC).(toggle).TPR;
        R.(toggle).FPR(ii,jj) =model_acc(M.index.eBIC).(toggle).FPR;
        R.(toggle).ACC(ii,jj) =model_acc(M.index.eBIC).(toggle).ACC;
    end
    R.bias(ii,jj) =model_acc(M.index.eBIC).bias;
    fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n', ...
        mean(R.common.F1(ii,1:jj)),mean(R.common.MCC(ii,1:jj)),mean(R.common.ACC(ii,1:jj)),mean(R.common.FPR(ii,1:jj)),mean(R.common.TPR(ii,1:jj)))
    
end
save([performance_path,'T150_CGN_CVX_result_K',int2str(K)],'R')
save([performance_path,'T150_CGN_CVX_ALL_RESULT_K',int2str(K)],'ALL_RESULT')