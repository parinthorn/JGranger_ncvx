clear
clc
inpath = './data_compare/';


resultpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx\';
% resultpath = 'D:\parinthorn_thesis\formulation_D_result\';

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
realization = 97;
for jj=1:realization
    fprintf('(%d,%d)\n',ii,jj)
    GTmodel = E{2,jj};
    fname = [resultpath,'resultT150_adaptive_formulationD_',mname{ii},'percent_lag3_K',int2str(K),'_',int2str(jj)];
    load(fname)
    M = augment_score(M,T,'llh_full');
    model_acc = performance_eval(M,GTmodel);
    toggle_list = {'total','common','differential'};
    %         M.index.bic=best_index(jj);
    ALL_RESULT(ii,jj).model_acc = model_acc;
    for kk=1:length(name_list)
        R.index(ii,jj).(name_list{kk}) = M.index.(name_list{kk});
    end
    for tt = 1:length(toggle_list)
        toggle = toggle_list{tt};
        R.(toggle).F1(ii,jj) =model_acc(M.index.bic).(toggle).F1;
        R.(toggle).MCC(ii,jj) =model_acc(M.index.bic).(toggle).MCC;
        R.(toggle).TPR(ii,jj) =model_acc(M.index.bic).(toggle).TPR;
        R.(toggle).FPR(ii,jj) =model_acc(M.index.bic).(toggle).FPR;
        R.(toggle).ACC(ii,jj) =model_acc(M.index.bic).(toggle).ACC;
    end
    fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n', ...
        mean(R.total.F1(ii,1:jj)),mean(R.total.MCC(ii,1:jj)),mean(R.total.ACC(ii,1:jj)),mean(R.total.FPR(ii,1:jj)),mean(R.total.TPR(ii,1:jj)))
    
end
save([resultpath,'formulation_DT150_adaptive_result_K',int2str(K)],'R')
save([resultpath,'formulation_DT150_adaptive_ALL_RESULT_K',int2str(K)],'ALL_RESULT')
%%
realization = 100;
for jj=1:realization
    fprintf('(%d,%d)\n',ii,jj)
    GTmodel = E{2,jj};
    fname = [resultpath,'resultT150_cvx_adaptive_formulationD_',mname{ii},'percent_lag3_K',int2str(K),'_',int2str(jj)];
    load(fname)
    M = augment_score(M,T,'llh');
    model_acc = performance_eval(M,GTmodel);
    toggle_list = {'total','common','differential'};
    %         M.index.bic=best_index(jj);
    ALL_RESULT(ii,jj).model_acc = model_acc;
    for kk=1:length(name_list)
        R.index(ii,jj).(name_list{kk}) = M.index.(name_list{kk});
    end
    for tt = 1:length(toggle_list)
        toggle = toggle_list{tt};
        R.(toggle).F1(ii,jj) =model_acc(M.index.bic).(toggle).F1;
        R.(toggle).MCC(ii,jj) =model_acc(M.index.bic).(toggle).MCC;
        R.(toggle).TPR(ii,jj) =model_acc(M.index.bic).(toggle).TPR;
        R.(toggle).FPR(ii,jj) =model_acc(M.index.bic).(toggle).FPR;
        R.(toggle).ACC(ii,jj) =model_acc(M.index.bic).(toggle).ACC;
    end
    fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n', ...
        mean(R.total.F1(ii,1:jj)),mean(R.total.MCC(ii,1:jj)),mean(R.total.ACC(ii,1:jj)),mean(R.total.FPR(ii,1:jj)),mean(R.total.TPR(ii,1:jj)))
    
end
save([resultpath,'formulation_DT150_cvx_adaptive_result_K',int2str(K)],'R')
save([resultpath,'formulation_DT150_cvx_adaptive_ALL_RESULT_K',int2str(K)],'ALL_RESULT')