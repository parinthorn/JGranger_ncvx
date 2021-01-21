clear
clc
inpath = './data_compare/';


resultpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_S_result\';
% resultpath = 'D:\parinthorn_thesis\formulation_D_result\';

mname = {'1','5'};
dd = length(mname);
% dd=1;
% dd=2;
realization = 100;
load([inpath,'model_K5_p1'])
name_list = {'bic_lasso','bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
for ii=1:dd
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{3,3,ii,jj};
        fname = [resultpath,'result_adaptive_formulationS_',mname{ii},'percent_lag1_K5_',int2str(jj)];
        load(fname)
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
        
    end
end
save([resultpath,'formulation_S_adaptive_result_K5'],'R')
save([resultpath,'formulation_S_adaptive_ALL_RESULT_K5'],'ALL_RESULT')