clear
clc
inpath = './data_compare/';


resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
performance_path = './experiment/result_to_plot/';
mname = {'1','5'};
dd = length(mname);
% K=5;
K = 50;
realization = 100;
load([inpath,'model_K',int2str(K),'_p1'])
name_list = {'bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
for jj=1:realization
    for ii=1:dd
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{2,3,ii,jj};
        fname = [resultpath,'result_JSS_formulationD_',mname{ii},'percent_lag1_K',int2str(K),'_',int2str(jj)];
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
            R.(toggle).F1(ii,jj) =model_acc(M.index.eBIC).(toggle).F1;
            R.(toggle).MCC(ii,jj) =model_acc(M.index.eBIC).(toggle).MCC;
            R.(toggle).TPR(ii,jj) =model_acc(M.index.eBIC).(toggle).TPR;
            R.(toggle).FPR(ii,jj) =model_acc(M.index.eBIC).(toggle).FPR;
            R.(toggle).ACC(ii,jj) =model_acc(M.index.eBIC).(toggle).ACC;
        end
        fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n', ...
            mean(R.total.F1(ii,1:jj)),mean(R.total.MCC(ii,1:jj)),mean(R.total.ACC(ii,1:jj)),mean(R.total.FPR(ii,1:jj)),mean(R.total.TPR(ii,1:jj)))
        
    end
end
save([performance_path,'adaptive_formulation_D_JSS_result_K',int2str(K)],'R')
save([performance_path,'adaptive_formulation_D_JSS_ALL_RESULT_K',int2str(K)],'ALL_RESULT')