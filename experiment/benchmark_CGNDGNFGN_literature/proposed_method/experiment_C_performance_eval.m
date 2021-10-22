%% Experiment: Common GC estimation, performance evaluation for CGN
clear
clc
inpath = './experiment/model_parameters/';
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_C_magda/';
performance_path = './results2plot/';
mname = {'10','20'};
name_list = {'bic_lasso','bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
realization = 100;
load([inpath,'model_K5_p1'])
% R.F1 = zeros(length(mname),realization);
% R.MCC = zeros(length(mname),realization);
% R.TPR = zeros(length(mname),realization);
% R.FPR = zeros(length(mname),realization);
% R.ACC = zeros(length(mname),realization);
toggle_list = {'common'};
for ii=1:length(mname)
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{2,ii+2,2,jj};
        fname = [resultpath,'estim_CGN_',mname{ii},'percent_',int2str(jj)];
        load(fname)
        model_acc = performance_eval(M,GTmodel);
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
        
    end
end
save([performance_path,'CGN_result'],'R')
save([performance_path,'CGN_ALL_RESULT'],'ALL_RESULT')
%%
% clear
% clc
% inpath = './data_compare/';
% performance_path = './experiment/result_to_plot/';
% mname = {'10','20'};
% name_list = {'bic_lasso','bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
% load([performance_path,'adaptive_formulation_C_result'])
% load([performance_path,'adaptive_formulation_C_ALL_RESULT'])
% acc_list = {'TPR','FPR','ACC','F1','MCC'};
% realization = 100;
% ARR = zeros(5,length(name_list));
% for ii=1:length(mname)
%     for nn=1:length(name_list)
%         for kk=1:length(acc_list)
%             for jj=1:realization
%                 index_selected = R.index(ii,jj).(name_list{nn});
%                 
%                 summary.common.(acc_list{kk})(ii,jj) =ALL_RESULT(ii,jj).model_acc(index_selected).common.(acc_list{kk});
%             end
%             ARR(kk,nn) = mean(summary.common.(acc_list{kk})(ii,:));
%         end
% %         ARR(:,nn)= [mean(summary.common.F1(ii,:)),mean(summary.common.MCC(ii,:)),mean(summary.common.ACC(ii,:)),mean(summary.common.FPR(ii,:)),mean(summary.common.TPR(ii,:))]';
%     end
%     disp(['common density:',mname{ii}])
%     t = array2table(ARR,'VariableNames',name_list,'RowNames', acc_list);
%     t.Variables =  round(t.Variables*100,2);
%     disp(t)
% end
