clear
clc
clf;close all
type_acc = {'total','common','differential'};
acc_list = {'ACC','F1','MCC'};
name_list = {'bic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};%{'bic','aicc'};
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
load([resultpath,'formulation_D_adaptive_ALL_RESULT_K50'])
load([resultpath,'formulation_D_adaptive_result_K50'])
dd = size(ALL_RESULT,1);
realz = size(ALL_RESULT,2);
diff_den = {'1%','5%'};
for ii=1:dd
    for tt=1:length(type_acc)
        figure;
        text_title = type_acc{tt};
        sgtitle([text_title,', diff density:',diff_den{ii}])
        for ss=1:length(acc_list)
            val = zeros(30,30);
            for jj=1:realz
                tmp = [ALL_RESULT(ii,jj).model_acc.(type_acc{tt})];tmp = [tmp.(acc_list{ss})];val = val +reshape(tmp,30,30)/realz;
            end
            subplot(length(acc_list),1,ss)
            imagesc(val)
            title(sprintf('bestcase=%.3f',max(max(val))))
            axis('square')
            colormap((1-gray).^0.4)
            caxis([0,1])
            set(gca,'xticklabel',[],'yticklabel',[])
            ylabel(acc_list{ss})
            
%             hold on
%             for mm=1:length(name_list)
%                 h = zeros(30,30);
%                 for jj=1:realz
% %                     h(index.(name_list{mm})(ii,jj)) = 1;
%                     h(R.index(ii,jj).(name_list{mm})) = 1;
%                 end
%                 [ind_row,ind_col] = ind2sub([30,30],find(h));
%                 scatter(ind_col,ind_row);
%             end
%             hold off

        end
    end
    ARR = zeros(5,length(name_list));
    for nn=1:length(name_list)
        
        for jj=1:realz
            index_selected = R.index(ii,jj).(name_list{nn});
            summary.total.F1(ii,jj) =ALL_RESULT(ii,jj).model_acc(index_selected).total.F1;
            summary.total.MCC(ii,jj) =ALL_RESULT(ii,jj).model_acc(index_selected).total.MCC;
            summary.total.TPR(ii,jj) =ALL_RESULT(ii,jj).model_acc(index_selected).total.TPR;
            summary.total.FPR(ii,jj) =ALL_RESULT(ii,jj).model_acc(index_selected).total.FPR;
            summary.total.ACC(ii,jj) =ALL_RESULT(ii,jj).model_acc(index_selected).total.ACC;
        end
%         disp(['Index:',name_list{nn}])
%         fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n \n\n\n', ...
%             mean(summary.total.F1(ii,:)),mean(summary.total.MCC(ii,:)),mean(summary.total.ACC(ii,:)),mean(summary.total.FPR(ii,:)),mean(summary.total.TPR(ii,:)))
        ARR(:,nn)= [mean(summary.total.F1(ii,:)),mean(summary.total.MCC(ii,:)),mean(summary.total.ACC(ii,:)),mean(summary.total.FPR(ii,:)),mean(summary.total.TPR(ii,:))]';
    end
%     load('C:\Users\CU_EE_LAB408\Dropbox\0MASTER\MATLAB_MASTER\JGranger_ncvx\experiment\Experiment_compare\skrip_code\data_R_formulationS\skrip_formulationS_accuracy_K5_55realz')
%     acc_list = {'F1','MCC','ACC','FPR','TPR'};
%     for kk=1:length(acc_list)
%         %     for nn=1:length(name_list)
%         ARR(kk,length(name_list)+1) = mean(score(ii).total.(acc_list{kk}));
%     end
%     
%     
    disp(['density:',diff_den{ii}])
    t = array2table(ARR,'VariableNames',name_list,'RowNames', {'F1', 'MCC', 'ACC', 'FPR','TPR'});
    t.Variables =  round(t.Variables*100,2);
    disp(t)
end
