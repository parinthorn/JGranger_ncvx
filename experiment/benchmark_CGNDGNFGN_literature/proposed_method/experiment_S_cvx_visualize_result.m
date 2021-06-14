clear
clc
clf;close all
type_acc = {'total','common','differential'};
acc_list = {'ACC','F1','MCC'};
acc_list_2 = {'TPR','FPR','ACC','F1','MCC'};
name_list = {'bic_lasso','bic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};%{'bic','aicc'};
% name_list = {'bic','aic','aicc'};
performance_path = './results2plot/';
load([performance_path,'FGN_cvx_ALL_RESULT_K5'])
load([performance_path,'FGN_cvx_result_K5'])
% load([resultpath,'formulation_S_index'])
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
            max_val = 0;
            for jj=1:realz
                tmp = [ALL_RESULT(ii,jj).model_acc.(type_acc{tt})];tmp = [tmp.(acc_list{ss})];val = val +reshape(tmp,30,30)/realz;
                max_val = max_val + max((tmp(:)))/realz;
            end
            subplot(length(acc_list),1,ss)
            imagesc(val)
            title(sprintf('bestcase=%.3f',max_val))
            axis('square')
            colormap((1-gray).^0.4)
            caxis([0,1])
            set(gca,'xticklabel',[],'yticklabel',[])
            ylabel(acc_list{ss})
            
            hold on
            for mm=4
                h = zeros(30,30);
                for jj=1:realz
%                     h(index.(name_list{mm})(ii,jj)) = 1;
                    h(R.index(ii,jj).(name_list{mm})) = 1;
                end
                [ind_row,ind_col] = ind2sub([30,30],find(h));
                scatter(ind_col,ind_row);
            end
            hold off

        end
    end
    ARR = zeros(5,length(name_list));
     for nn=1:length(name_list)
        for kk=1:length(acc_list_2)
            for jj=1:realz
                index_selected = R.index(ii,jj).(name_list{nn});
                summary.total.(acc_list_2{kk})(ii,jj) =ALL_RESULT(ii,jj).model_acc(index_selected).total.(acc_list_2{kk});
            end
            ARR(kk,nn) = mean(summary.total.(acc_list_2{kk})(ii,:));
        end
     end
    
    disp(['density:',diff_den{ii}])
    t = array2table(ARR,'VariableNames',name_list,'RowNames', acc_list_2);
    t.Variables =  round(t.Variables*100,2);
    disp(t)
end

%% skrip

% name_list = {'bic_lasso','bic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};%{'bic','aicc'};

for ii=1:2
    ARR = zeros(length(acc_list),1);

%     end
    disp(['density:',diff_den{ii}])
    t = array2table(ARR,'VariableNames',{'Skrip'},'RowNames', acc_list);
    t.Variables =  round(t.Variables*100,2);
    disp(t)
end
