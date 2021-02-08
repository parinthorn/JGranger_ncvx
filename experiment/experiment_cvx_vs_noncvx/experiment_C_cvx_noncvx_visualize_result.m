clear
clc
clf;close all
K=5;
type_acc = {'total','common','differential'};
acc_list = {'ACC','F1','MCC'};
acc_list_2 = {'TPR','FPR','ACC','F1','MCC'};
name_list = {'bic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};%{'bic','aicc'};
resultpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx\';

path_to_ALL_RESULT = {[resultpath,'formulation_CT150_cvx_adaptive_ALL_RESULT_K',int2str(K)], ...
    [resultpath,'formulation_CT150_adaptive_ALL_RESULT_K',int2str(K)]};
path_to_index = {[resultpath,'formulation_CT150_cvx_adaptive_result_K',int2str(K)], ...
    [resultpath,'formulation_CT150_adaptive_result_K',int2str(K)]};

diff_den = {'1%','5%'};
cnt = 0 ;
TYPE_NAME = {'CVX','NON-CVX'};
subplot_index = [1,3,5,2,4,6];
ii=2;
%%
for tt=1:length(type_acc)
    cnt = cnt+1;
    figure(cnt);
    text_title = type_acc{tt};
    sgtitle([text_title,', formulation C convex(left), non-convex(right)'])
    subplot_cnt = 0;
    for type=1:2
        load(path_to_ALL_RESULT{type})
        load(path_to_index{type})
        %         dd = size(ALL_RESULT,1);
        realz = size(ALL_RESULT,2);
        for ss=1:length(acc_list)
            subplot_cnt = subplot_cnt+1;
            val = zeros(1,30);
            max_val = 0;
            for jj=1:realz
                tmp = [ALL_RESULT(ii,jj).model_acc.(type_acc{tt})];tmp = [tmp.(acc_list{ss})];val = val +reshape(tmp,1,30)/realz;
                max_val = max_val + max(tmp(:))/realz;
            end
            subplot(length(acc_list),2,subplot_index(subplot_cnt))
            plot(val,'LineWidth',2)
            title(sprintf('bestcase=%.3f',max_val))
            axis('square')
            ylim([0 1])
            set(gca,'xticklabel',[],'yticklabel',[])
            ylabel(acc_list{ss})
        end

        ARR = zeros(5,length(name_list));
        clear summary
        for nn=1:length(name_list)
            for kk=1:length(acc_list_2)
                for jj=1:realz
                    index_selected = R.index(ii,jj).(name_list{nn});
                    summary.total.(acc_list_2{kk})(ii,jj) =ALL_RESULT(ii,jj).model_acc(index_selected).total.(acc_list_2{kk});
                end
                ARR(kk,nn) = mean(summary.total.(acc_list_2{kk})(ii,:));
            end
        end
        disp(['density:',TYPE_NAME{type}])
        t = array2table(ARR,'VariableNames',name_list,'RowNames', acc_list_2);
        t.Variables =  round(t.Variables*100,2);
        disp(t)
    end
    
end

%% Box plot
clear ARR
ARR.cvx = zeros(length(acc_list_2),length(name_list),100);
ARR.ncvx = zeros(length(acc_list_2),length(name_list),44);
tp_name = {'cvx','ncvx'};
for type=1:2
    load(path_to_ALL_RESULT{type})
    load(path_to_index{type})
    realz = size(ALL_RESULT,2);

    clear summary
    for nn=1:length(name_list)
        for kk=1:length(acc_list_2)
            for jj=1:realz
                index_selected = R.index(ii,jj).(name_list{nn});
                ARR.(tp_name{type})(kk,nn,jj) = ALL_RESULT(ii,jj).model_acc(index_selected).common.(acc_list_2{kk});
            end
        end
    end    
end
%%
% acc_list_2 = {'TPR','FPR','ACC','F1','MCC'};
% name_list = {'bic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};%{'bic','aicc'};

for nn=3
%     tmp = squeeze(ARR(:,:,nn,:));
    figure(99)
%     sgtitle(name_list{nn})
    for kk=1:length(acc_list_2)
        subplot(1,length(acc_list_2),kk)
        
        tmp = [squeeze((ARR.cvx(kk,nn,:)));squeeze((ARR.ncvx(kk,nn,:)))];
        sz_1 = length(squeeze((ARR.cvx(kk,nn,:))));
        sz_2 = length(squeeze((ARR.ncvx(kk,nn,:))));
        g = [zeros(sz_1,1);ones(sz_2,1)];
        boxplot(100*tmp,g)
        title(acc_list_2{kk})
        set(gca,'xticklabel',TYPE_NAME,'fontsize',14)
        set(findobj(gca,'type','line'),'linew',2)
        if ~strcmp('FPR',acc_list_2{kk})
        ylim([0 100])
        else
            ylim([0 100])
        end
    end
end