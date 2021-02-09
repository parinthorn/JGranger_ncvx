clear
clc
clf;close all
K=5;
% SPEC
% ROW 1: ALL CRITERIA BOXPLOT (eBIC)
% ROW 2: Convex Grid
% ROW 3: Nonconvex Grid

type_acc = {'total','common','differential'};
acc_list = {'ACC','F1','MCC'};
acc_list_2 = {'TPR','FPR','ACC','F1','MCC'};
name_list = {'bic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};%{'bic','aicc'};
resultpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx\';

path_to_ALL_RESULT = {[resultpath,'formulation_DT150_cvx_adaptive_ALL_RESULT_K',int2str(K)], ...
    [resultpath,'formulation_DT150_adaptive_ALL_RESULT_K',int2str(K)]};
path_to_index = {[resultpath,'formulation_DT150_cvx_adaptive_result_K',int2str(K)], ...
    [resultpath,'formulation_DT150_adaptive_result_K',int2str(K)]};


diff_den = {'1%','5%'};
TYPE_NAME = {'CVX','NON-CVX'};
tp_name = {'cvx','ncvx'};
ii=2;
%%
for tp=1:2
    load(path_to_ALL_RESULT{tp})
    realz.(tp_name{tp}) = size(ALL_RESULT,2);
    tmp.(tp_name{tp}) = ALL_RESULT;
    load(path_to_index{tp})
    index.(tp_name{tp})=R.index;
    
end
clear ALL_RESULT
ALL_RESULT = tmp;
%%
clc
clf;
close all
selection = 'bic';
for type_index = 1:length(type_acc) % total common differential
    figure(type_index)
    cnt=0;
    for row_count=1:3
        
        for kk=1:length(acc_list_2) % TPR to MCC
            cnt=cnt+1;

            subplot(3,5,cnt)
            if row_count == 1
                for id_f = 1:length(tp_name)
                    for jj=1:realz.(tp_name{id_f})
                        index_selected = index.(tp_name{id_f})(ii,jj).(selection);
                        ARR.(tp_name{id_f})(kk,jj) = ALL_RESULT.(tp_name{id_f})(ii,jj).model_acc(index_selected).(type_acc{type_index}).(acc_list_2{kk});
                    end
                end
                
                tmp = [squeeze((ARR.cvx(kk,:))) squeeze((ARR.ncvx(kk,:)))];
                sz_1 = length(squeeze((ARR.cvx(kk,:))));
                sz_2 = length(squeeze((ARR.ncvx(kk,:))));
                g = [zeros(sz_1,1);ones(sz_2,1)];
                boxplot(100*tmp,g)
                title(acc_list_2{kk})
                set(gca,'xticklabel',TYPE_NAME,'fontsize',14)
                set(findobj(gca,'type','line'),'linew',2)
            else

                tp = row_count-1;
                model_selection = zeros(30,30);
                val = zeros(30,30);
                ay = zeros(realz.(tp_name{tp}),1);
                ax = zeros(realz.(tp_name{tp}),1);
                max_val = 0;
                for jj=1:realz.(tp_name{tp})
                    [ay(jj),ax(jj)] = ind2sub([30,30],index.(tp_name{tp})(ii,jj).(selection));
                    tmp = [ALL_RESULT.(tp_name{tp})(ii,jj).model_acc.(type_acc{type_index})];tmp = [tmp.(acc_list_2{kk})];val = val +reshape(tmp,30,30)/realz.(tp_name{tp});
                    max_val = max_val + max(tmp(:))/realz.(tp_name{tp});
                end
                imagesc(val)
                if ~((strcmp(acc_list_2{kk},'TPR')||strcmp(acc_list_2{kk},'FPR')))
                    title(sprintf('bestcase=%.3f',max_val))
                end
                
                axis('square')
                colormap((1-gray).^0.4)
                caxis([0,1])
                set(gca,'xticklabel',[],'yticklabel',[])
                ylabel(acc_list_2{kk})
                hold on
                scatter(ax,ay,'xr')
                hold off
            end
            
            
        end
        
    end
end




%% Box plot
clear ARR
ARR.cvx = zeros(length(acc_list_2),length(name_list),100);
ARR.ncvx = zeros(length(acc_list_2),length(name_list),44);

for type=1:2
    load(path_to_ALL_RESULT{type})
    load(path_to_index{type})
    realz = size(ALL_RESULT,2);
    
    clear summary
    for nn=1:length(name_list)
        for kk=1:length(acc_list_2)
            for jj=1:realz
                index_selected = R.index(ii,jj).(name_list{nn});
                ARR.(tp_name{type})(kk,nn,jj) = ALL_RESULT(ii,jj).model_acc(index_selected).total.(acc_list_2{kk});
            end
        end
    end
end
%%

for nn=3
    figure(99)
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
            ylim([0 10])
        end
        if kk==1
            ylabel('Percent')
        end
    end
end