%% This experiment estimate VAR with formulation D by ADMM
clear
clc
clf
close all
inpath = './data_compare/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
mkdir(outpath)

type = 2; %D type

cd = 3;
T = 30;
p_true = 1;
p_est = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p_est,K);
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = 1;
r_list = [1:53];
GridSize = 30;
mname = {'1','5'};
type_acc = {'total','common','differential'};
acc_list = {'ACC','F1','MCC'};
name_list = {'bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
type_name = {'cvx','ncvx'};
for test_itr=1:length(r_list)
    jj= r_list(test_itr);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
%         M_cvx = test_cvxformulation_D(y,p_est,GridSize);
        
%         M_ncvx = formulation_D(y,p_est,GridSize);
%         clc
        load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\resultT30_adaptive_cvx_formulationD_',mname{ii},'percent_lag1_K5_',int2str(jj)])
        [M_cvx] = augment_score(M,T,'sse');
        load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\resultT30_adaptive_formulationD_',mname{ii},'percent_lag1_K5_',int2str(jj)])
        [M_ncvx] = augment_score(M,T,'sse');
        
        ALL_RESULT.cvx(ii,jj).model_acc = performance_eval(M_cvx,model);
        ALL_RESULT.ncvx(ii,jj).model_acc = performance_eval(M_ncvx,model);
        for kk=1:length(name_list)
            R.cvx.index(ii,jj).(name_list{kk}) = M_cvx.index.(name_list{kk});
            R.ncvx.index(ii,jj).(name_list{kk}) = M_ncvx.index.(name_list{kk});
        end
        %%% CVX
        ARR = zeros(3,length(name_list));
        for nn=1:length(name_list)
            for kk=1:length(acc_list)
                index_selected = R.cvx.index(ii,jj).(name_list{nn});
                summary_1.total.(acc_list{kk})(nn,ii,jj) =ALL_RESULT.cvx(ii,jj).model_acc(index_selected).total.(acc_list{kk});
                ARR(kk,nn) = mean(summary_1.total.(acc_list{kk})(nn,ii,r_list(1:test_itr)));
            end
        end
        disp(['CVX, density:',mname{ii},'%'])
        t = array2table(ARR,'VariableNames',name_list,'RowNames', acc_list);
        t.Variables =  round(t.Variables*100,2);
        disp(t)
        % NCVX
        ARR = zeros(3,length(name_list));
        for nn=1:length(name_list)
            for kk=1:length(acc_list)
                index_selected = R.ncvx.index(ii,jj).(name_list{nn});
                summary_2.total.(acc_list{kk})(nn,ii,jj) =ALL_RESULT.ncvx(ii,jj).model_acc(index_selected).total.(acc_list{kk});
                ARR(kk,nn) = mean(summary_2.total.(acc_list{kk})(nn,ii,r_list(1:test_itr)));
            end
        end
        disp(['NCVX, density:',mname{ii},'%'])
        t = array2table(ARR,'VariableNames',name_list,'RowNames', acc_list);
        t.Variables =  round(t.Variables*100,2);
        disp(t)
        figure(ii);
        sgtitle(['total, diff density:',mname{ii},'%'])
        subplot_cnt = 0;
        for fm=1:length(type_name)
            for ss=1:length(acc_list)
                subplot_cnt = subplot_cnt+1;
                val = zeros(30,30);
                max_val = 0;
                for sample=1:test_itr
                    tmp = [ALL_RESULT.(type_name{fm})(ii,r_list(sample)).model_acc.total];tmp = [tmp.(acc_list{ss})];val = val +reshape(tmp,GridSize,GridSize)/test_itr;
                    max_val = max_val + max(tmp(:))/test_itr;
                end
                subplot(2,length(acc_list),subplot_cnt)
                imagesc(val)
                title(sprintf('bestcase=%.3f',max_val))
                axis('square')
                colormap((1-gray).^0.4)
                caxis([0,1])
                set(gca,'xticklabel',[],'yticklabel',[])
                ylabel(acc_list{ss})
                pause(0.1)
            end
        end
    end
end
