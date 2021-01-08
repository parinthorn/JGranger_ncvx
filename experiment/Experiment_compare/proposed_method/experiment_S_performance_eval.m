clear
clc
inpath = './data_compare/';
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';
mname = {'1','5'};
realization = 62;
load([inpath,'model_K5_p1'])
name_list = {'bic','aicc'};
T=100;
p=1;
n= 20;
K=5;
for ii=1:length(mname)
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{3,3,ii,jj};
        fname = [resultpath,'result_fixdf_formulationS_',mname{ii},'percent_lag1_K5_',int2str(jj)];
        load(fname)
        model_acc = performance_eval(M,GTmodel);
        toggle_list = {'total','common','differential'};
        ALL_RESULT(ii,jj).model_acc = model_acc;
        
        tmp = [M.model];tmp=[tmp.stat];tmp=[tmp.model_selection_score];
        L = reshape([tmp.L],M.GridSize,M.GridSize);
        df = zeros(M.GridSize,M.GridSize);
        for a1=1:M.GridSize
            for a2 = 1:M.GridSize
                df(a1,a2) = length(find(M.model(a1,a2).A(:)));
            end
        end
        bic_lasso = L+log(T-p).*df;
        
        gamma = log(n^2*p*K)/log(n*(T-p));
        kappa = 1.5*(1-1/(2*gamma));
%         binom_term = arrayfun(@(x) nchoosek(n^2*p*K,x), df);
        binom_term = arrayfun(@(x) log_stirling_approx(n^2*p*K)-log_stirling_approx(n^2*p*K-x)-log_stirling_approx(x) , df);
        eBIC = L+log(T-p).*df + 2*kappa.*binom_term;
        
        GIC_2 =  L+df.*(n^2*p*K)^(1/3);
        GIC_3 =  L+df.*(2*log(n^2*p*K));
        GIC_4 =  L+df.*(2*(log(n^2*p*K)+log(log(n^2*p*K))));
        GIC_5 = L+df.*log(log(T-p)).*log(n^2*p*K);
        GIC_6 = L+df.*log(T-p).*log(n^2*p*K);
        % https://jmlr.org/papers/volume13/kim12a/kim12a.pdf
        
        [~,R.index(ii,jj).bic_lasso] = min(bic_lasso(:));
        [~,R.index(ii,jj).GIC_2] = min(GIC_2(:));
        [~,R.index(ii,jj).GIC_3] = min(GIC_3(:));
        [~,R.index(ii,jj).GIC_4] = min(GIC_4(:));
        [~,R.index(ii,jj).GIC_5] = min(GIC_5(:));
        [~,R.index(ii,jj).GIC_6] = min(GIC_6(:));
        [~,R.index(ii,jj).eBIC] = min(eBIC(:));

        index_selected = R.index(ii,jj).GIC_5;

        for kk=1:length(name_list)
            R.index(ii,jj).(name_list{kk}) = M.index.(name_list{kk});
        end
        
        for tt = 1:length(toggle_list)
            toggle = toggle_list{tt};
            R.(toggle).F1(ii,jj) =model_acc(index_selected).(toggle).F1;
            R.(toggle).MCC(ii,jj) =model_acc(index_selected).(toggle).MCC;
            R.(toggle).TPR(ii,jj) =model_acc(index_selected).(toggle).TPR;
            R.(toggle).FPR(ii,jj) =model_acc(index_selected).(toggle).FPR;
            R.(toggle).ACC(ii,jj) =model_acc(index_selected).(toggle).ACC;
        end
        fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n', ...
            mean(R.total.F1(ii,1:jj)),mean(R.total.MCC(ii,1:jj)),mean(R.total.ACC(ii,1:jj)),mean(R.total.FPR(ii,1:jj)),mean(R.total.TPR(ii,1:jj)))
        
    end
end
save([resultpath,'formulation_S_result_fixdf'],'R')
save([resultpath,'formulation_S_ALL_RESULT_fixdf'],'ALL_RESULT')

%% summarize model selection score
clear
clc
inpath = './data_compare/';
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';
mname = {'1','5'};
realization = 100;
load([inpath,'model_K5_p1'])
for ii=1:length(mname)
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{3,3,ii,jj};
        fname = [resultpath,'result_formulationS_',mname{ii},'percent_lag1_K5_',int2str(jj)];
        load(fname)
        tmp = [M.model];tmp=[tmp.stat];tmp=reshape([tmp.model_selection_score],[M.GridSize,M.GridSize]);
        DBG(ii,jj).model_selection_score =tmp;
    end
end
% save([resultpath,'formulation_S_summary_model_selection_score.mat'],'DBG')
%% check divergence
clear
clc
inpath = './data_compare/';
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';
mname = {'1','5'};
realization = 100;
load([inpath,'model_K5_p1'])
for ii=1:length(mname)
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{3,3,ii,jj};
        fname = [resultpath,'result_formulationS_',mname{ii},'percent_lag1_K5_',int2str(jj)];
        load(fname)
%         tmp = [M.model];tmp=[tmp.stat];tmp=reshape([tmp.model_selection_score],[M.GridSize,M.GridSize]);
        DBG(ii,jj).err = sum(M.flag,'all');
        if DBG(ii,jj).err>0
            fprintf('(%d,%d):DIVERGE\n',ii,jj)
        end
    end
end
