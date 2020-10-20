%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = '../formulation_D_result/';

type = 2; %D type
cd = 1; %common density set to 10%
T = 1250;
p = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p,K);
load([inpath,'model_K',int2str(K)]) % struct E
[~,~,dd,m] = size(E);
GridSize = 30;
mname = {'1','5'};
for ii=1:dd
    for jj=1:m
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = formulation_D(y,P,p,GridSize);
        A_select = M.A(:,:,:,:,M.index.bic);
        
        % performance eval
        [M.stat.accuracy.common.confusion_matrix] = compare_sparsity(model.ind_common,M.ind_common,n,K,'commonROC');
        M.stat.accuracy.common.score = performance_score(M.stat.accuracy.common.confusion_matrix);
        [M.stat.accuracy.differential.confusion_matrix] = compare_sparsity(model.ind_differential,M.ind_differential,n,K,'differentialROC');
        M.stat.accuracy.differential.score = performance_score(M.stat.accuracy.differential.confusion_matrix);
        [M.stat.accuracy.selected_common] = compare_sparsity(model.ind_common,M.ind_common{M.index.bic},n,K,'single_common');
        [M.stat.accuracy.selected_differential] = compare_sparsity(model.ind_differential,M.ind_differential{M.index.bic},n,K,'single_differential');
        
        [M.stat.accuracy.total.confusion_matrix] = compare_sparsity(model.ind_nz,M.ind_nz,n,K,'differentialROC');
        M.stat.accuracy.total.score = performance_score(M.stat.accuracy.total.confusion_matrix);
        
        M.stat.accuracy.detail = {'TP','TN','FP','FN'};
        save([outpath,'result_formulationD_',int2str(ii),'percent','_',int2str(jj)])
    end
end
