%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = '../formulation_S_result/';

type = 3; %S type
cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
T = 100;
p = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p,K);
Dtmp = diffmat(n,p,K);
D = sparse(D*P);
load([inpath,'model_K',int2str(K)]) % struct E
[~,~,dd,m] = size(E);
GridSize = 30;
mname = {'1','5'};
for ii=1:dd
    for jj=1:m
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = formulation_S(y,P,D,p,GridSize);
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
        M.stat.bias = squeeze(sqrt(sum(bsxfun(@minus, M.A,model.A).^2,[1,2,3,4]))./sqrt(sum(model.A.^2,'all'));
        M.stat.accuracy.detail = {'TP','TN','FP','FN'};
        save([outpath,'result_formulationS_',int2str(ii),'percent','_',int2str(jj)])
    end
end
