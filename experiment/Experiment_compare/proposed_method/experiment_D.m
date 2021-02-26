%% Experiment: Common GC estimation, estimating VAR coefficients using DGN, convex DGN, K=5 and K=50.
clear
clc
inpath = './data_compare/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
mkdir(outpath)
type = 2; %D type
cd = 3;
T = 100;
p_true = 1;
p_est = 1;
K = 5;
% K = 50;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = m;
GridSize = 30;
mname = {'1','5'};
for jj=1:realz
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
%         M = formulation_D(y,p_est,GridSize);
%         save([outpath,'result_adaptive_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        M = test_cvxformulation_D(y,p_est,GridSize);
        save([outpath,'result_adaptive_cvx_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
