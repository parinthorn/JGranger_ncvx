%% Experiment: Common GC estimation, estimating VAR coefficients using DGN, convex DGN, K=5 and K=50.
clear
clc
inpath = './experiment/model_parameters/';
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
parameter.cvx.varorder = p_est;
parameter.cvx.formulation = 'dgn'; % cgn, dgn, fgn
parameter.cvx.penalty_weight = 'LS'; % LS, uniform
parameter.cvx.GridSize = GridSize;
parameter.cvx.data_concat = 0;
parameter.cvx.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.cvx.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER.cvx = gen_alg_params(parameter.cvx.qnorm, parameter.cvx.formulation);
parameter.ncvx = parameter.cvx;
parameter.ncvx.qnorm = 'ncvx';
ALG_PARAMETER.ncvx = gen_alg_params(parameter.ncvx.qnorm, parameter.ncvx.formulation);
for jj=1:realz
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = jointvargc_optimizedPD(y,parameter.ncvx,ALG_PARAMETER.ncvx);
%         save([outpath,'estim_DGN_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        M = jointvargc_optimizedPD(y,parameter.cvx,ALG_PARAMETER.cvx);
%         save([outpath,'estim_DGN_cvx_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
