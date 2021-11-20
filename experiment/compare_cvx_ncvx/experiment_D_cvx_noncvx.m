%%
clear
clc
inpath = './experiment/model_parameters/';
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx\';
mkdir(outpath)

T = 150;
p_true = 3;
p_est = 3;
K = 5;
load([inpath,'compare_convex_model_K',int2str(K),'_p',int2str(p_true)]) % struct E
m= size(E,2);
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
%% Non-CVX
ii=2;
for jj=1:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
    M = jointvargc(y,parameter.ncvx,ALG_PARAMETER.ncvx);
    
    save([outpath,'estim_DGN_T150_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
end
%% CVX
ii=2;
for jj=1:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
    M = jointvargc(y,parameter.ncvx,ALG_PARAMETER.ncvx);
    save([outpath,'estim_DGN_cvx_T150_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

end
