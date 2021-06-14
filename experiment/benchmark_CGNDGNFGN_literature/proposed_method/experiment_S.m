%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './experiment/model_parameters/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';
mkdir(outpath)
type = 3; %S type
cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]

p_true = 1;
K = 5;
n = 20; % time-series channels
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
T = 100;
p_est = 1;
[P,~] = offdiagJSS(n,p_est,K);
Dtmp = diffmat(n,p_est,K);
D = sparse(Dtmp*P);
[~,~,dd,m] = size(E);
realz = m;
GridSize = 30;
mname = {'1','5'};
parameter.cvx.varorder = p_est;
parameter.cvx.formulation = 'fgn'; % cgn, dgn, fgn
parameter.cvx.penalty_weight = 'LS'; % LS, uniform
parameter.cvx.GridSize = GridSize;
parameter.cvx.data_concat = 0;
parameter.cvx.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.cvx.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER.cvx = gen_alg_params(parameter.cvx.qnorm, parameter.cvx.formulation);
parameter.ncvx = parameter.cvx;
parameter.ncvx.qnorm = 'ncvx';
ALG_PARAMETER.ncvx = gen_alg_params(parameter.ncvx.qnorm, parameter.ncvx.formulation);
for jj=1:m
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = jointvargc_optimizedPD(y,parameter.ncvx,ALG_PARAMETER.ncvx);
%         save([outpath,'estim_FGN_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        M = jointvargc_optimizedPD(y,parameter.cvx,ALG_PARAMETER.cvx);
%         save([outpath,'estim_FGN_cvx_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
