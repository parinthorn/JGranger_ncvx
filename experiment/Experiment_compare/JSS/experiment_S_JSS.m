%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
% outpath = '../formulation_S_result/';
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
% Dtmp = diffmat(n,p_est,K);
Dtmp = genDIFFmatJSS(n,p_est,K);
D = sparse(Dtmp*P);
[~,~,dd,m] = size(E);
realz = m;
GridSize = 30;
mname = {'1','5'};

parameter.varorder = p_est;
parameter.formulation = 'fgn'; % cgn, dgn, fgn
parameter.penalty_weight = 'uniform'; % LS, uniform
parameter.GridSize = GridSize;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER = gen_alg_params(parameter.qnorm, parameter.formulation);
ALG_PARAMETER.PRINT_RESULT = 0;
for jj=1:m
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = jointvargc(y,parameter,ALG_PARAMETER);
        save([outpath,'estim_FGN_JSS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
