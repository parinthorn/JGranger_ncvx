%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\';

type = 2; %D type
% cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
T = 100;
p_true = 1;
p_est = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p_est,K);
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,~,m] = size(E);
dd = 2;
GridSize = 30;
mname = {'10','20'};
cnt = 0;
parameter.varorder = p_est;
parameter.formulation = 'cgn'; % cgn, dgn, fgn
parameter.penalty_weight = 'uniform'; % LS, uniform
parameter.GridSize = GridSize;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER = gen_alg_params(parameter.qnorm, parameter.formulation);
ALG_PARAMETER.PRINT_RESULT = 1;
for ii=3:4
    cnt = cnt+1;
    for jj=1:m
        % generate data from given seed
        model = E{type,ii,dd,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = jointvargc(y,parameter,ALG_PARAMETER);
        save([outpath,'estim_CGN_JSS_',mname{cnt},'percent','_',int2str(jj)],'M')
    end
end
