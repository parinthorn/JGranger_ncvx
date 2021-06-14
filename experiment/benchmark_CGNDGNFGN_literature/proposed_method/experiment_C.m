%% Experiment: Common GC estimation, estimating VAR coefficients using CGN, convex CGN.
clear
clc
inpath = './experiment/model_parameters/';
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\';

type = 2; %D type
T = 100;
p = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p,K);
load([inpath,'model_K',int2str(K),'_p1']) % struct E
[~,~,~,m] = size(E);
dd = 2;
GridSize = 30;
mname = {'10','20'};
cnt = 0;
parameter.cvx.varorder = p;
parameter.cvx.formulation = 'cgn'; % cgn, dgn, fgn
parameter.cvx.penalty_weight = 'LS'; % LS, uniform
parameter.cvx.GridSize = GridSize;
parameter.cvx.data_concat = 0;
parameter.cvx.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.cvx.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER.cvx = gen_alg_params(parameter.cvx.qnorm, parameter.cvx.formulation);
parameter.ncvx = parameter.cvx;
parameter.ncvx.qnorm = 'ncvx';
ALG_PARAMETER.ncvx = gen_alg_params(parameter.ncvx.qnorm, parameter.ncvx.formulation);
for ii=3:4
    cnt = cnt+1;
    for jj=1:m
        % generate data from given seed
        model = E{type,ii,dd,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        
        M = jointvargc_optimizedPD(y,parameter.ncvx,ALG_PARAMETER.ncvx);
        % save([outpath,'estim_CGN_',mname{cnt},'percent','_',int2str(jj)],'M')
        M = jointvargc_optimizedPD(y,parameter.cvx,ALG_PARAMETER.cvx);
%         save([outpath,'estim_CGN_cvx_',mname{cnt},'percent','_',int2str(jj)],'M')
    end
end
