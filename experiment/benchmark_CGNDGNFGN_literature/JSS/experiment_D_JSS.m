%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './experiment/model_parameters/';
% outpath = '../formulation_D_result/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
mkdir(outpath)

type = 2; %D type
% cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
cd = 3;
T = 100;
p_true = 1;
p_est = 1;
% K = 5;
K = 50;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = m;
GridSize = 30;
mname = {'1','5'};
parameter.varorder = p_est;
parameter.formulation = 'dgn'; % cgn, dgn, fgn
parameter.penalty_weight = 'uniform'; % LS, uniform
parameter.GridSize = GridSize;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER = gen_alg_params(parameter.qnorm, parameter.formulation);
ALG_PARAMETER.PRINT_RESULT = 0;
for jj=1:realz
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = jointvargc(y,parameter,ALG_PARAMETER);
        save([outpath,'estim_DGN_JSS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
