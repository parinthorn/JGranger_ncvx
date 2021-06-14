%% This demo program show the process of, data generation and estimation
clear
clc

%% demo data preparation
inpath = './experiment/model_parameters/';
type = 2; %D type
cd = 3;
T = 100;
p_true = 1;
K = 5;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)])
[~,~,dd,m] = size(E);
realz = m;
mname = {'1','5'};
model = E{type,cd,1,1};

%% The time-series 'y' has dimension [n,T,K], [n channel x T time-points x K sets]
y = sim_VAR(model.A,T,1,model.seed,0); % time-series generation, or can be the array with same structure
parameter.varorder = 1;
parameter.formulation = 'cgn'; % cgn, dgn, fgn
parameter.penalty_weight = 'LS'; % LS, uniform
parameter.GridSize = 30;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER = gen_alg_params(parameter.qnorm, parameter.formulation);
ALG_PARAMETER.PRINT_RESULT = 0;

M = jointvargc(y,parameter,ALG_PARAMETER);
plot_group_GC(M.model(M.index.eBIC).GC)

switch parameter.formulation
    case 'cgn'
        tmp = [M.model.stat]; tmp = [tmp.model_selection_score];eBIC = reshape([tmp.eBIC],parameter.GridSize,1); figure; plot(eBIC);
    case 'dgn'
        tmp = [M.model.stat]; tmp = [tmp.model_selection_score];eBIC = reshape([tmp.eBIC],parameter.GridSize,parameter.GridSize); imagesc(eBIC)
    case 'fgn'
        tmp = [M.model.stat]; tmp = [tmp.model_selection_score];eBIC = reshape([tmp.eBIC],parameter.GridSize,parameter.GridSize); imagesc(eBIC)
end