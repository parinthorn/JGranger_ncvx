%% This experiment estimate VAR with formulation S by Spectral ADMM
clear
clc
inpath = './data_compare/';
outpath = '../formulation_D_result/';
mkdir(outpath)
type = 2; %USE S FORMULATION ON D TYPE DATA
cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]

p_true = 1;
K = 50;
n = 20; % time-series channels
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
T = 100;
p_est = 1;
[P,~] = offdiagJSS(n,p_est,K);
Dtmp = diffmat(n,p_est,K);
D = sparse(Dtmp*P);
[~,~,dd,m] = size(E);
realz = 10;
GridSize = 30;
mname = {'1','5'};
for ii=1:dd
    for jj=1:realz
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = test_cvxformulation_S(y,P,D,p_est,GridSize);
       save([outpath,'cvx_result_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
