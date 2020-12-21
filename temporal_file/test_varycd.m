%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
mkdir(outpath)

type = 2; %D type
T = 100;
p_true = 1;
p_est = 1;
K = 50;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p_est,K);
load([inpath,'vary_CD_model',int2str(K),'_p',int2str(p_true)]) % struct E
[~,cc,dd,m] = size(E);
realz = m;
GridSize = 30;
cname = {'2','4','5','6','8'};
dname = {'8','6','5','4','2'};
% for cd=1:cc
for ii=2:dd
    for jj=1:realz
        % generate data from given seed
        model = E{type,ii,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = test_cvxformulation_D(y,P,p_est,GridSize);
        save([outpath,'cvx_CDvary_result_formulationD_c',cname{ii},'d',dname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
%     result_formulationD_1percent_lag1_K5_12
    end
end
% end
