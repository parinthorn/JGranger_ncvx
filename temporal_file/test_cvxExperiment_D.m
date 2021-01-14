%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
mkdir(outpath)

type = 2; %D type
% cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
cd = 3;
T = 1000;
p_true = 1;
p_est = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p_est,K);
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = 10;
GridSize = 30;
mname = {'1','5'};
for jj=1:realz
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = test_cvxformulation_D(y,p_est,GridSize);
%         save([outpath,'result_fixAdaptive_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        %     result_formulationD_1percent_lag1_K5_12
    end
end
