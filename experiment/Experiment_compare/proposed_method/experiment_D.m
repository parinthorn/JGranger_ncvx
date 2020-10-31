%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = '../formulation_D_result/';

type = 2; %D type
% cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
cd = 2;
T = 100;
p_true = 1;
p_est = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p_est,K);
% load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
load([inpath,'test_model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = 5;
GridSize = 30;
mname = {'1','5'};
for ii=1:dd
    for jj=1:realz
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = formulation_D(y,P,p_est,GridSize);
%         save([outpath,'result_formulationD_',mname{ii},'percent','_',int2str(jj),'_p',int2str(p)],'M')
    end
end
