%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = '../formulation_D_result/';

type = 2; %D type
cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
T = 100;
p = 3;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p,K);
load([inpath,'model_K',int2str(K)]) % struct E
[~,~,dd,m] = size(E);
GridSize = 30;
mname = {'1','5'};
for ii=1:dd
    for jj=1:m
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = formulation_D(y,P,p,GridSize);
        save([outpath,'result_formulationD_',mname{ii},'percent','_',int2str(jj),'_p',int2str(p)],'M')
    end
end
