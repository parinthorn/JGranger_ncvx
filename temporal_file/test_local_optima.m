%% This temporal file test if K=50 case suffered from poor local optimum
clear
clc
global x_cheat
inpath = './data_compare/';
outpath = '../formulation_D_result/';
mkdir(outpath)

type = 2; %D type
% cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
cd = 3;
T = 100;
p_true = 1;
p_est = 1;
K = 50;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p_est,K);
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
dd=1;
realz = 1;
GridSize = 30;
mname = {'1','5'};
idx = efficient_vect([n,p_est,K]);
for ii=1:dd
    for jj=1:realz
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        x_cheat = model.A(idx);
        M = formulation_D(y,P,p_est,GridSize);
        save([outpath,'cheat_result_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
%     result_formulationD_1percent_lag1_K5_12
    end
end