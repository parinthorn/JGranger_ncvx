%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
load([inpath,'model_K5']) % struct E
[~,~,dd,m] = size(E);
type = 2; %D type
cd = 1; %common density set to 10%
T = 100;
p = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p,K);
GridSize = 30;
for ii=1:dd
    for jj=1:1
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = formulation_D(y,P,p,GridSize);
    end
end
