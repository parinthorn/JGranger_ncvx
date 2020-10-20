


%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = './experiment/Experiment_compare/skrip_code_2/data_R_formulationD/';
mkdir(outpath)


type = 2; %D type
cd = 1; %common density set to 10%
T = 100;
p = 1;
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
        y = reshape(permute(y,[1,3,2]),[n*K,T]);
        writematrix(y,[outpath,'data_',mname{ii},'percent_',int2str(jj),'.csv']) 
    end
end
