%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = './experiment/Experiment_compare/Wilms_code/data_R_formulationS/';
mkdir(outpath)


type = 3; %S type
cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
T = 100;
p = 1;
K = 5;
n = 20; % time-series channels
load([inpath,'model_K',int2str(K),'_p',int2str(p)]) % struct E
[~,~,dd,m] = size(E);
mname = {'1','5'};
dd=1;m=1;
for ii=1:dd
    for jj=1:m
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0); % dim [n,T,K]
        y = reshape(permute(y,[2,1,3]),[T,n*K]); % must be [T,n,K]
%         writematrix(y,[outpath,'K',int2str(K),'_data_',mname{ii},'percent_',int2str(jj),'.csv']) 
    end
end
