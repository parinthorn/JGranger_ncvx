%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = '../formulation_S_result/';

type = 3; %S type
cd = 2; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
T = 100;
p = 3;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p,K);
Dtmp = diffmat(n,p,K);
D = sparse(Dtmp*P);
load([inpath,'model_K',int2str(K),'_p',int2str(p)]) % struct E
[~,~,dd,m] = size(E);
realz = 10;
GridSize = 30;
mname = {'1','5'};
for ii=1:dd
    for jj=1:realz
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = formulation_S(y,P,D,p,GridSize);
        
        save([outpath,'result_formulationS_',int2str(ii),'percent','_',int2str(jj)])
    end
end
