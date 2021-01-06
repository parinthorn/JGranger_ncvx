%% This experiment estimate VAR with formulation D by ADMM
clear
clc
inpath = './data_compare/';
outpath = '../formulation_S_result/';
mkdir(outpath)
type = 3; %S type
cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]

p_true = 1;
K = 5;
n = 20; % time-series channels
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
T = 100;
p_est = 1;
[P,~] = offdiagJSS(n,p_est,K);
Dtmp = diffmat(n,p_est,K);
D = sparse(Dtmp*P);
[~,~,dd,m] = size(E);
realz = m;
GridSize = 30;
mname = {'1','5'};
for ii=1:1
    for jj=1:1
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        M = formulation_S(y,P,D,p_est,GridSize);
%        save([outpath,'result_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
