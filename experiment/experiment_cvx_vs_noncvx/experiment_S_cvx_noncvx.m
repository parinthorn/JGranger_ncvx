%%
clear
clc
inpath = './data_compare/';
% outpath = '../formulation_D_result/';
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx\';
mkdir(outpath)

T = 150;
p_true = 3;
p_est = 3;
K = 5;
load([inpath,'compare_convex_model_K',int2str(K),'_p',int2str(p_true)]) % struct E
m= size(E,2);
realz = m;
GridSize = 30;
mname = {'1','5'};
%% Non-CVX
ii=2;
for jj=76:80
    % generate data from given seed
    model = E{3,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
    M = formulation_S(y,p_est,GridSize);
    
    save([outpath,'resultT150_adaptive_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
end
%% CVX
ii=2;
for jj=1:realz
    % generate data from given seed
    model = E{3,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
    M = test_cvxformulation_S(y,p_est,GridSize);
    save([outpath,'resultT150_cvx_adaptive_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

end
