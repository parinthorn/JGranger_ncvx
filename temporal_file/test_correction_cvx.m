%% CGN fix
clear
clc
inpath = './data_compare/';
% outpath = '../formulation_D_result/';
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx\';
mkdir(outpath)
T=150;
p_true = 3;
p_est = 3;
K = 5;
load([inpath,'compare_convex_model_K',int2str(K),'_p',int2str(p_true)]) % struct E
m= size(E,2);
realz = m;
GridSize = 30;
mname = {'1','5'};
ii=2;
for jj=1:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
%     M = noncvx_CGN(y,p_est,GridSize);
    load([outpath,'resultT150_adaptive_formulationC_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

    toggle = 'LLH_full';
    data_concat = 0;
    M = correction_C(y,M,p_est,GridSize,toggle,data_concat);
    save([outpath,'LLHcorrected_resultT150_adaptive_formulationC_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
end
for jj=1:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
%     M = cvx_CGN(y,p_est,GridSize);
    load([outpath,'resultT150_cvx_adaptive_formulationC_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

    toggle = 'LLH_full';
    data_concat = 0;
    M = correction_C(y,M,p_est,GridSize,toggle,data_concat);
    save([outpath,'LLHcorrected_resultT150_cvx_adaptive_formulationC_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

end

%% DGN fix
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
ii=2;
for jj=1:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
    %     M = noncvx_DGN(y,p_est,GridSize);
    load([outpath,'resultT150_adaptive_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    toggle = 'LLH_full';
    data_concat = 0;
    M = correction_S(y,M,p_est,GridSize,toggle,data_concat);
    save([outpath,'LLHcorrected_resultT150_adaptive_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
end
ii=2;
for jj=1:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
%     M = cvx_DGN(y,p_est,GridSize);
    load([outpath,'resultT150_cvx_adaptive_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

    toggle = 'LLH_full';
    data_concat = 0;
    M = correction_S(y,M,p_est,GridSize,toggle,data_concat);
    save([outpath,'LLHcorrected_resultT150_cvx_adaptive_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

end

%% FGN fix
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
ii=2;
for jj=1:realz    % generate data from given seed
    model = E{3,jj}; % type S
    y = sim_VAR(model.A,T,1,model.seed,0);
%     M = noncvx_FGN(y,p_est,GridSize);
    load([outpath,'resultT150_adaptive_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

        toggle = 'LLH_full';
    data_concat = 0;
    M = correction_S(y,M,p_est,GridSize,toggle,data_concat);
    save([outpath,'LLHcorrected_resultT150_adaptive_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
end
ii=2;
for jj=1:realz
    % generate data from given seed
    model = E{3,jj}; % type S
    y = sim_VAR(model.A,T,1,model.seed,0);
%     M = cvx_FGN(y,p_est,GridSize);
    load([outpath,'resultT150_cvx_adaptive_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
            toggle = 'LLH_full';
    data_concat = 0;
    M = correction_S(y,M,p_est,GridSize,toggle,data_concat);
    save([outpath,'LLHcorrected_resultT150_cvx_adaptive_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

end
