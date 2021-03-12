%% Experiment: Common GC estimation, estimating VAR coefficients using CGN, convex CGN.
clear
clc
inpath = './data_compare/';
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\';

type = 2; %D type
T = 100;
p = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p,K);
load([inpath,'model_K',int2str(K),'_p1']) % struct E
[~,~,~,m] = size(E);
dd = 2;
GridSize = 30;
mname = {'10','20'};
cnt = 0;
for ii=3:4
    cnt = cnt+1;
    for jj=1:m
        % generate data from given seed
        model = E{type,ii,dd,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        
        
        
        % M = noncvx_CGN(y,p,GridSize);
        load([outpath,'result_JSS_formulationC_',mname{cnt},'percent','_',int2str(jj)],'M')
        toggle = 'LLH_full';
        data_concat = 0;
        M = correction_C(y,M,p,GridSize,toggle,data_concat);
        save([outpath,'LLHcorrected_result_JSS__formulationC_',mname{cnt},'percent','_',int2str(jj)],'M')
        
        
        
        % M = cvx_CGN(y,p,GridSize);
        load([outpath,'result_JSS_formulationC_',mname{cnt},'percent','_',int2str(jj)],'M')
        toggle = 'LLH_full';
        data_concat = 0;
        M = correction_C(y,M,p,GridSize,toggle,data_concat);
        save([outpath,'LLHcorrected_result_JSS_formulationC_',mname{cnt},'percent','_',int2str(jj)],'M')
    end
end




%% This experiment estimate VAR with FGN
clear
clc
inpath = './data_compare/';
% outpath = '../formulation_S_result/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';
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
for jj=1:m
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        
        
        %         M = noncvx_FGN(y,p_est,GridSize);
        load([outpath,'result_JSS_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        toggle = 'LLH_full';
        data_concat = 0;
        M = correction_S(y,M,p_est,GridSize,toggle,data_concat);
        save([outpath,'LLHcorrected_result_JSS_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        
        
        %         M = cvx_FGN(y,p_est,GridSize);
        load([outpath,'result_JSS_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        toggle = 'LLH_full';
        data_concat = 0;
        M = correction_S(y,M,p_est,GridSize,toggle,data_concat);
        save([outpath,'LLHcorrected_result_JSS_formulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
%% Experiment: GC estimation, estimating VAR coefficients using DGN, convex DGN, K=5 and K=50.
% clear
clc
inpath = './data_compare/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
mkdir(outpath)
type = 2; %D type
cd = 3;
T = 100;
p_true = 1;
p_est = 1;
K_list = [50];
% K = 5;
% K = 50;
for K=K_list
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = m;
GridSize = 30;
mname = {'1','5'};
for jj=64:realz
    for ii=1:dd
        % generate data from given seed
        model = E{type,cd,ii,jj};
        y = sim_VAR(model.A,T,1,model.seed,0);
        %         M = noncvx_DGN(y,p_est,GridSize);
        load([outpath,'result_JSS_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        toggle = 'LLH_full';
        data_concat = 0;
        M = correction_S(y,M,p_est,GridSize,toggle,data_concat);
        save([outpath,'LLHcorrected_result_JSS_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        %         M = cvx_DGN(y,p_est,GridSize);
        load([outpath,'result_JSS_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        toggle = 'LLH_full';
        data_concat = 0;
        M = correction_S(y,M,p_est,GridSize,toggle,data_concat);
        save([outpath,'LLHcorrected_result_JSS_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
    end
end
end
