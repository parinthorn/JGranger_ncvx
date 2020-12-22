clear
clc
inpath = './data_compare/';
outpath = '../formulation_D_result/';
mkdir(outpath)

type = 2; %D type
% cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
cd = 3;
T = 100;
p_true = 1;
p_est = 1;
% K = 50;
K_list = [5,15,25,35,50];

n = 20; % time-series channels

load([inpath,'model_K50_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = 10;
GridSize = 30;
mname = {'1','5'};
for ii=2:dd
    for jj=1:realz
        % generate data from given seed
        for K=K_list
            [P,~] = offdiagJSS(n,p_est,K);
            model = E{type,cd,ii,jj};
            y = sim_VAR(model.A(:,:,:,1:K),T,1,model.seed,0);
            M = formulation_D(y,P,p_est,GridSize);
            save([outpath,'varyK_result_formulationD_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        end
        %     result_formulationD_1percent_lag1_K5_12
    end
end
%% This experiment estimate VAR with formulation S by ADMM on type D groundtruth
clear
clc
inpath = './data_compare/';
outpath = '../formulation_D_result/';
mkdir(outpath)

type = 2; %D type
% cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
cd = 3;
T = 100;
p_true = 1;
p_est = 1;
% K = 50;
K_list = [5,15,25,35,50];

n = 20; % time-series channels

load([inpath,'model_K50_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = 10;
GridSize = 30;
mname = {'1','5'};
for ii=1:dd
    for jj=1:realz
        % generate data from given seed
        for K=K_list
            [P,~] = offdiagJSS(n,p_est,K);
            Dtmp = diffmat(n,p_est,K);
            D = sparse(Dtmp*P);
            model = E{type,cd,ii,jj};
            y = sim_VAR(model.A(:,:,:,1:K),T,1,model.seed,0);
            M = test_cvxformulation_S(y,P,D,p_est,GridSize);
            save([outpath,'varyK_result_cvxformulationS_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
        end
        %     result_formulationD_1percent_lag1_K5_12
    end
end