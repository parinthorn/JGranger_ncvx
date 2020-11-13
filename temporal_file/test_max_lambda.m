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
p_est = 3;
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
        for kk=1:K
            [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p_est);
        end
        disp('vectorizing model')
        [yc,gc] = vectorize_VAR(Y,H,[n,p_est,K,T]);
        disp('calculating Lambda max')
        Lmax(ii,jj) = lambdamax_grouplasso(gc,yc,[n ,p_est ,K]);
        
        Lmax_v2(ii,jj) = lambdamax_grouplasso_v2(gc,yc,p_est,[n ,p_est ,K]);
        M = test_lambdamax(y,P,p_est,GridSize);
    end
end

% the result suggested that Lambdamax from lambdamax_grouplasso_v2.m gives
% a tight bound for sparsest solution at lambda = Lmax_v2

