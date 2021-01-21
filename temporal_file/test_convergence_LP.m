%% This script test convergence of formulation S, D for a fixed value of rho
% to test adaptive_L and adaptive_P.
clear
clc
%% Generate data
inpath = './data_compare/';
outpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
type = 2; %D type
% cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
cd = 3;
T = 10;
p_true = 1;
p_est = 1;
K = 5;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
model =  E{type,cd,1,1};
y = sim_VAR(model.A,T,1,model.seed,0);
%% Initializing algorithm
qq=0.5;
[n,T,K] = size(y);
GridSize = 30;

p = p_est;
H = zeros(n*p,T-p,K);
Y = zeros(n,T-p,K);
disp('Generating H matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc] = vectorize_VAR(Y,H,[n,p,K,T]);
if size(gc,1)>size(gc,2)
    xLS = gc\yc;
else
    xLS = ridge(yc,gc,0.01);
end

%%
ALG_PARAMETER.PRINT_RESULT=1;

ALG_PARAMETER.dim = [n,p,K,p,p*K];

ALG_PARAMETER.epscor = 0.1;
ALG_PARAMETER.Ts = 100;
ALG_PARAMETER.is_chol = 1;
ALG_PARAMETER.multiplier = 2;
ALG_PARAMETER.toggle = 'formulationD';
ALG_PARAMETER.gamma = 1; % for adaptive case
x0 = xLS;
ii=20;
jj=20;
%% 

toggle = 'adaptive_L';
[Lambda_1,Lambda_2,opt] = grid_generation(gc,yc,GridSize,ALG_PARAMETER,qq,toggle);
ALG_PARAMETER.L1 = opt.L1;
ALG_PARAMETER.L2 = opt.L2;
ALG_PARAMETER.IS_ADAPTIVE =1;
ALG_PARAMETER.rho_init = 1;
[x_L, ~,~, history_L] = spectral_ADMM_adaptive(gc, yc, Lambda_1(:,ii), Lambda_2(:,jj),2,qq, ALG_PARAMETER,x0);
%%
toggle = 'adaptive_P';
[Lambda_1,Lambda_2,opt] = grid_generation(gc,yc,GridSize,ALG_PARAMETER,qq,toggle);
ALG_PARAMETER.L1 = opt.L1;
ALG_PARAMETER.L2 = opt.L2;
ALG_PARAMETER.IS_ADAPTIVE =1;
ALG_PARAMETER.rho_init = 1;
[x_P, ~,~, history_P] = spectral_ADMM_adaptive(gc, yc, Lambda_1(:,ii), Lambda_2(:,jj),2,qq, ALG_PARAMETER,x0);

%%
fprintf(' relative difference:%f \n L-P:%d \n P-L: %d \n objval diff:%f \n',norm(x_L-x_P)/norm(x_L), ...
length(setdiff(find(x_L),find(x_P))), ...
length(setdiff(find(x_P),find(x_L))), ...
(history_L.objval(end)-history_P.objval(end))/history_L.objval(end))
%% plot differences
