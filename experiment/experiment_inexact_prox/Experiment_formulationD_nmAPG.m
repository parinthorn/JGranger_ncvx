clear
clc
load('G:\Shared drives\MASTER_DRIVE\THESIS\THESIS_DATA_Sparse\DATA_BANK_5_S2.mat')
trial_ind = 1;
y = data.y(:,:,:,trial_ind);
[n,Num,K] = size(y);

p=2;
[P,~] = offdiagJSS(n,p,K);
Dmat = diffmat(n,p,K);
D = sparse(Dmat*P);
GridSize = 30;
Lambda = logspace(-3,0,GridSize);
PARAMETER_nmAPG_BB.proj_ind = PInd;
PARAMETER_nmAPG_BB.p = pp;
PARAMETER_nmAPG_BB.q = qq;
PARAMETER_nmAPG_BB.dim = [n,p,K];
PARAMETER_nmAPG_BB.delta = 0.8;
PARAMETER_nmAPG_BB.eta = 0.1;
PARAMETER_nmAPG_BB.rho = 0.5;
PARAMETER_nmAPG_BB.IS_LINESEARCH = 1;

PARAMETER_nmAPG = PARAMETER_nmAPG_BB;
PARAMETER_nmAPG.IS_LINESEARCH = 0;

PARAMETER_D.PRINT_RESULT=1;
PARAMETER_D.IS_ADAPTIVE =1;
PARAMETER_D.L1 = P;
PARAMETER_D.L2 = P;
PARAMETER_D.dim = [n,p,K,p,p*K];
PARAMETER_D.rho_init = 1;
PARAMETER_D.epscor = 0.1;
PARAMETER_D.Ts = 100;

PARAMETER_S = PARAMETER_D;
PARAMETER_S.L2 = D;
PARAMETER_S.dim = [n,p,K,p,p];

tmp(K)= struct();
disp('Generating H matrix')
for kk=1:K
    [tmp(kk).H,tmp(kk).Y] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc,~] = veccoefmatgroup(reshape([tmp.Y],[n,Num-p,K]),reshape([tmp.H],[n*p,Num-p,K]),[],[n,p,K,Num-p]);
disp('calculating Lambda max')
Lmax = lambdamax_grouplasso(gc,yc,[n ,p ,K]);
Lambda = [0 Lambda*Lmax];
xLS = gc\yc;

a1 = Lambda(10);
a2 = Lambda(9);
PInd = P*((1:1:n^2*p*K)');
pp=2;qq=0.5;
objfun = @(x) 0.5*norm(gc*x-yc)^2 + a1*normpq(x(PInd),pp,qq,p)+a2*normpq(x(PInd),pp,qq,p*K);
%%
[x_S_cvx, ~, ~,history_S_cvx] = spectral_ADMM(gc, yc,a1, a2,pp,1,PARAMETER_S,xLS);
[x_D_cvx, ~, ~,history_D_cvx] = spectral_ADMM(gc, yc,a1, a2,pp,1,PARAMETER_D,xLS);

if size(gc,1)>size(gc,2)
    x0 = xLS;
else
    x0 = x_D_cvx;
end
%% high rho
[x_S, L1x_S,L2x_S, history_S] = spectral_ADMM(gc, yc, P,D,a1, a2,pp,qq, [n,p,K,p,p], 1,0.1,100,x0);
[x_D, L1x_D,L2x_D, history_D] = spectral_ADMM(gc, yc, P,P,a1, a2,pp,qq, [n,p,K,p,p*K], 1,0.1,100,x0);


[x_S2, L1x_S2,L2x_S2, history_S2] = spectral_ADMM(gc, yc, P,D,a1, a2,pp,qq, [n,p,K,p,p], 0.001,0.1,100,x0);
[x_D2, L1x_D2,L2x_D2, history_D2] = spectral_ADMM(gc, yc, P,P,a1, a2,pp,qq, [n,p,K,p,p*K], 0.001,0.1,100,x0);
fprintf(' highS\t highD\t lowS\t lowD\n')
fprintf(' %10.4f\t %10.4f\t %10.4f\t %10.4f\n',history_S.objval(end),history_D.objval(end),history_S2.objval(end),history_D2.objval(end))
%%
t1 =  tic;
[x_D, L1x_D,L2x_D, history_D] = spectral_ADMM(gc, yc, P,P,a1, a2,pp,qq, [n,p,K,p,p*K], 1,0.1,100,x0);
t_spectral = toc(t1);

t1 = tic;

[x_nmAPG,history_nmAPG]= nmAPG_BB_FormulationD(gc,yc,PARAMETER_nmAPG,x0);
t_nmAPG = toc(t1);

t1 = tic;
[x_nmAPG_BB,history_nmAPG_BB]= nmAPG_BB_FormulationD(gc,yc,PARAMETER_nmAPG_BB,x0);
t_nmAPG_BB = toc(t1);


fprintf(' sADMM\t\t nmAPG_BB\t\t nmAPG\n')
fprintf(' %.2f\t %.2f\t\t %.2f\n',objfun(x_D),objfun(x_nmAPG_BB),objfun(x_nmAPG))
%
figure(1)
plot(0:1:length(history_D.objval)-1,history_D.objval,'color',[0 0.5 0]);
hold on
plot(0:1:length(history_nmAPG_BB.objval)-1,history_nmAPG_BB.objval,'r');
plot(0:1:length(history_nmAPG.objval)-1,history_nmAPG.objval,'b');
hold off
legend('spectral ADMM','nmAPG + BB','nmAPG')
xlabel('iteration')
ylabel('loss')
%
figure(2)
plot(history_D.tpi*(0:1:length(history_D.objval)-1),history_D.objval,'color',[0 0.5 0]);
hold on
plot(mean(history_nmAPG_BB.tpi)*(0:1:length(history_nmAPG_BB.objval)-1),history_nmAPG_BB.objval,'r');
plot(mean(history_nmAPG.tpi)*(0:1:length(history_nmAPG.objval)-1),history_nmAPG.objval,'b');
hold off

legend('spectral ADMM','nmAPG + BB','nmAPG')
xlabel('computation time (s)')
ylabel('loss')
% [x_proxgrad,history3] = prox_grad_D(gc,yc,PInd,a1,a2,pp,qq,[n,p,K],x0);
%% GPU vs CPU
% clc
% A = sparse(gc'*gc);
% for ii =1:100
% A = sprand(A);
% gA = gpuArray(single(full(A)));
% tmp = randn(size(A,1),1);
% 
% gb = gpuArray(gA*tmp);
% b = A*tmp;
% 
% t1 = tic;gx = gA\gb;tg(ii) = toc(t1);
% t1 = tic;x = A\b;ts(ii) = toc(t1);
% end % Result: GPU is much slower than sparse array computation
%% Indexing vs sparse projection matrix multiplication
% PInd = full(P*((1:1:n^2*p*K)'));
% P = sparse(P);
% for ii=1:100000
%     x = randn(size(P,2),1);
%     t1 = tic;Px = P*x;tg(ii) = toc(t1);
%     t1 = tic;xPInd = x(PInd);ts(ii) = toc(t1);
% end Result: avoid multiplication
%% constrained Least square
% formulation S
Dtilde = D(L2x_S==0,:);
Ptilde = P(L1x_S==0,:);
nz_ind = find((x_S)~=0);
G = gc(:,nz_ind);
Dtilde = Dtilde(:,nz_ind);
x_constrLS_S = refit_fused(gc,yc,D,L2x_S,P,L1x_S);
% formulation D
z_ind = find((x_D)==0);
nz_ind = setdiff(1:(n^2)*p*K,z_ind);
G = gc;
b = yc;
G(:,z_ind) = [];
xhat = G\b;
x_hat_all = zeros(n^2*p*K,1);
x_hat_all(nz_ind) = xhat;
x_constrLS_D = refit(gc,yc,P,L1x_D);