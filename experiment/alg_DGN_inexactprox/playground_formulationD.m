clear
clc
load('G:\Shared drives\MASTER_DRIVE\THESIS\THESIS_DATA_Sparse\DATA_BANK_5_S2.mat')
trial_ind = 5;
T = 300;

y = data.y(:,1:T,:,trial_ind);
[n,Num,K] = size(y);

p=2;
[P,~] = offdiagJSS(n,p,K);
PInd = P*((1:1:n^2*p*K)');

Dmat = diffmat(n,p,K);
D = sparse(Dmat*P);
GridSize = 30;
Lambda = logspace(-3,0,GridSize);
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
pp=2;qq=0.5;
objfun = @(x) 0.5*norm(gc*x-yc)^2 + a1*normpq(x(PInd),pp,qq,p)+a2*normpq(x(PInd),pp,qq,p*K);
PARAMETER_nmAPG_BB.proj_ind = PInd;
PARAMETER_nmAPG_BB.p = pp;
PARAMETER_nmAPG_BB.q = qq;
PARAMETER_nmAPG_BB.dim = [n,p,K];
PARAMETER_nmAPG_BB.delta = 1;
PARAMETER_nmAPG_BB.eta = 0.7;
PARAMETER_nmAPG_BB.rho = 0.9;
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

%%
[x_D_cvx, ~, ~,history_D_cvx] = spectral_ADMM(gc, yc,a1, a2,pp,1,PARAMETER_D,xLS);
if size(gc,1)>size(gc,2)
    x0 = xLS;
else
    x0 = x_D_cvx;
end
%%
t1 =  tic;
[x_D, L1x_D,L2x_D, history_D] = spectral_ADMM(gc, yc, a1, a2,pp,qq, PARAMETER_D,x0);
t_spectral = toc(t1);

t1 = tic;

[x_nmAPG,history_nmAPG]= nmAPG_BB_FormulationD(gc,yc,a1,a2,PARAMETER_nmAPG,x0);
t_nmAPG = toc(t1);

t1 = tic;
[x_nmAPG_BB,history_nmAPG_BB]= nmAPG_BB_FormulationD(gc,yc,a1,a2,PARAMETER_nmAPG_BB,x0);
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

figure(3)
subplot(1,3,1)
semilogy(abs(diff(history_D.objval)))
subplot(1,3,2)
semilogy(abs(diff(history_nmAPG_BB.objval)))
subplot(1,3,3)
semilogy(abs(diff(history_nmAPG.objval)))
%% perturb spectral ADMM to initialize nmAPG
x_perturb = x_D+randn(n^2*p*K,1)*0.0001/norm(x_D);
[x_nmAPG_perturb,history_nmAPG_perturb]= nmAPG_BB_FormulationD(gc,yc,a1,a2,PARAMETER_nmAPG_BB,x_perturb);

fprintf(' sADMM\t\t nmAPG_BB\t\t loss\n')
fprintf(' %.2f\t %.2f\n',objfun(x_D),objfun(x_nmAPG_perturb))
fprintf(' relative parameter difference\n')
fprintf(' %f\n',norm(x_D-x_nmAPG_perturb)/norm(x_D))
%% initial spectral ADMM using high init rho with x_nmAPG as initial point
PARAMETER_test = PARAMETER_D;
PARAMETER_test.IS_ADAPTIVE = 0;
PARAMETER_test.rho_init = 1e4;
[x_test, ~,~, history_test] = spectral_ADMM(gc, yc, a1, a2,pp,qq, PARAMETER_test,x_nmAPG_BB);
fprintf(' sADMM\t\t nmAPG_BB\t\t loss\n')
fprintf(' %.2f\t %.2f\n',objfun(x_test),objfun(x_nmAPG))
figure(4)
plot(history_test.objval)
%% initial spectral ADMM using high init rho with perturbed x_D as initial point
[x_test, ~,~, history_test] = spectral_ADMM(gc, yc, a1, a2,pp,qq, PARAMETER_test,x_perturb);
fprintf(' sADMM\t\t original\t\t loss\n')
fprintf(' %.2f\t %.2f\n',objfun(x_test),objfun(x_D))
figure(4)
plot(history_test.objval)