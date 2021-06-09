%% Algorithm performance

clear
clc
inpath = './data_compare/';
% outpath = './results/';
% mkdir(outpath)
type = 2; %D type
cd = 3;
T = 100;
p_true = 1;
p_est = 2;
K = 5;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = m;
mname = {'1','5'};
model = E{type,cd,1,1};
y = sim_VAR(model.A,T,1,model.seed,0);
n = size(y,1);
parameter.varorder = p_est;
parameter.formulation = 'dgn'; % cgn, dgn, fgn
parameter.penalty_weight = 'LS'; % LS, uniform
parameter.GridSize = 30;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'

eff_T = T-p_est;
H = zeros(n*p_est,eff_T,K);
Y = zeros(n,eff_T,K);
disp('Generating H, Y matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p_est);
end
disp('vectorizing model')
[b,G] = vectorize_VAR(Y,H,[n,p_est,K,eff_T]);
xLS = G\b;
[P,~] = offdiagJSS(n,p_est,K);
Dtmp = diffmat(n,p_est,K);
D = Dtmp*P;
L1 = sparse(P);

if strcmp(parameter.formulation,'fgn')
    MODEL_DIM = [n, p_est,K,p_est,p_est];
    L2 = sparse(D);
else
    MODEL_DIM = [n, p_est,K,p_est,p_est*K];
    L2 = sparse(P);
end






[Lambdacrit_1.cvx,Lambdacrit_2.cvx] = gen_critical_lambdas(G,b, xLS,MODEL_DIM,'cvx',parameter.penalty_weight,parameter.formulation);
[Lambdacrit_1.ncvx,Lambdacrit_2.ncvx] = gen_critical_lambdas(G,b, xLS,MODEL_DIM,'ncvx',parameter.penalty_weight,parameter.formulation);
Lambda = logspace(-6,0,parameter.GridSize);
ii = 21;
jj = 23;
a1.cvx = Lambdacrit_1.cvx*Lambda(ii)*0;
a2.cvx = Lambdacrit_2.cvx*Lambda(jj);
a1.ncvx = Lambdacrit_1.ncvx*Lambda(ii)*0;
a2.ncvx = Lambdacrit_2.ncvx*Lambda(jj);

ALG_PARAMETER.cvx = gen_alg_params('cvx', parameter.formulation);
ALG_PARAMETER.cvx.PRINT_RESULT = 0;
ALG_PARAMETER.cvx.dim = MODEL_DIM;
ALG_PARAMETER.cvx.x0 = xLS;
[x.cvx,~,~,history.cvx] = spectral_ADMM(G, b,a1.cvx,a2.cvx,L1,L2,ALG_PARAMETER.cvx);
ALG_PARAMETER.ncvx = gen_alg_params('ncvx', parameter.formulation);
ALG_PARAMETER.ncvx.PRINT_RESULT = 0;
ALG_PARAMETER.ncvx.dim = MODEL_DIM;
ALG_PARAMETER.ncvx.x0 = xLS;

[x.ncvx,~,~,history.ncvx] = adaptive_ADMM(G, b,a1.ncvx,a2.ncvx,L1,L2,ALG_PARAMETER.ncvx);
%%
clf
figurepath = './plotting/figures/';
formulation = parameter.formulation;

% tt=tiledlayout(2,1);
figure(1);
% tmp = diff(history.cvx.objval);
semilogy(history.cvx.reldiff_norm,'linewidth',2)
grid on
xlabel('iterations')
h=ylabel('$\frac{\Vert x^{+}-x \Vert_{2}}{\Vert x \Vert_{2}}$','Interpreter','latex');
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',36)
set(h, 'FontSize', 60) 
set(gcf,'WindowState','fullscreen')
ylim([10^-7,2*sqrt(2)])
print([figurepath,'algperf_xreldiff_cvx_',formulation],'-dsvg','-r300')


figure(2);
semilogy(abs(((history.cvx.objval)-history.cvx.objval(end))/history.cvx.objval(end)),'linewidth',2)
grid on
xlabel('iterations')
h=ylabel('$\frac{F(x^{+})-p^{*}}{p^{*}}$','Interpreter','latex');
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',36)
set(h, 'FontSize', 60) 
set(gcf,'WindowState','fullscreen')
ylim([10^-11,2*sqrt(2)])
print([figurepath,'algperf_objreldiff_cvx_',formulation],'-dsvg','-r300')


figure(3);
% tmp = diff(history.ncvx.objval);
semilogy(history.ncvx.reldiff_norm,'linewidth',2)
grid on
xlabel('iterations')
h=ylabel('$\frac{\Vert x^{+}-x \Vert_{2}}{\Vert x \Vert_{2}}$','Interpreter','latex');
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',36)
set(h, 'FontSize', 60) 
set(gcf,'WindowState','fullscreen')
ylim([10^-7,2*sqrt(2)])
print([figurepath,'algperf_xreldiff_ncvx_',formulation],'-dsvg','-r300')


figure(4);
semilogy(abs(((history.ncvx.objval)-history.ncvx.objval(end))/history.ncvx.objval(end)),'linewidth',2)
grid on
xlabel('iterations')
h=ylabel('$\frac{F(x^{+})-p^{*}}{p^{*}}$','Interpreter','latex');
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',36)
set(h, 'FontSize', 60) 
set(gcf,'WindowState','fullscreen')
ylim([10^-8,2*sqrt(2)])
print([figurepath,'algperf_objreldiff_ncvx_',formulation],'-dsvg','-r300')

%% convex CGN comparison
PARAMETER.proj_ind = P*(1:1:n^2*p_est*K)';
PARAMETER.p=2;
PARAMETER.q=1;
PARAMETER.dim(1)=n;
PARAMETER.dim(2)=p_est;
PARAMETER.dim(3)=K;
PARAMETER.delta=0.7;
PARAMETER.eta=0.5;
PARAMETER.rho=0.5;
PARAMETER.IS_LINESEARCH = 1;
[x.cvx_nmAPG,history.cvx_nmAPG]= nmAPG_BB_cgn(G,b,a2.cvx,PARAMETER,xLS);
PARAMETER.q=0.5;
[x.ncvx_nmAPG,history.ncvx_nmAPG]= nmAPG_BB_cgn(G,b,a2.ncvx,PARAMETER,xLS);

PARAMETER.IS_LINESEARCH = 0;
PARAMETER.q=1;
[x.cvx_nmAPG_noBB,history.cvx_nmAPG_noBB]= nmAPG_BB_cgn(G,b,a2.cvx,PARAMETER,xLS);
PARAMETER.q=0.5;
[x.ncvx_nmAPG_noBB,history.ncvx_nmAPG_noBB]= nmAPG_BB_cgn(G,b,a2.ncvx,PARAMETER,xLS);

ALG_PARAMETER.cvx = gen_alg_params('cvx', 'cgn');
% ALG_PARAMETER.cvx.epscor = 0.01;
ALG_PARAMETER.cvx.PRINT_RESULT = 0;
ALG_PARAMETER.cvx.dim = MODEL_DIM;
ALG_PARAMETER.cvx.x0 = xLS;
ALG_PARAMETER.ncvx = gen_alg_params('ncvx', 'cgn');
ALG_PARAMETER.ncvx.PRINT_RESULT = 0;
ALG_PARAMETER.ncvx.dim = MODEL_DIM;
ALG_PARAMETER.ncvx.x0 = xLS;

[x.cvx,~,~,history.cvx] = spectral_ADMM(G, b,a1.cvx*0,a2.cvx,L1,L2,ALG_PARAMETER.cvx);

[x.ncvx,~,~,history.ncvx] = adaptive_ADMM(G, b,a1.ncvx*0,a2.ncvx,L1,L2,ALG_PARAMETER.ncvx);

ALG_PARAMETER.static = ALG_PARAMETER.cvx;
ALG_PARAMETER.static.IS_ADAPTIVE = 0;
ALG_PARAMETER.static.rho_init = 600;
[x.cvx_static,~,~,history.cvx_static] = spectral_ADMM(G, b,a1.cvx*0,a2.cvx,L1,L2,ALG_PARAMETER.static);
%%
figure(1)
tt=tiledlayout(2,1);
nexttile;
semilogy(history.cvx.reldiff_norm,'linewidth',3)
hold on
semilogy(history.cvx_nmAPG.reldiff_norm,'linewidth',3)
semilogy(history.cvx_nmAPG_noBB.reldiff_norm,'linewidth',3)
hold off
xlabel('iterations')
ylabel('$\frac{\Vert x^{+}-x \Vert_{2}}{\Vert x \Vert_{2}}$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Spectral ADMM', 'nmAPG with BB', 'nmAPG no BB','location','eastoutside')
nexttile;
semilogy(((history.cvx.objval)-history.cvx.objval(end))/history.cvx.objval(end),'linewidth',3)
grid on
hold on
semilogy(((history.cvx_nmAPG.objval)-history.cvx_nmAPG.objval(end))/history.cvx_nmAPG.objval(end),'linewidth',3)
semilogy(((history.cvx_nmAPG_noBB.objval)-history.cvx_nmAPG_noBB.objval(end))/history.cvx_nmAPG_noBB.objval(end),'linewidth',3)
grid on
hold off
xlabel('iterations')
ylabel('$\frac{f(x^{+})-p^{*}}{p^{*}}$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Spectral ADMM', 'nmAPG with BB', 'nmAPG no BB','location','eastoutside')

figure(2)
tt=tiledlayout(2,1);
nexttile;
semilogy(history.ncvx.reldiff_norm,'linewidth',3)
grid on
hold on
semilogy(history.ncvx_nmAPG.reldiff_norm,'linewidth',3)
semilogy(history.ncvx_nmAPG_noBB.reldiff_norm,'linewidth',3)
grid on
hold off
xlabel('iterations')
ylabel('$\frac{\Vert x^{+}-x \Vert_{2}}{\Vert x \Vert_{2}}$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Adaptive ADMM', 'nmAPG with BB', 'nmAPG no BB','location','eastoutside')
nexttile;
semilogy(((history.ncvx.objval)-history.ncvx.objval(end))/history.ncvx.objval(end),'linewidth',3)
grid on
hold on
semilogy(((history.ncvx_nmAPG.objval)-history.ncvx_nmAPG.objval(end))/history.ncvx_nmAPG.objval(end),'linewidth',3)
semilogy(((history.ncvx_nmAPG_noBB.objval)-history.ncvx_nmAPG_noBB.objval(end))/history.ncvx_nmAPG_noBB.objval(end),'linewidth',3)
grid on
hold off
xlabel('iterations')
ylabel('$\frac{f(x^{+})-p^{*}}{p^{*}}$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Adaptive ADMM', 'nmAPG with BB', 'nmAPG no BB','location','eastoutside')
%% time usage nonconvex

figure(3)
% plot([0,history.ncvx.t],history.ncvx.objval,'linewidth',3)
% hold on
% plot([0;history.ncvx_nmAPG.t],history.ncvx_nmAPG.objval,'linewidth',3)
% plot([0;history.ncvx_nmAPG_noBB.t],history.ncvx_nmAPG_noBB.objval,'linewidth',3)
% hold off
plot(history.ncvx.objval,'linewidth',3)
hold on
plot(history.ncvx_nmAPG.objval,'linewidth',3)
plot(history.ncvx_nmAPG_noBB.objval,'linewidth',3)
hold off
xlabel('iteration')
ylabel('objective')
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
grid on
legend('Adaptive ADMM', 'nmAPG with BB', 'nmAPG no BB')
print([figurepath,'algperf_CGN'],'-dsvg','-r300')
%% time usage convex

figure(4)
% h1=plot([0,history.cvx.t],history.cvx.objval,'linewidth',3);
% hold on
% h2=plot([0,history.cvx_static.t],history.cvx_static.objval,'linewidth',3);
% plot([0;history.cvx_nmAPG.t],history.cvx_nmAPG.objval,'linewidth',3)
% plot([0;history.cvx_nmAPG_noBB.t],history.cvx_nmAPG_noBB.objval,'linewidth',3)
% hold off

h1=plot(history.cvx.objval,'linewidth',3);
hold on
h2=plot(history.cvx_static.objval,'linewidth',3);
h3=plot(history.cvx_nmAPG.objval,'linewidth',3);
h4=plot(history.cvx_nmAPG_noBB.objval,'linewidth',3);
hold off

xlabel('iteration')
ylabel('objective')
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Spectral ADMM', 'fixed ADMM', 'nmAPG with BB', 'nmAPG no BB')
uistack(h4,'top')
uistack(h3,'top')
uistack(h2,'top')
uistack(h1,'top')

print([figurepath,'algperf_cvxCGN'],'-dsvg','-r300')