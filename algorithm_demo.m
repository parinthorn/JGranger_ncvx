%% This script is to run an example on the algorithm comparison
% The problem is to solve min_{x} (1/2)||Gx-b||_2^2 + a_1h(Px) + a_2h(Px) [Formulation DGN]
% using inexact nmAPG algorithm
% The description of h_1, h_2 is provided in the manuscript

%% Data generation
clc
inpath = './experiment/model_parameters/';
figurepath = './results2plot/figures/';
% outpath = './results2plot/';
% mkdir(outpath)
type = 2; %D type
cd = 3;
T = 100;
p_true = 1;

K = 5;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = m;
mname = {'1','5'};
model = E{type,cd,1,1};
y = sim_VAR(model.A,T,1,model.seed,0);
n = size(y,1);
%% data processing & config estimation parameter

p_est = 2;

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
% 
L2 = sparse(P);
MODEL_DIM = [n, p_est,K,p_est,p_est*K];

% L2 = sparse(D);
% MODEL_DIM = [n, p_est,K,p_est,p_est];
%% ADMM vs Spectral ADMM
% Setting 
% We compare Spectral ADMM with fixed penalty ADMM
% fixed penalty (rho = 0.1*rho_S, 10*rho_S, 10*max(a1,a2))
parameter.varorder = p_est;
parameter.formulation = 'dgn'; % cgn, dgn
parameter.penalty_weight = 'LS'; % LS, uniform
parameter.GridSize = 30;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'

[Lambdacrit_1.cvx,Lambdacrit_2.cvx] = gen_critical_lambdas(G,b, xLS,MODEL_DIM,'cvx',parameter.penalty_weight,parameter.formulation);
Lambda = logspace(-6,0,parameter.GridSize);
ii = 21;
jj = 23;
a1.cvx = Lambdacrit_1.cvx*Lambda(ii);
a2.cvx = Lambdacrit_2.cvx*Lambda(jj);





ALG_PARAMETER.cvx = gen_alg_params('cvx', parameter.formulation);
ALG_PARAMETER.cvx.PRINT_RESULT = 0;
ALG_PARAMETER.cvx.dim = MODEL_DIM;
ALG_PARAMETER.cvx.x0 = xLS;
[x.adaptive_rho,~,~,history.adaptive_rho] = spectral_ADMM(G, b,a1.cvx,a2.cvx,L1,L2,ALG_PARAMETER.cvx);

ALG_PARAMETER.low_rho = ALG_PARAMETER.cvx;
ALG_PARAMETER.low_rho.IS_ADAPTIVE = 0;
ALG_PARAMETER.low_rho.rho_init = history.adaptive_rho.rho(end)*0.1;
[x.low_rho,~,~,history.low_rho] = spectral_ADMM(G, b,a1.cvx,a2.cvx,L1,L2,ALG_PARAMETER.low_rho);

ALG_PARAMETER.equal_rho = ALG_PARAMETER.cvx;
ALG_PARAMETER.equal_rho.IS_ADAPTIVE = 0;
ALG_PARAMETER.equal_rho.rho_init = history.adaptive_rho.rho(end);
[x.equal_rho,~,~,history.equal_rho] = spectral_ADMM(G, b,a1.cvx,a2.cvx,L1,L2,ALG_PARAMETER.equal_rho);

ALG_PARAMETER.high_rho = ALG_PARAMETER.cvx;
ALG_PARAMETER.high_rho.IS_ADAPTIVE = 0;
ALG_PARAMETER.high_rho.rho_init = history.adaptive_rho.rho(end)*10;
[x.high_rho,~,~,history.high_rho] = spectral_ADMM(G, b,a1.cvx,a2.cvx,L1,L2,ALG_PARAMETER.high_rho);

ALG_PARAMETER.rule_rho = ALG_PARAMETER.cvx;
ALG_PARAMETER.rule_rho.IS_ADAPTIVE = 0;
ALG_PARAMETER.rule_rho.rho_init = full(10*max([a1.cvx;a2.cvx],[],'all'));
[x.rule_rho,~,~,history.rule_rho] = spectral_ADMM(G, b,a1.cvx,a2.cvx,L1,L2,ALG_PARAMETER.rule_rho);

%% plot objective
type_list = {'adaptive_rho', 'low_rho', 'equal_rho', 'high_rho', 'rule_rho'};
name_list = {'spectral', sprintf('rho=%.2f',ALG_PARAMETER.low_rho.rho_init),  ...
    sprintf('rho=%.2f',ALG_PARAMETER.equal_rho.rho_init), ...
    sprintf('rho=%.2f',ALG_PARAMETER.high_rho.rho_init), ...
    sprintf('rho=%.2f',ALG_PARAMETER.rule_rho.rho_init)};
name_list = {'spectral rule', ... 
    ['0.1 $\rho_s$ ($\rho$=', sprintf('%.2f)',ALG_PARAMETER.low_rho.rho_init)], ...
    ['1 $\rho_s$ ($\rho$=', sprintf('%.2f)',ALG_PARAMETER.equal_rho.rho_init)], ...
    ['10 $\rho_s$ ($\rho$=', sprintf('%.2f)',ALG_PARAMETER.high_rho.rho_init)], ...
    ['(Songsiri,15) ($\rho$=', sprintf('%.2f)',ALG_PARAMETER.rule_rho.rho_init)]};

% tt = tiledlayout(5,1);
for ii=1:length(type_list)
    figure(ii);
    % yyaxis left
    semilogy(abs(((history.(type_list{ii}).objval)-history.(type_list{ii}).objval(end))/history.(type_list{ii}).objval(end)),'linewidth',2)
    yyaxis right
    plot(history.(type_list{ii}).rho,'linewidth',2)
    ylabel('$\rho$','Interpreter','latex')
    yyaxis left
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
%     h=ylabel(name_list{ii},'Interpreter','none');
    h=ylabel('$\frac{F(x_{k})-p^{*}}{p^{*}}$','Interpreter','latex');
    set(gca,'fontsize',60)
    set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
    print([figurepath,'ADMM_comparison_objdiff_',type_list{ii}],'-depsc','-r300')
end
% linkaxes(ax,'x')
close all
for ii=1:length(type_list)
    figure(ii);
    % yyaxis left
    semilogy(history.(type_list{ii}).reldiff_norm,'linewidth',2)
%     yyaxis right
%     plot(history.(type_list{ii}).rho,'linewidth',2)
%     ylabel('$\rho$','Interpreter','latex')
%     yyaxis left
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
%     h=ylabel(name_list{ii},'Interpreter','none');
h=ylabel('$\frac{\Vert x_{k+1}-x_{k} \Vert_{2}}{\Vert x_{k} \Vert_{2}}$','Interpreter','latex');
    set(gca,'fontsize',60)
    set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
    print([figurepath,'ADMM_comparison_reldiff_',type_list{ii}],'-depsc','-r300')
end

figure(ii+1);
for ii=1:length(type_list)
%     plot(history.(type_list{ii}).objval,'linewidth',2)
semilogy(abs(((history.(type_list{ii}).objval)-history.(type_list{ii}).objval(end))/history.(type_list{ii}).objval(end)),'linewidth',2)
    if ii==1
        hold on
    elseif ii==length(type_list)
        hold off
    end
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
    ylabel('$\frac{F(x_{k})-p^{*}}{p^{*}}$','Interpreter','latex');
%     h=ylabel(type_list{ii},'Interpreter','none');
    set(gca,'fontsize',60)
%     set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
end
legend(name_list,'Interpreter','latex','fontsize',32)
% set(gca,'Interpreter','latex')
print([figurepath,'ADMM_comparison_relative_obj'],'-depsc','-r300')


%% Non-convex CGN (Adaptive ADMM vs nmAPG)



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
PARAMETER.q=0.5;




parameter.varorder = p_est;
parameter.formulation = 'cgn'; % cgn, dgn
parameter.penalty_weight = 'LS'; % LS, uniform
parameter.GridSize = 30;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'
% for nonconvex
[Lambdacrit_1.ncvx,Lambdacrit_2.ncvx] = gen_critical_lambdas(G,b, xLS,MODEL_DIM,'ncvx',parameter.penalty_weight,parameter.formulation);
Lambda = logspace(-6,0,parameter.GridSize);
ii = 21;
jj = 23;
a1.ncvx = Lambdacrit_1.ncvx*Lambda(ii);
a2.ncvx = Lambdacrit_2.ncvx*Lambda(jj);

if strcmp(parameter.formulation,'cgn')
    a1.cvx = a1.cvx*0;
end

ALG_PARAMETER = gen_alg_params('ncvx', parameter.formulation);
ALG_PARAMETER.dim = MODEL_DIM;
ALG_PARAMETER.x0 = xLS;
[x.ADMM,~,~,history.ADMM] = adaptive_ADMM(G, b,a1.ncvx*0,a2.ncvx,L1,L2,ALG_PARAMETER);
[x.nmAPG,history.nmAPG]= nmAPG_BB_cgn(G,b,a2.ncvx,PARAMETER,xLS);
%% plot 
close all
type_list = {'ADMM', 'nmAPG'};
name_list = {sprintf('Adaptive ADMM (obj.=%.2f)',history.ADMM.objval(end)),sprintf('nmAPG (obj.=%.2f)',history.nmAPG.objval(end))};
figure(1)
for ii=1:length(type_list)
%     plot(history.(type_list{ii}).objval,'linewidth',2)
semilogy(abs(((history.(type_list{ii}).objval)-history.(type_list{ii}).objval(end))/history.(type_list{ii}).objval(end)),'linewidth',2)
    if ii==1
        hold on
    elseif ii==length(type_list)
        hold off
    end
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
    h=ylabel('$\frac{F(x_{k})-p^{*}}{p^{*}}$','Interpreter','latex');
    set(gca,'fontsize',60)
    set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
end
legend(name_list,'Interpreter','latex', 'FontSize', 36)
print([figurepath,'ADMM_comparison_nmAPG_CGN_objdiff'],'-depsc','-r300')
figure(2)
for ii=1:length(type_list)
%     plot(history.(type_list{ii}).objval,'linewidth',2)
    semilogy(history.(type_list{ii}).reldiff_norm,'linewidth',2)
    if ii==1
        hold on
    elseif ii==length(type_list)
        hold off
    end
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
    h=ylabel('$\frac{\Vert x_{k+1}-x_{k} \Vert_{2}}{\Vert x_{k} \Vert_{2}}$','Interpreter','latex');
    set(gca,'fontsize',60)
    set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
end
legend(name_list,'Interpreter','latex', 'FontSize', 36)
print([figurepath,'ADMM_comparison_nmAPG_CGN_reldiff'],'-depsc','-r300')

figure(3)
for ii=1:length(type_list)
%     plot(history.(type_list{ii}).objval,'linewidth',2)
plot(history.(type_list{ii}).objval,'linewidth',2)
    if ii==1
        hold on
    elseif ii==length(type_list)
        hold off
    end
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
        h=ylabel('$F(x_{k})$','Interpreter','latex');
    set(gca,'fontsize',60)
    set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
end
legend(name_list,'Interpreter','latex', 'FontSize', 36)
print([figurepath,'ADMM_comparison_nmAPG_CGN_obj'],'-depsc','-r300')

%% Non-convex DGN (Adaptive ADMM vs inexact nmAPG)
% We reduce K = 5 to K = 3 to show

K_new = 5;


[b_trunc,G_trunc] = vectorize_VAR(Y(:,:,1:K_new),H(:,:,1:K_new),[n,p_est,K_new,eff_T]);
xLS_trunc = G_trunc\b_trunc;

MODEL_DIM = [n, p_est,K_new,p_est,p_est*K_new];

parameter.varorder = p_est;
parameter.formulation = 'dgn'; % cgn, dgn
parameter.penalty_weight = 'LS'; % LS, uniform
parameter.GridSize = 30;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'

[P,~] = offdiagJSS(n,p_est,K_new);
Dtmp = diffmat(n,p_est,K_new);
D = Dtmp*P;
L1 = sparse(P);
L2 = sparse(P);

% for nonconvex
[Lambdacrit_1.ncvx,Lambdacrit_2.ncvx] = gen_critical_lambdas(G_trunc,b_trunc, xLS_trunc,MODEL_DIM,'ncvx',parameter.penalty_weight,parameter.formulation);
Lambda = logspace(-6,0,parameter.GridSize);
ii = 21;
jj = 23;
a1.ncvx = full(Lambdacrit_1.ncvx*Lambda(ii));
a2.ncvx = full(Lambdacrit_2.ncvx*Lambda(jj));

ALG_PARAMETER = gen_alg_params('ncvx', parameter.formulation);
ALG_PARAMETER.dim = MODEL_DIM;
ALG_PARAMETER.x0 = xLS_trunc;
[x.ADMM,~,~,history.ADMM] = adaptive_ADMM(G_trunc, b_trunc,a1.ncvx,a2.ncvx,L1,L2,ALG_PARAMETER);

PInd = P*((1:1:n^2*p_est*K_new)');


GridSize = 30;
PARAMETER_nmAPG_BB.proj_ind = PInd;
PARAMETER_nmAPG_BB.p = 2;
PARAMETER_nmAPG_BB.q = 0.5;
PARAMETER_nmAPG_BB.dim = [n,p_est,K_new];
PARAMETER_nmAPG_BB.delta = 0.8;
PARAMETER_nmAPG_BB.eta = 0.1;
PARAMETER_nmAPG_BB.rho = 0.5;
PARAMETER_nmAPG_BB.IS_LINESEARCH = 1;
x0 = xLS_trunc;

% [x.inexact_nmAPG,history.inexact_nmAPG]= nmAPG_BB_FormulationD(G_trunc,b_trunc,a1.ncvx,a2.ncvx,PARAMETER_nmAPG_BB,x0);
[x.inexact_nmAPG,history.inexact_nmAPG]= nmAPG_BB_DGN(G_trunc,b_trunc,a1.ncvx,a2.ncvx,PARAMETER_nmAPG_BB,x0);
%% validate when constant [these two functions are necessary to generate same sequence] [VALIDATED]
% b1 = a1.ncvx;
% b2 = a2.ncvx;
% b1(:) = max(b1,[],'all');
% b2(:) = max(b2,[],'all');
% [x.a,history.a]= nmAPG_BB_FormulationD(G_trunc,b_trunc,b1(1),b2(1),PARAMETER_nmAPG_BB,x0);
% [x.b,history.b]= nmAPG_BB_DGN(G_trunc,b_trunc,b1,b2,PARAMETER_nmAPG_BB,x0);
%%

close all
type_list = {'ADMM', 'inexact_nmAPG'};
name_list = {sprintf('Adaptive ADMM (obj.=%.2f)',history.ADMM.objval(end)),sprintf('inexact nmAPG (obj.=%.2f)',history.inexact_nmAPG.objval(end))};
figure(1)
for ii=1:length(type_list)
%     plot(history.(type_list{ii}).objval,'linewidth',2)
semilogy(abs(((history.(type_list{ii}).objval)-history.(type_list{ii}).objval(end))/history.(type_list{ii}).objval(end)),'linewidth',2)
    if ii==1
        hold on
    elseif ii==length(type_list)
        hold off
    end
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
    h=ylabel('$\frac{F(x_{k})-p^{*}}{p^{*}}$','Interpreter','latex');
%     h=ylabel(type_list{ii},'Interpreter','none');
    set(gca,'fontsize',60)
    set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
end
legend(name_list,'Interpreter','latex','location','southeast', 'FontSize', 36)
print([figurepath,'ADMM_comparison_inexact_DGN_objdiff'],'-depsc','-r300')
figure(2)
for ii=1:length(type_list)
%     plot(history.(type_list{ii}).objval,'linewidth',2)
semilogy(history.(type_list{ii}).reldiff_norm,'linewidth',2)
    if ii==1
        hold on
    elseif ii==length(type_list)
        hold off
    end
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
h=ylabel('$\frac{\Vert x_{k+1}-x_{k} \Vert_{2}}{\Vert x_{k} \Vert_{2}}$','Interpreter','latex');
%     h=ylabel(type_list{ii},'Interpreter','none');
    set(gca,'fontsize',60)
    set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
end
legend(name_list,'Interpreter','latex','location','south', 'FontSize', 36)
print([figurepath,'ADMM_comparison_inexact_DGN_reldiff'],'-depsc','-r300')
figure(3)
for ii=1:length(type_list)
%     plot(history.(type_list{ii}).objval,'linewidth',2)
plot(history.(type_list{ii}).objval,'linewidth',2)
    if ii==1
        hold on
    elseif ii==length(type_list)
        hold off
    end
    grid on
    xlabel('iterations $(k)$','Interpreter','latex')
    h=ylabel('$F(x_{k})$','Interpreter','latex');
%     h=ylabel(type_list{ii},'Interpreter','none');
set(gca,'fontsize',60)
    set(h, 'FontSize', 60)
    set(gcf,'WindowState','fullscreen')
end
legend(name_list,'Interpreter','latex', 'FontSize', 36)
print([figurepath,'ADMM_comparison_inexact_DGN_compare'],'-depsc','-r300')
%% old version

[x.inexact_nmAPG_from_ADMM,history.inexact_nmAPG_from_ADMM]= nmAPG_BB_DGN(G_trunc,b_trunc,a1.ncvx,a2.ncvx,PARAMETER_nmAPG_BB,x.ADMM);
ALG_PARAMETER_tmp = ALG_PARAMETER;
ALG_PARAMETER_tmp.x0 = x.inexact_nmAPG;
ALG_PARAMETER_tmp.IS_ADAPTIVE = 0;
ALG_PARAMETER_tmp.rho_init = history.ADMM.rho(end);
[x.ADMM_from_inexact_nmAPG,~,~,history.ADMM_from_inexact_nmAPG] = adaptive_ADMM(G_trunc, b_trunc,a1.ncvx,a2.ncvx,L1,L2,ALG_PARAMETER_tmp);
