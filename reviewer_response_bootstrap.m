clear
clc

outpath = "G:\My Drive\0FROM_SHARED_DRIVE\THESIS\reviewer_response_results\";
inpath = './experiment/model_parameters/';
type = 1; %D type
cd = 3;

p_true = 1;
p_est = 1;
K = 5;
% K = 50;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = m;
GridSize = 30;
mname = {'1','5'};
parameter.cvx.varorder = p_est;
parameter.cvx.formulation = 'cgn'; % cgn, dgn, fgn
parameter.cvx.penalty_weight = 'LS'; % LS, uniform
parameter.cvx.GridSize = GridSize;
parameter.cvx.data_concat = 0;
parameter.cvx.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.cvx.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER.cvx = gen_alg_params(parameter.cvx.qnorm, parameter.cvx.formulation);

model = E{type,3,2,55};
T = 150;
T_bootstrap = 75;
n=20;

num_blocks = 2;
repetitions = 100;

y = sim_VAR(model.A,T,1,model.seed,1);


for ii=0:T-T_bootstrap-1
    bootstrap_sample{ii+1} = (1:T_bootstrap) + ii;
%     bootstrap_index = randi([1 100],1,5);
end
s = RandStream('mlfg6331_64');

common_GC_bootstrap = zeros(n,n, repetitions);
for tt=1:repetitions

bootstrap_index = randsample(s, length(bootstrap_sample),num_blocks,false);
% y_bootstrap = zeros(n, num_blocks*T_bootstrap+p_true*num_blocks, K);

cnt = 0;
eff_T = (T_bootstrap-p_true)*num_blocks;
H = zeros(n*p_true,eff_T,K);
Y = zeros(n,eff_T,K);
for bb=1:length(bootstrap_index)
    cnt = cnt+1;
    y_bootstrap(:,p_true*cnt+(T_bootstrap*(cnt-1)+1:T_bootstrap*cnt),:) = y(:,bootstrap_sample{bootstrap_index(bb)},:);
    
    

    disp('Concatenating H, Y matrix')
    for kk=1:K
    [H(:,(T_bootstrap-p_true)*(bb-1)+1:(T_bootstrap-p_true)*(bb),kk),Y(:,(T_bootstrap-p_true)*(bb-1)+1:(T_bootstrap-p_true)*(bb),kk)] = H_gen(y(:,bootstrap_sample{bootstrap_index(bb)},kk),p_true);
    end
end
parameter.cvx.is_YH = 1;
parameter.cvx.Y= Y;
parameter.cvx.H=H;
parameter.cvx.eff_T=eff_T;

M_bootstrap = jointvargc(y_bootstrap,parameter.cvx,ALG_PARAMETER.cvx);
result = (M_bootstrap.model(M_bootstrap.index.eBIC).GC~=0);
GC = ones(n,n)-eye(n);
diag_ind = 1:n+1:n^2;
for kk=1:K
    tmp = result(:,:,kk);
    tmp(diag_ind) = 0;
    GC = GC & (tmp~=0);
end


common_GC_bootstrap(:,:,tt) = GC;
end
% save('H:\Pycharm_projects\bootstrapping_commonGC_DGN_100runs.mat', 'common_GC_bootstrap')

%% visualize
n=20;
K=5;
figurepath = './results2plot/figures/';

ground_truth_common_network = ones(n,n)-eye(n);
diag_ind = 1:n+1:n^2;
for kk=1:K
    tmp = model.GC(:,:,kk);
    tmp(diag_ind) = 0;
    ground_truth_common_network = ground_truth_common_network & (tmp~=0);
end


figure(1)
tt =tiledlayout(1,2);
nexttile;
load('H:\Pycharm_projects\bootstrapping_commonGC_CGN_100runs.mat')
imagesc(sum(common_GC_bootstrap,3))
hold on
spy(ground_truth_common_network, 'r',48)
hold off
title("CGN")


nexttile;
load('H:\Pycharm_projects\bootstrapping_commonGC_DGN_100runs.mat')
imagesc(sum(common_GC_bootstrap,3))
hold on
spy(ground_truth_common_network, 'r',48)
hold off
title("DGN")

sgtitle("Detected GC network in 100 bootstrap samples", 'FontSize', 28)

colormap((1-gray).^(0.7))
colorbar


set(findall(gcf,'-property','FontSize'),'FontSize',18)
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
print([figurepath,'reviewer_response_bootstrap'],'-dpng','-r300')
print([figurepath,'reviewer_response_bootstrap'],'-depsc','-r300')