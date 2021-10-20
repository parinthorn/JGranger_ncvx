clear
clc

outpath = "G:\My Drive\0FROM_SHARED_DRIVE\THESIS\reviewer_response_results\";
inpath = './experiment/model_parameters/';
type = 2; %D type
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
parameter.cvx.formulation = 'dgn'; % cgn, dgn, fgn
parameter.cvx.penalty_weight = 'LS'; % LS, uniform
parameter.cvx.GridSize = GridSize;
parameter.cvx.data_concat = 0;
parameter.cvx.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.cvx.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER.cvx = gen_alg_params(parameter.cvx.qnorm, parameter.cvx.formulation);

model = E{type,3,2,1};
T = 150;
T_bootstrap = 100;
n=20;

num_blocks = 5;
repetitions = 5;

y = sim_VAR(model.A,T,1,model.seed,1);


for ii=0:T-T_bootstrap-1
    bootstrap_sample{ii+1} = (1:T_bootstrap) + ii;
%     bootstrap_index = randi([1 100],1,5);
end
s = RandStream('mlfg6331_64');

common_GC_bootstrap = zeros(n,n, repetitions);
for tt=1:repetitions

bootstrap_index = randsample(s, length(bootstrap_sample),num_blocks,false);
y_bootstrap = zeros(n, num_blocks*T_bootstrap+p_true*num_blocks, K);

cnt = 0;
for bb=1:length(bootstrap_index)
    cnt = cnt+1;
    y_bootstrap(:,p_true*cnt+(T_bootstrap*(cnt-1)+1:T_bootstrap*cnt),:) = y(:,bootstrap_sample{bootstrap_index(bb)},:);
end

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


%% visualize
n=20;
K=5;


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
imagesc(mean(common_GC_bootstrap,3))
hold on
spy(ground_truth_common_network,'r')
hold off
nexttile;
load('H:\Pycharm_projects\bootstrapping_commonGC_DGN_5runs.mat')
imagesc(mean(common_GC_bootstrap,3))
hold on
spy(ground_truth_common_network,'r')
hold off

colormap(1-gray)

