%% This script is to address reviewers' comments


%% No. 2 Visualize moderate networks

% Generate model

% Visualize model

%% No. 12 Calculate SNR using SNR = max(A)/sigma
clear
clc
inpath = './experiment/model_parameters/';
sigma = 1;
K = 5;
m = 100;

SNR.exp_CGN_K5 = zeros(2,100);
SNR.exp_DGN_K50 = zeros(2,100);
SNR.exp_FGN_K5 = zeros(2,100);
load([inpath,'model_K5','_p1']) % struct E
n=20;
p=1;
K=5;
sigmasq = 1;
arr1 = [];
arr2 = [];
arr3 = [];

for ii=3:4
    for jj=1:m
        % generate data from given seed
        model = E{2,ii,2,jj};
        for kk =1:K
            Ak = model.A(:,:,:,kk);
            A = [reshape(Ak,n,n*p) ; [eye(n*(p-1)) zeros(n*(p-1),n)]];
            B = [eye(n); zeros(n*(p-1),n)]; C = [eye(n) zeros(n,n*(p-1))];
            Se = sigmasq*eye(n); Sx = dlyap(A,B*Se*B'); Sy = C*Sx*C'; % sigmasq is the noise variance you generated in VAR system
            SNR_new = 10*log10( trace(Sy)/trace(Se));
            disp(SNR_new)
            arr1 = [arr1 SNR_new];
        end
        SNR.exp_CGN_K5(ii-2,jj) = max(max(abs(model.A(:))));
    end
end
for jj=1:100
    for ii=1:2
        % generate data from given seed
        model = E{3,3,ii,jj};
        for kk =1:K
            Ak = model.A(:,:,:,kk);
            A = [reshape(Ak,n,n*p) ; [eye(n*(p-1)) zeros(n*(p-1),n)]];
            B = [eye(n); zeros(n*(p-1),n)]; C = [eye(n) zeros(n,n*(p-1))];
            Se = sigmasq*eye(n); Sx = dlyap(A,B*Se*B'); Sy = C*Sx*C'; % sigmasq is the noise variance you generated in VAR system
            SNR_new = 10*log10( trace(Sy)/trace(Se));
            disp(SNR_new)
            arr2 = [arr2 SNR_new];
        end
        SNR.exp_FGN_K5(ii,jj) = max(max(abs(model.A(:))));
    end
end


load([inpath,'model_K50','_p1']) % struct E
for jj=1:100
    for ii=1:2
        % generate data from given seed
        model = E{2,3,ii,jj};
        for kk =1:K
            Ak = model.A(:,:,:,kk);
            A = [reshape(Ak,n,n*p) ; [eye(n*(p-1)) zeros(n*(p-1),n)]];
            B = [eye(n); zeros(n*(p-1),n)]; C = [eye(n) zeros(n,n*(p-1))];
            Se = sigmasq*eye(n); Sx = dlyap(A,B*Se*B'); Sy = C*Sx*C'; % sigmasq is the noise variance you generated in VAR system
            SNR_new = 10*log10( trace(Sy)/trace(Se));
            disp(SNR_new)
            arr3 = [arr3 SNR_new];
        end
        SNR.exp_DGN_K50(ii,jj) = max(max(abs(model.A(:))));
    end
end

x =  [SNR.exp_CGN_K5(:) SNR.exp_FGN_K5(:) SNR.exp_DGN_K50(:)];


histogram(x, 'Normalization','probability')
%%
close all
subplot(1,3,1)
histogram(arr1, 'Normalization','probability')
xlabel("SNR (dB)")
ylabel("Occurrence frequency")
title("Data for CGN")
subplot(1,3,2)
histogram(arr2, 'Normalization','probability')
xlabel("SNR (dB)")
title("Data for FGN")
subplot(1,3,3)
histogram(arr3, 'Normalization','probability')
xlabel("SNR (dB)")
title("Data for DGN")
%% No. 21 Bootstrapping
clear
clc
inpath = './experiment/model_parameters/';
type = 2; %D type
cd = 3;
T = 100;
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
y = sim_VAR(model.A,T,1,model.seed,0);
% M = jointvargc(y,parameter.cvx,ALG_PARAMETER.cvx);
for ii = 1:50
    t1 = tic;
    y_bootstrap = y(:,ii:T-50+ii-1,:);
    disp(length(ii:T-50+ii-1))
    M_bootstrap(ii) = jointvargc(y_bootstrap,parameter.cvx,ALG_PARAMETER.cvx);
    t2 = toc(t1);
    disp(fprintf("Time elapsed: %.2f seconds", t2))
end
%% No. 21 Bootstraping
GTmodel = model;
PERF_INDEX = zeros(50, 5);
for bb=1:50
    
    model_acc = performance_eval(M_bootstrap(bb),GTmodel);
    toggle_list = {'total','common','differential'};
    %         M.index.bic=best_index(jj);
    name_list = {'eBIC'};
    for kk=1:length(name_list)
        R.index.(name_list{kk}) = M_bootstrap(bb).index.(name_list{kk});
    end
    toggle = toggle_list{1};
    R.(toggle).F1 =model_acc(M_bootstrap(bb).index.eBIC).(toggle).F1;
    R.(toggle).MCC =model_acc(M_bootstrap(bb).index.eBIC).(toggle).MCC;
    R.(toggle).TPR =model_acc(M_bootstrap(bb).index.eBIC).(toggle).TPR;
    R.(toggle).FPR =model_acc(M_bootstrap(bb).index.eBIC).(toggle).FPR;
    R.(toggle).ACC =model_acc(M_bootstrap(bb).index.eBIC).(toggle).ACC;
    PERF_INDEX(bb, :) = [R.(toggle).F1 R.(toggle).MCC R.(toggle).TPR R.(toggle).FPR R.(toggle).ACC];
    fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n', ...
        mean(R.total.F1),mean(R.total.MCC),mean(R.total.ACC),mean(R.total.FPR),mean(R.total.TPR))
end
%% No. 21 Bootstraping
% clear
% clc
load("H:\Pycharm_projects\bootstrapping_results.mat")
tt = tiledlayout(5, 1);
for kk=1:5
    nexttile;
    histogram(PERF_INDEX(:,kk), 'Normalization', 'Probability')
end

%% load data
clear
clc
load("H:\Pycharm_projects\reviewer_bootstrapping_data.mat")
%% No. 21 Bootstraping
% cnt = 0;
% figure(1)
% selected_model = M_bootstrap(1).model(M_bootstrap(1).index.eBIC);
% plot_group_GC(selected_model.GC)
% sgtitle("Bootstrapping sample: #1")
%
% figure(2)
% selected_model = M_bootstrap(25).model(M_bootstrap(25).index.eBIC);
% plot_group_GC(selected_model.GC)
% sgtitle("Bootstrapping sample: #25")
%
% figure(3)
% selected_model = M_bootstrap(50).model(M_bootstrap(50).index.eBIC);
% plot_group_GC(selected_model.GC)
% sgtitle("Bootstrapping sample: #50")
model_ind = [1, 25, 50];
cnt = 0;
clf
close all
for ii=1:3
    selected_model = M_bootstrap(model_ind(ii)).model(M_bootstrap(model_ind(ii)).index.eBIC);
    GC = selected_model.GC;
    [n,~,K] = size(GC);
    commonNZ = ones(n,n)-eye(n);
    diag_ind = 1:n+1:n^2;
    for kk=1:K
        tmp = GC(:,:,kk);
        tmp(diag_ind) = 0;
        commonNZ = commonNZ & (tmp~=0);
    end
    s = [];
    for kk =1:K
        cnt=cnt+1;
        tmp = GC(:,:,kk);
        tmp(diag_ind) = 0;
        subplot(3,K,cnt)
        spy(tmp,'r')
        hold on
        spy(commonNZ,'k')
        hold off
        axis('square')
        set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
    end
    title(sprintf("Bootstrapping sample: #%d", model_ind(ii)))
    hold off
end

%% No. 32 check total density of ncvx vs cvx model

clear
clc
inpath = './experiment/model_parameters/';

T = 150;
p_true = 1;
p_est = 1;
K = 5;
n = 20;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
m= size(E,2);
realz = m;
GridSize = 30;
density.DGN = zeros(realz, 1);
density.FGN = zeros(realz, 1);


for jj=1:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    nz_number = 0;
    for kk = 1:5
        nz_number = nz_number +length(model.ind{kk});
    end
    density.DGN(jj) = nz_number/((n^2-n)*K);
    
end




for jj=1:realz
    % generate data from given seed
    model = E{3,jj}; % type S
    nz_number = 0;
    for kk = 1:5
        nz_number = nz_number +length(model.ind{kk});
    end
    density.FGN(jj) = nz_number/((n^2-n)*K);
    
end

tt = tiledlayout(1, 2);
nexttile;
histogram(density.DGN, 'Normalization', 'Probability')
title('n=20, p=3, K=5 (DGN)')
xlabel("model density")
ylabel("Occurrence")


nexttile;
histogram(density.FGN, 'Normalization', 'Probability')
title('n=20, p=3, K=5 (FGN)')
xlabel("model density")
ylabel("Occurrence")

%% ROC


clear
clc
clf
close all
type_list = {'total','common','differential'};
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
method_path = {'CGN_ALL_RESULT.mat','CGN_CVX_ALL_RESULT.mat'};


load([resource_path,'CGN_result.mat'])
R1 = R;

load([resource_path,'CGN_CVX_result.mat'])
R2 = R;

DATA_INDEX = 1:1;
ii =2;
for mm=1:2
    load([resource_path,method_path{mm}])
    
    
    tmp_avg = zeros(2,30);
    max1_avg = 0;
    max2_avg = 0;
    for jj=DATA_INDEX
        sample_acc = ALL_RESULT(1,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(1,:)=[metrics.F1];max1 = {max(tmp(1,:))};
        
        sample_acc = ALL_RESULT(2,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(2,:)=[metrics.F1];max2 = {max(tmp(2,:))};
        
        tmp_avg = tmp_avg+tmp/length(DATA_INDEX);
        max1_avg = max1_avg+max1{1}/length(DATA_INDEX);
        max2_avg = max2_avg+max2{1}/length(DATA_INDEX);
    end
    F1_array{mm} = tmp_avg;
    
    tmp_avg = zeros(2,30);
    max1_avg = 0;
    max2_avg = 0;
    for jj=DATA_INDEX
        sample_acc = ALL_RESULT(1,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(1,:)=[metrics.FPR];max1 = {max(tmp(1,:))};
        
        sample_acc = ALL_RESULT(2,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(2,:)=[metrics.FPR];max2 = {max(tmp(2,:))};
        
        tmp_avg = tmp_avg+tmp/length(DATA_INDEX);
        max1_avg = max1_avg+max1{1}/length(DATA_INDEX);
        max2_avg = max2_avg+max2{1}/length(DATA_INDEX);
    end
    FPR_array{mm} = tmp_avg;
    
    tmp_avg = zeros(2,30);
    max1_avg = 0;
    max2_avg = 0;
    for jj=DATA_INDEX
        sample_acc = ALL_RESULT(1,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(1,:)=[metrics.TPR];max1 = {max(tmp(1,:))};
        
        sample_acc = ALL_RESULT(2,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(2,:)=[metrics.TPR];max2 = {max(tmp(2,:))};
        
        tmp_avg = tmp_avg+tmp/length(DATA_INDEX);
        max1_avg = max1_avg+max1{1}/length(DATA_INDEX);
        max2_avg = max2_avg+max2{1}/length(DATA_INDEX);
    end
    TPR_array{mm} = tmp_avg;
end
plot(100*[FPR_array{1};FPR_array{2}]',100*[TPR_array{1};TPR_array{2}]','Linewidth',4)
hold on
index = R1.index(1, DATA_INDEX).eBIC;
scatter(100*FPR_array{1}(1,index),100*TPR_array{1}(1,index),72,[0 0.7 0], 'filled');
a1 = scatter(100*FPR_array{1}(1,:),100*TPR_array{1}(1,:),36,[0 0 1]);
% text(100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), cellstr(string(100*F1_array{1}(1,:)) ))
% datatip(a1,100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), dataTipTextRow('F1',100*F1_array{1}(1,:)) )

row = dataTipTextRow('F1',100*F1_array{1}(1,:),'%+4.4g');
a1.DataTipTemplate.DataTipRows(end+1) = row;


index = R1.index(2, DATA_INDEX).eBIC;
scatter(100*FPR_array{1}(2,index),100*TPR_array{1}(2,index),72,[0 0.7 0], 'filled');
a2 = scatter(100*FPR_array{1}(2,:),100*TPR_array{1}(2,:),36,[0 0 1]);
% text(100*FPR_array{1}(2,:),100*TPR_array{1}(2,:), cellstr(string(100*F1_array{1}(2,:)) ))
% datatip(a1,100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), 'F1', 100*F1_array{1}(1,:))
row = dataTipTextRow('F1',100*F1_array{1}(2,:),'%+4.4g');
a2.DataTipTemplate.DataTipRows(end+1) = row;


index = R2.index(1, DATA_INDEX).eBIC;
scatter(100*FPR_array{2}(1,index),100*TPR_array{2}(1,index),72,[0 0.7 0], 'filled');
a3 = scatter(100*FPR_array{2}(1,:),100*TPR_array{2}(1,:),36,[0 0 1]);
% text(100*FPR_array{2}(1,:),100*TPR_array{2}(1,:), cellstr(string(100*F1_array{2}(1,:)) ))
% datatip(a1,100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), 'F1', 100*F1_array{1}(1,:))
row = dataTipTextRow('F1',100*F1_array{2}(1,:),'%+4.4g');
a3.DataTipTemplate.DataTipRows(end+1) = row;


index = R2.index(2, DATA_INDEX).eBIC;
scatter(100*FPR_array{2}(2,index),100*TPR_array{2}(2,index),72,[0 0.7 0], 'filled');
a4 = scatter(100*FPR_array{2}(2,:),100*TPR_array{2}(2,:),36,[0 0 1]);
% text(100*FPR_array{2}(2,:),100*TPR_array{2}(2,:), cellstr(string(100*F1_array{2}(2,:)) ))
% datatip(a1,100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), 'F1', 100*F1_array{1}(1,:))
row = dataTipTextRow('F1',100*F1_array{2}(2,:),'%+4.4g');
a4.DataTipTemplate.DataTipRows(end+1) = row;



hold off


legend('CGN on common density 10%','CGN on common density 20%','cvx-CGN on common density 10%','cvx-CGN on common density 20%','location','southeast')
% u = gca;
% color_list = {[0,0,1],[0,1,0],[0,0,0],[0,0,0]};
% for kk=1:4
% u.Children(kk).Color = color_list{kk};
% end
% legend
axis([0 30 70 100])
axis([0, 100, 0, 130])
% legend('Common density 10%','Common density 20%','location','northeast')
% title(sprintf('10%%, F1 bestcase: %2.1f, 20%%,  F1 bestcase: %2.1f',100*max1_avg{1},100*max2_avg{1}))
% ylabel('F1 score')
% xlabel('denser $\leftarrow\lambda\rightarrow$ sparser','Interpreter','latex')
xlabel('FPR (%)')
ylabel('TPR (%)')
grid on
set(gca,'FontSize',28)
% pp = get(0, 'Screensize');
% pp(3) = pp(3)*0.75;
% set(gcf, 'Position', pp);
% saveas(gcf,[figurepath,'exp_CGN_ROC'])
% print([figurepath,'exp_CGN_ROC'],'-painters','-depsc','-r300')

%% Dickey-Fuller test for real data
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0); % load data in format (n,T,K)
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
y_total = cat(3,y_TDC,y_ADHD_C);
y_total = y_total-mean(y_total,2); % detrend in time axis



%%
for ii=1:116
    subplot(2,1,1)
    autocorr(y_total(ii,:, 1))
    subplot(2,1,2)
    parcorr(y_total(ii,:, 1))
    pause()
end


%%
% clear
% clc
% inpath = './experiment/model_parameters/';
% type = 2; %D type
% cd = 3;
% T = 1000;
% p_true = 1;
% p_est = 1;
% K = 5;
% K = 50;
% load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E

% model = E{type,3,2,1};
% y = sim_VAR(model.A,T,1,model.seed,1);
% t = 1:1000;
% z = 1/2*t + t.*randn(size(t));
% plot(t,z)

for ii=1:116
    subplot(2,1,1)
    autocorr(y_total(ii,:))
    subplot(2,1,2)
    parcorr(y_total(ii,:))
    pause()
end

%% ROC curve



clear
clc
clf
close all
type_list = {'total','common','differential'};
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
method_path = {'CGN_ALL_RESULT.mat','CGN_CVX_ALL_RESULT.mat'};


load([resource_path,'CGN_result.mat'])
R1 = R;

load([resource_path,'CGN_CVX_result.mat'])
R2 = R;
sample_all = [1, 31, 61, 91];
for oo = sample_all
DATA_INDEX = oo:oo;
ii =2;
for mm=1:2
    load([resource_path,method_path{mm}])


    tmp_avg = zeros(2,30);
    max1_avg = 0;
    max2_avg = 0;
    for jj=DATA_INDEX
        sample_acc = ALL_RESULT(1,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(1,:)=[metrics.F1];max1 = {max(tmp(1,:))};
        
        sample_acc = ALL_RESULT(2,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(2,:)=[metrics.F1];max2 = {max(tmp(2,:))};
        
        tmp_avg = tmp_avg+tmp/length(DATA_INDEX);
        max1_avg = max1_avg+max1{1}/length(DATA_INDEX);
        max2_avg = max2_avg+max2{1}/length(DATA_INDEX);
    end
    F1_array{mm} = tmp_avg;
    
    tmp_avg = zeros(2,30);
    max1_avg = 0;
    max2_avg = 0;
    for jj=DATA_INDEX
        sample_acc = ALL_RESULT(1,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(1,:)=[metrics.FPR];max1 = {max(tmp(1,:))};
        
        sample_acc = ALL_RESULT(2,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(2,:)=[metrics.FPR];max2 = {max(tmp(2,:))};
        
        tmp_avg = tmp_avg+tmp/length(DATA_INDEX);
        max1_avg = max1_avg+max1{1}/length(DATA_INDEX);
        max2_avg = max2_avg+max2{1}/length(DATA_INDEX);
    end
    FPR_array{mm} = tmp_avg;
    
    tmp_avg = zeros(2,30);
    max1_avg = 0;
    max2_avg = 0;
    for jj=DATA_INDEX
        sample_acc = ALL_RESULT(1,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(1,:)=[metrics.TPR];max1 = {max(tmp(1,:))};
        
        sample_acc = ALL_RESULT(2,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(2,:)=[metrics.TPR];max2 = {max(tmp(2,:))};
        
        tmp_avg = tmp_avg+tmp/length(DATA_INDEX);
        max1_avg = max1_avg+max1{1}/length(DATA_INDEX);
        max2_avg = max2_avg+max2{1}/length(DATA_INDEX);
    end
    TPR_array{mm} = tmp_avg;
end

figure(1)
% plot(100*[FPR_array{1};FPR_array{2}]',100*[TPR_array{1};TPR_array{2}]','Linewidth',4)

subplot(2,2,1)
plot(100*FPR_array{1}(1,:),100*TPR_array{1}(1, :),'Linewidth',4)
hold on
index = R1.index(1, DATA_INDEX).eBIC;
scatter(100*FPR_array{1}(1,index),100*TPR_array{1}(1,index),72,[0 0.7 0], 'filled');
a1 = scatter(100*FPR_array{1}(1,:),100*TPR_array{1}(1,:),36,[0 0 1]);
% text(100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), cellstr(string(100*F1_array{1}(1,:)) ))
% datatip(a1,100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), dataTipTextRow('F1',100*F1_array{1}(1,:)) )

row = dataTipTextRow('F1',100*F1_array{1}(1,:),'%+4.4g');
a1.DataTipTemplate.DataTipRows(end+1) = row;
axis([0 30 70 100])
axis('square')
xlabel('FPR (%)')
ylabel('TPR (%)')
grid on
title(sprintf("CGN, d=10%%, F1=%.2f %%", 100*F1_array{1}(1,index)))
hold off


subplot(2,2,2)
plot(100*FPR_array{1}(2,:),100*TPR_array{1}(2, :),'Linewidth',4)
hold on
index = R1.index(2, DATA_INDEX).eBIC;
scatter(100*FPR_array{1}(2,index),100*TPR_array{1}(2,index),72,[0 0.7 0], 'filled');
a2 = scatter(100*FPR_array{1}(2,:),100*TPR_array{1}(2,:),36,[0 0 1]);
% text(100*FPR_array{1}(2,:),100*TPR_array{1}(2,:), cellstr(string(100*F1_array{1}(2,:)) ))
% datatip(a1,100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), 'F1', 100*F1_array{1}(1,:))
row = dataTipTextRow('F1',100*F1_array{1}(2,:),'%+4.4g');
a2.DataTipTemplate.DataTipRows(end+1) = row;
axis([0 30 70 100])
axis('square')
xlabel('FPR (%)')
ylabel('TPR (%)')
grid on
title(sprintf("CGN, d=20%%, F1=%.2f %%", 100*F1_array{1}(2,index)))
hold off


subplot(2,2,3)
plot(100*FPR_array{2}(1,:),100*TPR_array{2}(1, :),'Linewidth',4)
hold on
index = R2.index(1, DATA_INDEX).eBIC;
scatter(100*FPR_array{2}(1,index),100*TPR_array{2}(1,index),72,[0 0.7 0], 'filled');
a3 = scatter(100*FPR_array{2}(1,:),100*TPR_array{2}(1,:),36,[0 0 1]);
% text(100*FPR_array{2}(1,:),100*TPR_array{2}(1,:), cellstr(string(100*F1_array{2}(1,:)) ))
% datatip(a1,100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), 'F1', 100*F1_array{1}(1,:))
row = dataTipTextRow('F1',100*F1_array{2}(1,:),'%+4.4g');
a3.DataTipTemplate.DataTipRows(end+1) = row;
axis([0 30 70 100])
axis('square')
xlabel('FPR (%)')
ylabel('TPR (%)')
grid on
title(sprintf("cvx-CGN, d=10%%, F1=%.2f %%", 100*F1_array{2}(1,index)))
hold off


subplot(2,2,4)
plot(100*FPR_array{2}(2,:),100*TPR_array{2}(2, :),'Linewidth',4)
hold on
index = R2.index(2, DATA_INDEX).eBIC;
scatter(100*FPR_array{2}(2,index),100*TPR_array{2}(2,index),72,[0 0.7 0], 'filled');
a4 = scatter(100*FPR_array{2}(2,:),100*TPR_array{2}(2,:),36,[0 0 1]);
% text(100*FPR_array{2}(2,:),100*TPR_array{2}(2,:), cellstr(string(100*F1_array{2}(2,:)) ))
% datatip(a1,100*FPR_array{1}(1,:),100*TPR_array{1}(1,:), 'F1', 100*F1_array{1}(1,:))
row = dataTipTextRow('F1',100*F1_array{2}(2,:),'%+4.4g');
a4.DataTipTemplate.DataTipRows(end+1) = row;
axis([0 30 70 100])
axis('square')
xlabel('FPR (%)')
ylabel('TPR (%)')
grid on
title(sprintf("cvx-CGN, d=20%%, F1=%.2f %%", 100*F1_array{2}(2,index)))



hold off


% legend('CGN on common density 10%','CGN on common density 20%','cvx-CGN on common density 10%','cvx-CGN on common density 20%','location','southeast')
% u = gca;
% color_list = {[0,0,1],[0,1,0],[0,0,0],[0,0,0]};
% for kk=1:4
% u.Children(kk).Color = color_list{kk};
% end
% legend

% set(gca,'FontSize',28)
set(findall(gcf,'-property','FontSize'),'FontSize',18)
pp = get(0, 'Screensize');
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
saveas(gcf,[figurepath,'exp_CGN_ROC'])
print([figurepath,'reviewer_response_exp_CGN_ROC_sample', int2str(oo)],'-dpng','-r300')

end

%% fMRI ACF, PACF
figurepath = './results2plot/figures/';

cnt = 0;
clear fMRI
for ii=1:116
    for jj=1:36
        cnt=cnt+1;
        fMRI.pacf(cnt,:) = parcorr(y_total(ii,:,jj), 'NumLags', 10);
        fMRI.acf(cnt,:) = xcov(y_total(ii,:,jj), 10, 'coeff');

    end
end

tt = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

patient_index = 5;
ax1=nexttile;


stem(abs(fMRI.pacf(patient_index,2:end)))
xlabel('Lags')
ylabel('value (scaled to 1)')
title('Partial autocorrelation (sample)')
ax2=nexttile;
stem(abs(fMRI.acf(patient_index,12:end)))
xlabel('Lags')
ylabel('value (scaled to 1)')
title('Autocorrelation (sample)')

ax3=nexttile;
boxplot(abs(fMRI.pacf(:,2:end)))
xlabel('Lags')
ylabel('value (scaled to 1)')
title('Partial autocorrelation (distribution)')
ax4=nexttile;
boxplot(abs(fMRI.acf(:,12:end)))
xlabel('Lags')
ylabel('value (scaled to 1)')
title('Autocorrelation (distribution)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
pp = get(0, 'Screensize');
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);

linkaxes([ax1,ax2, ax3, ax4],'x');


set(findall(gcf,'-property','FontSize'),'FontSize',18)
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
print([figurepath,'reviewer_response_fMRI_ACF_PACF'],'-dpng','-r300')
print([figurepath,'reviewer_response_fMRI_ACF_PACF'],'-depsc','-r300')


