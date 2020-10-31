% clear
clc

load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\result_test_formulationD_5percent_5_p1.mat')
load('D:\JGranger_ncvx\data_compare\test_model_K5_p1.mat')
type = 2;
cd = 2;
dd = 2;
realz = 5;


GTmodel = E{type,cd,dd,realz};
n = GTmodel.dim(1);
p = GTmodel.dim(2);
K = GTmodel.dim(3);

bic_index = M.index.bic;
estimated_model = M.model(bic_index);

GC_estim = estimated_model.GC;
TP=0;
TN=0;
FP=0;
FN=0;

for kk=1:K
    ind{kk} = find(GC_estim(:,:,kk));
    ind{kk} = setdiff(ind{kk},1:n+1:n^2);
    disp(length(setdiff(ind{kk},estimated_model.ind{1}{kk})))
    ind_true{kk} = find(GTmodel.GC(:,:,kk));
    TP = TP+ length(intersect(ind{kk},ind_true{kk}));
    FN = FN+ length(setdiff(ind_true{kk},ind{kk}));
    FP = FP+ length(setdiff(ind{kk},ind_true{kk}));
end
TN = (n^2-n)*K-TP-FN-FP;


score.FPR = FP./(FP+TN);
score.TPR = TP./(TP+FN);
score.ACC = (TP+TN)./(TP+FN+TN+FP);
score.F1 = 2*TP./(2*TP+FP+FN);


[stat] = compare_sparsity(ind_true,ind,n,K,'single_differential');

performance_score(squeeze(stat))
%%
T=100;

GTmodel = E{type,cd,dd,realz};
y = sim_VAR(GTmodel.A,T,1,GTmodel.seed,0);
[P,~] = offdiagJSS(n,p,K);
% M_test = formulation_D(y,P,p,30);
% save('./temporal_file/estimated_model_demo','M_test')
load('./temporal_file/estimated_model_demo')
%%


bic_index = M_test.index.bic;


for ii=1:30
    for jj=1:30
        
        estimated_model = M_test.model(ii,jj);

GC_estim = estimated_model.GC;
TP=0;
TN=0;
FP=0;
FN=0;

for kk=1:K
    ind{kk} = find(GC_estim(:,:,kk));
    ind{kk} = setdiff(ind{kk},1:n+1:n^2);
    disp(length(setdiff(ind{kk},estimated_model.ind{1}{kk})))
    ind_true{kk} = find(GTmodel.GC(:,:,kk));
    TP = TP+ length(intersect(ind{kk},ind_true{kk}));
    FN = FN+ length(setdiff(ind_true{kk},ind{kk}));
    FP = FP+ length(setdiff(ind{kk},ind_true{kk}));
end
TN = (n^2-n)*K-TP-FN-FP;


score.FPR = FP./(FP+TN);
score.TPR = TP./(TP+FN);
score.ACC = (TP+TN)./(TP+FN+TN+FP);
score.F1 = 2*TP./(2*TP+FP+FN);


[stat] = compare_sparsity(ind_true,ind,n,K,'single_differential');

result = performance_score(squeeze(stat));
FPR(ii,jj) = result.FPR;
TPR(ii,jj) = result.TPR;
    end
end
%%
figure;
plot(FPR,TPR)
hold on
scatter(FPR(bic_index),TPR(bic_index),'xr')
hold off
%%

plot_group_GC(GC_estim)
sgtitle('estimated GC (bic selected)')
plot_group_GC(GTmodel.GC)
sgtitle('ground truth GC')