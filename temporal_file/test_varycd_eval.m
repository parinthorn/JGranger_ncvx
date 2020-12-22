clear
clc
total_r = 50;
F1_avg = 0;
n= 20;p=1;
type = 2;
ii=1;
K=50;
model_name = {'2','4','5','6','8'};
% model_name = {'5','15','25'};
for ii=1:5
for realz = 1:total_r
    fprintf('realization:%d \n',realz)
load(['G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/cvx_CDvary_result_formulationD_c',model_name{ii},'d',model_name{5-ii+1},'percent_lag1_K',int2str(K),'_',int2str(realz),'.mat'])
load(['./data_compare/vary_CD_model',int2str(K),'_p1.mat'])

stat = [0,0,0,0];

% realz = 2;

for kk=1:K
    TP = length(intersect(E{type,ii,ii,realz}.ind{kk},M.model(M.index.bic).ind{1}{kk}));
%     FP = intersect(M.model(M.index.bic).ind{kk},E{2,1,1,1}.ind{kk});
    FP = length(setdiff(M.model(M.index.bic).ind{1}{kk},E{type,ii,ii,realz}.ind{kk}));
    FN = length(setdiff(E{type,ii,ii,realz}.ind{kk},M.model(M.index.bic).ind{1}{kk}));
    TN = n^2-n-TP-FP-FN;
    stat = stat+[TP FP FN TN];
end
TP = stat(1);
FP = stat(2);
FN = stat(3);
TN = stat(4);

score(ii).TPR(realz) = TP/(TP+FN);
score(ii).ACC(realz) = (TP+TN)/(TP+FN+TN+FP);
score(ii).F1(realz) = 2*TP/(2*TP+FP+FN);
tmp = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
score(ii).MCC(realz) = (TP*TN-FP*FN)/tmp;

end
end