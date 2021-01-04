[data,header,raw] = tsvread('C:\Users\CU_EE_LAB408\Desktop\sfnwmrda0010001_session_1_rest_1_aal_TCs.1D');
load('./experiment/experiment_real_data/AAL_116.mat')

y = data(2:end,3:end)';



y_DMN = y(AAL_116.DMN,:);
n = size(y_DMN,1);
p=1;
K=size(y_DMN,3);
[P,~] = offdiagJSS(size(y,1),p,size(y,3));
[P_DMN,~] = offdiagJSS(n,p,K);
M = formulation_C(y_DMN,P_DMN,p,30);
% AAL_ATLAS = 
M_original = formulation_C(y,P,p,30);