clf
clc
clear
close all
realz = 76;
load('C:\Users\CU_EE_LAB408\Dropbox\0MASTER\MATLAB_MASTER\JGranger_ncvx\data_compare\model_K5_p1.mat')
GT_model = E(2,3,2,realz);

load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\formulationC_10percent_',int2str(realz),'.mat'])

M_c = M;
load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\magda_formulationC_10percent_',int2str(realz),'.mat'])
M_m =M;
load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\result_formulationD_5percent_lag1_K5_',int2str(realz),'.mat'])
M_d = M;

h = M_c.model(M_c.index.bic).ind{1}{1};
for ii=2:5
h = intersect(h,M_c.model(M_c.index.bic).ind{1}{ii});
end
all((M_c.model(M_c.index.bic).ind_common{1}-h)==0,'all')

plot_group_GC(M_c.model(M_c.index.bic).GC);
plot_group_GC(M_d.model(M_d.index.bic).GC);