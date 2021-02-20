%% Consensus Result on Real data
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0);
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0);
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
y_total(:,:,1) = reshape(y_TDC,[116,size(y_TDC,2)*size(y_TDC,3)]);
y_total(:,:,2) = reshape(y_ADHD_C,[116,size(y_ADHD_C,2)*size(y_ADHD_C,3)]);
% map
AAL_116.name = {'PreCG_L','PreCG_R','SFGdor_L','SFGdor_R','ORBsup_L','ORBsup_R','MFG_L','MFG_R','ORBmid_L','ORBmid_R','IFGoperc_L','IFGoperc_R','IFGtriang_L','IFGtriang_R','ORBinf_L','ORBinf_R','ROL_L','ROL_R','SMA_L','SMA_R','OLF_L','OLF_R','SFGmed_L','SFGmed_R','ORBsupmed_L','ORBsupmed_R','REC_L','REC_R','INS_L','INS_R','ACG_L','ACG_R','MCG_L','MCG_R','PCG_L','PCG_R','HIP_L','HIP_R','PHG_L','PHG_R','AMYG_L','AMYG_R','CAL_L','CAL_R','CUN_L','CUN_R','LING_L','LING_R','SOG_L','SOG_R','MOG_L','MOG_R','IOG_L','IOG_R','FFG_L','FFG_R','PoCG_L','PoCG_R','SPG_L','SPG_R','IPG_L','IPG_R','SMG_L','SMG_R','ANG_L','ANG_R','PCUN_L','PCUN_R','PCL_L','PCL_R','CAU_L','CAU_R','PUT_L','PUT_R','PAL_L','PAL_R','THA_L','THA_R','HES_L','HES_R','STG_L','STG_R','TPOsup_L','TPOsup_R','MTG_L','MTG_R','TPOmid_L','TPOmid_R','ITG_L','ITG_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R','Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R','Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R','Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R','Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R','Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8','Vermis_9','Vermis_10'};
%% K=2
% formulation D
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered')
M = fix_loglikelihood(M,y_total);
M = augment_score(M,size(y_total,2),'llh_hetero');
result.TDC_index.D2K = M.model(M.index.eBIC).ind{1}{1};
result.ADHD_index.D2K = M.model(M.index.eBIC).ind{1}{2};
clear M
%% K=18
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered')
M.TDC = fix_loglikelihood(M.TDC,y_TDC);
M.TDC = augment_score(M.TDC,size(y_TDC,2),'llh_hetero');
M.ADHD_C = fix_loglikelihood(M.ADHD_C,y_ADHD_C);
M.ADHD_C = augment_score(M.ADHD_C,size(y_ADHD_C,2),'llh_hetero');

result.TDC_index.C18K = M.TDC.model(M.TDC.index.eBIC).ind_common{1};
result.ADHD_index.C18K = M.ADHD_C.model(M.ADHD_C.index.eBIC).ind_common{1};
clear M
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_D_unfiltered')
M.TDC = fix_loglikelihood(M.TDC,y_TDC);
M.TDC = augment_score(M.TDC,size(y_TDC,2),'llh_full');
M.ADHD_C = fix_loglikelihood(M.ADHD_C,y_ADHD_C);
M.ADHD_C = augment_score(M.ADHD_C,size(y_ADHD_C,2),'llh_full');

result.TDC_index.D18K = M.TDC.model(M.TDC.index.eBIC).ind_common{1};
result.ADHD_index.D18K = M.ADHD_C.model(M.ADHD_C.index.eBIC).ind_common{1};

save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\summary_real_fixdf','result')
%%
TDC_common_index = setdiff(1:1:116^2,1:116+1:116^2);
ADHD_common_index = setdiff(1:1:116^2,1:116+1:116^2);
namelist = {'D2K','C18K','D18K'};
for ii=1:3
    TDC_common_index = intersect(TDC_common_index,result.TDC_index.(namelist{ii}));
    ADHD_common_index = intersect(ADHD_common_index,result.ADHD_index.(namelist{ii}));
end

extra_link = setdiff(ADHD_common_index,TDC_common_index);
missing_link = setdiff(TDC_common_index,ADHD_common_index);

missing_grid = zeros(116);
extra_grid = zeros(116);
missing_grid(missing_link)=1;
extra_grid(extra_link)=1;

AAL_grid = zeros(116);
AAL_grid(missing_link) = -1;
AAL_grid(extra_link) = 1;


%% temporal K=2
ii=26;
jj=27;
% found common Cerebellum_9_L -> Vermis_1_2 in both ADHD and TDC
figure(1);
atlas_index = 1:116;
M.TDC.GC = (M.model(ii,jj).GC(:,:,1)).*(1-eye(116));  
M.ADHD_C.GC = (M.model(ii,jj).GC(:,:,2)).*(1-eye(116));  
nametag = {'TDC','ADHD_C'};
titletag = {'TDC','ADHD'};

for ii=1:2

subplot(1,2,ii)
imagesc(M.(nametag{ii}).GC(atlas_index,atlas_index))
grid on
set(gca,'xtick',1:1:length(atlas_index),'ytick',1:1:length(atlas_index), ...
    'xticklabel',AAL_116.name(atlas_index),'yticklabel',AAL_116.name(atlas_index))
axis('square')
title(titletag{ii})
colormap(jet)
caxis([0 1])
xtickangle(30)
end

%% K=18