%% Group analysis

% LOAD DATA
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu');
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu');
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered')
M.TDC = augment_score(M.TDC,size(y_TDC,2),'llh_hetero');
M.ADHD_C = augment_score(M.ADHD_C,size(y_TDC,2),'llh_hetero');
%% PLOT AVERAGE GC
clf
close all
AAL_116.name = {'PreCG_L','PreCG_R','SFGdor_L','SFGdor_R','ORBsup_L','ORBsup_R','MFG_L','MFG_R','ORBmid_L','ORBmid_R','IFGoperc_L','IFGoperc_R','IFGtriang_L','IFGtriang_R','ORBinf_L','ORBinf_R','ROL_L','ROL_R','SMA_L','SMA_R','OLF_L','OLF_R','SFGmed_L','SFGmed_R','ORBsupmed_L','ORBsupmed_R','REC_L','REC_R','INS_L','INS_R','ACG_L','ACG_R','MCG_L','MCG_R','PCG_L','PCG_R','HIP_L','HIP_R','PHG_L','PHG_R','AMYG_L','AMYG_R','CAL_L','CAL_R','CUN_L','CUN_R','LING_L','LING_R','SOG_L','SOG_R','MOG_L','MOG_R','IOG_L','IOG_R','FFG_L','FFG_R','PoCG_L','PoCG_R','SPG_L','SPG_R','IPG_L','IPG_R','SMG_L','SMG_R','ANG_L','ANG_R','PCUN_L','PCUN_R','PCL_L','PCL_R','CAU_L','CAU_R','PUT_L','PUT_R','PAL_L','PAL_R','THA_L','THA_R','HES_L','HES_R','STG_L','STG_R','TPOsup_L','TPOsup_R','MTG_L','MTG_R','TPOmid_L','TPOmid_R','ITG_L','ITG_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R','Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R','Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R','Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R','Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R','Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8','Vermis_9','Vermis_10'};
AAL_116.DMN = [23,24,31,32,35,36,67,68,65,66];
AAL_116.FPN = [65,66,7,8,11,12,13,14,61,62];
AAL_116.CC = [31,32,33,34,35,36];
AAL_116.sensorimotor_cortex = [1,2,57,58,19,20,77,78];
% atlas_index = [7,8,11,12,13,14,15,16,31,32,35,36,67,68,19,20];
% atlas_index = union(AAL_116.DMN,AAL_116.FPN,'stable');
% atlas_index = [7:10, 12,14,16]; %MFG and IFG_R extra link [NOT FOUND]
% atlas_index = [32,35,36]; % R_dACC and PCC Decrease (32 -> 35,36) [FOUND MISSING LINK PCC -> R_dACC]
% atlas_index = [19,20,7:10]; %SMA -> MFG Increase (19,20 -> 7:10) [found extra links from both SMA_L,R to MFG Orbital]
% atlas_index = [11:16,19,20]; %IFG -> SMA Increase (11:16 -> 19,20)

index_list{1} = AAL_116.DMN;
index_list{2} = AAL_116.FPN;
index_list{3} = union(AAL_116.DMN,AAL_116.FPN,'stable');
index_list{4} = AAL_116.CC;
index_list{5} = [7:10, 11,12,13,14,15,16]; %MFG and IFG_R extra link
index_list{6} = [32,35,36,67,68]; % R_dACC and PCC/Precuneus Decrease (32 -> 35,36)
index_list{7} = [19,20,7:10]; %SMA -> MFG Increase (19,20 -> 7:10)
index_list{8} = [11:16,19,20]; %IFG -> SMA Increase (11:16 -> 19,20)
index_list{9} = union(AAL_116.DMN,AAL_116.sensorimotor_cortex,'stable'); % found abnormal connectivity in ADHD type C
index_list{10} = [60, 9];  %Extra SPG_R -> ORB_mid_L
index_list{11} = [18, 93, 97, 49]; % Extra SOG_L -> Cerebellum 4,5_L
index_list{12} = [30, 47];  % missing   Lingual L -> Insular R
index_list{13} = [46, 55];% Missing FFG_L -> Cuneus R
index_list{14} = [61, 63]; % no
index_list{15} = [67, 89]; % no
index_list{16} = [85, 115];% Missing MTG_L -> Vermis 9
index_list{17} = [67, 113];% Missing Precuneus L -> Vermis 7
network_name = {'Default Mode Network','Fronto-Parietal Network','DMN&FPN','Cingulate Cortex','Lit1','Lit2','Lit3','Lit4','Lit5'};

% figure(1)
% tt= tiledlayout(1,2);
M.TDC.GC_avg = mean(M.TDC.model(M.TDC.index.eBIC).GC,3).*(1-eye(116));  
M.ADHD_C.GC_avg = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3).*(1-eye(116));  
nametag = {'TDC','ADHD_C'};
titletag = {'TDC','ADHD'};
for id = 1:length(index_list)
    atlas_index = index_list{id};
    figure(id)
for ii=1:2
% nexttile;

subplot(1,2,ii)
imagesc(M.(nametag{ii}).GC_avg(atlas_index,atlas_index))
grid on
set(gca,'xtick',1:1:length(atlas_index),'ytick',1:1:length(atlas_index), ...
    'xticklabel',AAL_116.name(atlas_index),'yticklabel',AAL_116.name(atlas_index))
axis('square')
title(titletag{ii})
colormap(jet)
caxis([0 0.3])
xtickangle(30)
end
end


atlas_index = index_list{2};
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
TDC_GC_print = M.TDC.GC_avg(atlas_index,atlas_index);
ADHD_GC_print = M.ADHD_C.GC_avg(atlas_index,atlas_index);
tmp = AAL_116.name((atlas_index));
figure(id+1)
subplot(1,2,1)
tt = circularGraph(TDC_GC_print','Label',tmp);
for ii=1:length(tt.Node)
    tt.Node(ii).Visible = false;
end
title('TDC')
% figure(id+2)
subplot(1,2,2)
tt= circularGraph(ADHD_GC_print','Label',tmp);
for ii=1:length(tt.Node)
    tt.Node(ii).Visible = false;
end
title('ADHD')