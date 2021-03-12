% A Physarum Centrality Measure of the Human Brain Network
clear
clc
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
% AAL_116.name = {'PreCG_L','PreCG_R','SFGdor_L','SFGdor_R','ORBsup_L','ORBsup_R','MFG_L','MFG_R','ORBmid_L','ORBmid_R','IFGoperc_L','IFGoperc_R','IFGtriang_L','IFGtriang_R','ORBinf_L','ORBinf_R','ROL_L','ROL_R','SMA_L','SMA_R','OLF_L','OLF_R','SFGmed_L','SFGmed_R','ORBsupmed_L','ORBsupmed_R','REC_L','REC_R','INS_L','INS_R','ACG_L','ACG_R','MCG_L','MCG_R','PCG_L','PCG_R','HIP_L','HIP_R','PHG_L','PHG_R','AMYG_L','AMYG_R','CAL_L','CAL_R','CUN_L','CUN_R','LING_L','LING_R','SOG_L','SOG_R','MOG_L','MOG_R','IOG_L','IOG_R','FFG_L','FFG_R','PoCG_L','PoCG_R','SPG_L','SPG_R','IPG_L','IPG_R','SMG_L','SMG_R','ANG_L','ANG_R','PCUN_L','PCUN_R','PCL_L','PCL_R','CAU_L','CAU_R','PUT_L','PUT_R','PAL_L','PAL_R','THA_L','THA_R','HES_L','HES_R','STG_L','STG_R','TPOsup_L','TPOsup_R','MTG_L','MTG_R','TPOmid_L','TPOmid_R','ITG_L','ITG_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R','Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R','Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R','Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R','Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R','Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8','Vermis_9','Vermis_10'};
% AAL_116.name_full = {'Precentral gyrus_L','Precentral gyrus_R','Superior frontal gyrus (dorsolateral)_L','Superior frontal gyrus (dorsolateral)_R','Superior frontal gyrus (orbital)_L','Superior frontal gyrus (orbital)_R','Middle frontal gyrus_L','Middle frontal gyrus_R','Middle frontal gyrus (orbital)_L','Middle frontal gyrus (orbital)_R','Inferior frontal gyrus (opercular)_L','Inferior frontal gyrus (opercular)_R','Inferior frontal gyrus (triangular)_L','Inferior frontal gyrus (triangular)_R','Inferior frontal gyrus (orbital)_L','Inferior frontal gyrus (orbital)_R','Rolandic operculum_L','Rolandic operculum_R','Supplementary motor area_L','Supplementary motor area_R','Olfactroy cortex_L','Olfactroy cortex_R','Superior frontal gyrus (medial)_L','Superior frontal gyrus (medial)_R','Superior frontal gyrus (medial orbital)_L','Superior frontal gyrus (medial orbital)_R','Rectus gyrus_L','Rectus gyrus_R','Insula_L','Insula_R','Anterior cingulate gyrus_L','Anterior cingulate gyrus_R','Median cingulate gyrus_L','Median cingulate gyrus_R','Posterior cingulate gyrus_L','Posterior cingulate gyrus_R','Hippocampus_L','Hippocampus_R','Parahippocampal gyrus_L','Parahippocampal gyrus_R','Amygdala_L','Amygdala_R','Calcarine cortex_L','Calcarine cortex_R','Cuneus_L','Cuneus_R','Lingual gyrus_L','Lingual gyrus_R','Superior occipital gyrus_L','Superior occipital gyrus_R','Middle occipital gyrus_L','Middle occipital gyrus_R','Inferior occipital gyrus_L','Inferior occipital gyrus_R','Fusiform gyrus_L','Fusiform gyrus_R','Postcentral gyrus_L','Postcentral gyrus_R','Superior parietal gyrus_L','Superior parietal gyrus_R','Inferior parietal gyrus_L','Inferior parietal gyrus_R','Supramarginal gyrus_L','Supramarginal gyrus_R','Angular gyrus_L','Angular gyrus_R','Precuneus_L','Precuneus_R','Paracentral lobule_L','Paracentral lobule_R','Caudate_L','Caudate_R','Putamen_L','Putamen_R','Pallidum_L','Pallidum_R','Thalamus_L','Thalamus_R','Heschl gyrus_L','Heschl gyrus_R','Superior temporal gyrus_L','Superior temporal gyrus_R','Temporal pole (superior)_L','Temporal pole (superior)_R','Middle temporal gyrus_L','Middle temporal gyrus_R','Temporal pole (middle)_L','Temporal pole (middle)_R','Inferior temporal gyrus_L','Inferior temporal gyrus_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R','Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R','Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R','Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R','Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R','Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8','Vermis_9','Vermis_10'};
%% D2K
clc
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
G.TDC = digraph(GC.total(:,:,1)',AAL_116.name);
G.ADHD = digraph(GC.total(:,:,2)',AAL_116.name);
% 
C1 = (centrality(G.TDC,'betweenness','Cost',G.TDC.Edges.Weight));
C2 = (centrality(G.ADHD,'betweenness','Cost',G.ADHD.Edges.Weight));

% C1 = (centrality(G.TDC,'betweenness'));
% C2 = (centrality(G.ADHD,'betweenness'));

% C1 = (centrality(G.TDC,'hubs','Importance',G.TDC.Edges.Weight));
% C2 = (centrality(G.ADHD,'hubs','Importance',G.ADHD.Edges.Weight));
% C1 = (centrality(G.TDC,'pagerank','Importance',G.TDC.Edges.Weight));
% C2 = (centrality(G.ADHD,'pagerank','Importance',G.ADHD.Edges.Weight));

C_in = (C1-C2);
C=C_in;
THRESH =2*std(C_in);

C(abs(C_in-mean(C_in))>THRESH)=1;
C(abs(C_in-mean(C_in))<=THRESH)=0;

figure(1)
subplot(1,2,1)
p1 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
p1.NodeCData = C;
p1.MarkerSize = 7;
colormap(1-gray)

subplot(1,2,2)
p2 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
p2.NodeCData = C;
p2.MarkerSize = 7;
colormap(1-gray)
selected_index = find(C);
[~,I] = sort(C_in(selected_index),'descend');
fprintf('Difference nodes are,\n')
for ii=1:length(selected_index)
    fprintf([AAL_116.name_full{selected_index(I(ii))},'\n'])
end
%% S2K
clc
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
G.TDC = digraph(GC.total(:,:,1)',AAL_116.name);
G.ADHD = digraph(GC.total(:,:,2)',AAL_116.name);
C1 = (centrality(G.TDC,'betweenness','Cost',G.TDC.Edges.Weight));
C2 = (centrality(G.ADHD,'betweenness','Cost',G.ADHD.Edges.Weight));

% C1 = (centrality(G.TDC,'betweenness'));
% C2 = (centrality(G.ADHD,'betweenness'));

% C1 = (centrality(G.TDC,'hubs','Importance',G.TDC.Edges.Weight));
% C2 = (centrality(G.ADHD,'hubs','Importance',G.ADHD.Edges.Weight));
% C1 = (centrality(G.TDC,'pagerank','Importance',G.TDC.Edges.Weight));
% C2 = (centrality(G.ADHD,'pagerank','Importance',G.ADHD.Edges.Weight));

C_in = (C1-C2);
C=C_in;
THRESH =2*std(C_in);

C(abs(C_in-mean(C_in))>THRESH)=1;
C(abs(C_in-mean(C_in))<=THRESH)=0;

figure(1)
subplot(1,2,1)
p1 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
p1.NodeCData = C;
p1.MarkerSize = 7;
colormap(1-gray)

subplot(1,2,2)
p2 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
p2.NodeCData = C;
p2.MarkerSize = 7;
colormap(1-gray)
selected_index = find(C);
[~,I] = sort(C_in(selected_index),'descend');
fprintf('Difference nodes are,\n')
for ii=1:length(selected_index)
    fprintf([AAL_116.name_full{selected_index(I(ii))},'\n'])
end
%% C18K
clc
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered_LLHcorrection','M')
GC.total(:,:,1) = mean(M.TDC.model(M.TDC.index.eBIC).GC,3);
GC.total(:,:,2) = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3);
GC.total(GC.total>0) = 1./GC.total(GC.total>0);

G.TDC = digraph(GC.total(:,:,1)',AAL_116.name);
G.ADHD = digraph(GC.total(:,:,2)',AAL_116.name);
C1 = (centrality(G.TDC,'betweenness','Cost',G.TDC.Edges.Weight));
C2 = (centrality(G.ADHD,'betweenness','Cost',G.ADHD.Edges.Weight));

% C1 = (centrality(G.TDC,'betweenness'));
% C2 = (centrality(G.ADHD,'betweenness'));

% C1 = (centrality(G.TDC,'hubs','Importance',G.TDC.Edges.Weight));
% C2 = (centrality(G.ADHD,'hubs','Importance',G.ADHD.Edges.Weight));
% C1 = (centrality(G.TDC,'pagerank','Importance',G.TDC.Edges.Weight));
% C2 = (centrality(G.ADHD,'pagerank','Importance',G.ADHD.Edges.Weight));

C_in = (C1-C2);
C=C_in;
THRESH =2*std(C_in);

C(abs(C_in-mean(C_in))>THRESH)=1;
C(abs(C_in-mean(C_in))<=THRESH)=0;

figure(1)
subplot(1,2,1)
p1 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
p1.NodeCData = C;
p1.MarkerSize = 7;
colormap(1-gray)

subplot(1,2,2)
p2 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
p2.NodeCData = C;
p2.MarkerSize = 7;
colormap(1-gray)
selected_index = find(C);
[~,I] = sort(C_in(selected_index),'descend');
fprintf('Difference nodes are,\n')
for ii=1:length(selected_index)
    fprintf([AAL_116.name_full{selected_index(I(ii))},'\n'])
end