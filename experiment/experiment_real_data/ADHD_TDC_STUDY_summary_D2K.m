%% D2K result analysis
clear
clc
clf
close all

% load estimated model

load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected_LLHcorrected','M')
GCmat.D2K = M.model(M.index.eBIC).GC;

load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected_LLHcorrected','M')
GCmat.S2K = M.model(M.index.eBIC).GC;

load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered_LLHcorrection','M')
GCmat.C18K(:,:,1) = mean(M.TDC.model(M.TDC.index.eBIC).GC,3);
GCmat.C18K(:,:,2) = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3);

load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\summary_real_DS2K_C18K_timecorrected_LLHcorrected')


% Left Cerebellum 3 (95) -> Left Anterior cingulate gyrus (31)
% index = [95,31];
index.DMN = [59,60,61,62,85,86];
index.DMNa = [29,30,31,32,87,88];
index.DMNv = [35,36,37,38,39,40,55,56,65,66,67,68];
index.SM = [1,2,7,8,19,20,57,58,63,64,69,70];
index.Visual = [43,44,45,46,47,48,49,50,51,52,53,54];
index.SN = [7,8,9,10, 29,30, 23,24, 19,20];

index_1 = index.DMN;
index_2 = index.SN;

color_vector = [repmat([1 0 0],[length(index_1),1]);repmat([0 0 1],[length(index_2),1])];

selected_index= [index_1 index_2];
selected_ROI = AAL_116.name(selected_index);


% rr = figure(1);
rr = tiledlayout(3,2);
f1 = [1,3,5];
f2 = [2,4,6];
name_list = {'D2K','S2K','C18K'};
for kk=1:3
% subplot(3,2,f1(kk))
nexttile;
tt = circularGraph(GCmat.(name_list{kk})(selected_index,selected_index,1)','Label', cellstr(num2str(selected_index')),'ColorMap',color_vector);
% for ii=1:length(tt.Node)
%     tt.Node(ii).Visible = false;
% end
title('TDC')
nexttile;
% figure(id+2)
% subplot(3,2,f2(kk))
tt= circularGraph(GCmat.(name_list{kk})(selected_index,selected_index,2)','Label', cellstr(num2str(selected_index')),'ColorMap',color_vector);
% for ii=1:length(tt.Node)
%     tt.Node(ii).Visible = false;
% end
title('ADHD')
end
rr.Children(5).Title.Visible = 'on';
rr.Children(6).Title.Visible = 'on';