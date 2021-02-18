% [data,header,raw] = tsvread('C:\Users\CU_EE_LAB408\Desktop\sfnwmrda0010001_session_1_rest_1_aal_TCs.1D');
% load('./experiment/experiment_real_data/AAL_116.mat')
% y = data(2:end,3:end)';

%% ADHD formulation C
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0); %unfilter=0
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0);
K = size(y_TDC,3);

%%
y_TDC_concat = reshape(y_TDC,[116,172*18]);
y_ADHD_C_concat = reshape(y_ADHD_C,[116,172*18]);
% %%
% figure(1)
% subplot(2,2,1)
% imagesc(fcon_TDC.*(abs(fcon_TDC.*(1-eye(116)))>0.2))
% axis('square')
% subplot(2,2,2)
% imagesc(fcon_ADHD_C.*(abs(fcon_ADHD_C.*(1-eye(116)))>0.2))
% axis('square')
% 
% subplot(2,2,3)
% imagesc(fcon_TDC_concat.*(abs(fcon_TDC_concat.*(1-eye(116)))>0.4))
% axis('square')
% subplot(2,2,4)
% imagesc(fcon_ADHD_C_concat.*(abs(fcon_ADHD_C_concat.*(1-eye(116)))>0.4))
% axis('square')
%%
y_TDC_concat = detrend(y_TDC_concat')';
y_ADHD_C_concat = detrend(y_ADHD_C_concat')';
y_total(:,:,1) = y_TDC_concat;
y_total(:,:,2) = y_ADHD_C_concat;
%%
imagesc(squeeze(std(y_total,[],2)))
%%
for ii=1:116
    figure(1)
    plot([y_total(ii,:,1)' y_total(ii,:,2)'])
    pause()
end
% all data are in the same scale, it can be concatenate
%%
M = test_cvxformulation_D(y_total,1,30);
M = augment_score(M,size(y_total,2),'llh');
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered','M')
%%
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D')
load('.\experiment\experiment_real_data\AAL_116.mat')
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
M = augment_score(M,size(y_total,2),'sse');
index = M.index.eBIC;
figure(2);
name_list = {'TDC','ADHD'};
for kk=1:2
subplot(1,2,kk)
imagesc(M.model(index).GC(AAL_116.DMN,AAL_116.DMN,kk).*(1-eye(length(AAL_116.DMN))))
grid on
set(gca,'xtick',1:1:length(AAL_116.DMN),'ytick',1:1:length(AAL_116.DMN), ...
    'xticklabel',AAL_116.name(AAL_116.DMN),'yticklabel',AAL_116.name(AAL_116.DMN))
title(name_list{kk})
colormap(1-gray)
xtickangle(30)
caxis([0 0.25])
end
%%
M.fcon(:,:,1) = partialcorr(y_TDC_concat');
M.fcon(:,:,2) = partialcorr(y_ADHD_C_concat');
figure(3)
for kk=1:2
subplot(1,2,kk)
imagesc(abs((M.fcon(AAL_116.DMN,AAL_116.DMN,kk).*(1-eye(length(AAL_116.DMN))))))
grid on
set(gca,'xtick',1:1:length(AAL_116.DMN),'ytick',1:1:length(AAL_116.DMN), ...
    'xticklabel',AAL_116.name(AAL_116.DMN),'yticklabel',AAL_116.name(AAL_116.DMN))
title(name_list{kk})
colormap(1-gray)
xtickangle(30)
caxis([0 0.9])
end
%%
fcon_TDC_avg = zeros(size(y_TDC,1));
fcon_ADHD_C_avg = zeros(size(y_ADHD_C,1));
parfor kk=1:K
    disp(kk)
    fcon_TDC_avg = fcon_TDC_avg+partialcorr(y_TDC(:,:,kk)')/K;
    fcon_ADHD_C_avg = fcon_ADHD_C_avg+partialcorr(y_ADHD_C(:,:,kk)')/K;
end
M.fcon_avg(:,:,1)=fcon_TDC_avg;
M.fcon_avg(:,:,2)=fcon_ADHD_C_avg;
%%
figure(4)
for kk=1:2
subplot(1,2,kk)
imagesc(abs((M.fcon_avg(AAL_116.DMN,AAL_116.DMN,kk).*(1-eye(length(AAL_116.DMN))))))
grid on
set(gca,'xtick',1:1:length(AAL_116.DMN),'ytick',1:1:length(AAL_116.DMN), ...
    'xticklabel',AAL_116.name(AAL_116.DMN),'yticklabel',AAL_116.name(AAL_116.DMN))
title(name_list{kk})
colormap(1-gray)
xtickangle(30)
% caxis([0 0.9])
end
%%
