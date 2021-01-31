% [data,header,raw] = tsvread('C:\Users\CU_EE_LAB408\Desktop\sfnwmrda0010001_session_1_rest_1_aal_TCs.1D');
% load('./experiment/experiment_real_data/AAL_116.mat')
% y = data(2:end,3:end)';

%% ADHD formulation C
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu');
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu');
K = size(y_TDC,3);
fcon_TDC = zeros(size(y_TDC,1));
fcon_ADHD_C = zeros(size(y_ADHD_C,1));
for kk=1:K
    disp(kk)
    fcon_TDC = fcon_TDC+partialcorr(y_TDC(:,:,kk)')/K;
    fcon_ADHD_C = fcon_ADHD_C+partialcorr(y_ADHD_C(:,:,kk)')/K;
end
%%
y_TDC_concat = reshape(y_TDC,[116,172*18]);
y_ADHD_C_concat = reshape(y_ADHD_C,[116,172*18]);
fcon_TDC_concat = partialcorr(y_TDC_concat');
fcon_ADHD_C_concat = partialcorr(y_ADHD_C_concat');

%%
figure(1)
subplot(2,2,1)
imagesc(fcon_TDC.*(abs(fcon_TDC.*(1-eye(116)))>0.2))
axis('square')
subplot(2,2,2)
imagesc(fcon_ADHD_C.*(abs(fcon_ADHD_C.*(1-eye(116)))>0.2))
axis('square')

subplot(2,2,3)
imagesc(fcon_TDC_concat.*(abs(fcon_TDC_concat.*(1-eye(116)))>0.4))
axis('square')
subplot(2,2,4)
imagesc(fcon_ADHD_C_concat.*(abs(fcon_ADHD_C_concat.*(1-eye(116)))>0.4))
axis('square')
%%
y_total(:,:,1) = y_TDC_concat;
y_total(:,:,2) = y_ADHD_C_concat;

M_total = test_cvxformulation_D(10*y_total,1,30);
%%

M_TDC = test_cvxformulation_C(10*y_TDC,1,30);
M_ADHD_C = test_cvxformulation_C(10*y_ADHD_C,1,30);