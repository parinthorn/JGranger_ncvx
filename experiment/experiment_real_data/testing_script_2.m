clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0);
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0);
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered')

M_in =M;
clear M
clf
close all
% M.TDC = fix_loglikelihood(M_in.TDC,y_TDC);
% M.ADHD_C = fix_loglikelihood(M_in.ADHD_C,y_ADHD_C);
M.TDC = augment_score(M_in.TDC,size(y_TDC,2),'LLH_hetero');
M.ADHD_C = augment_score(M_in.ADHD_C,size(y_TDC,2),'LLH_hetero');
%%
% M = fix_loglikelihood(M,y_total);


M = augment_score(M,size(y_total,2),'LLH_hetero','similar');

tmp = [M.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];

eBIC = reshape([tmp.eBIC],[30,30]);
df = reshape([tmp.df],[30,30]);
n=116;p=1;K=18;T=172;
gamma = log(n^2*p*K)/log(n*(T-p)*K);
kappa = min([1,1.5*(1-1/(2*gamma))]);
binom_term = arrayfun(@(x) log_stirling_approx(n^2*p*K)-log_stirling_approx(n^2*p*K-x)-log_stirling_approx(x) , df);
eBIC_complexity =log(T-p)*df + 2*kappa*binom_term;
LLH_full = reshape([tmp.LLH_full],[30,30]);
LLH_hetero = reshape([tmp.LLH_hetero],[30,30]);
LLH_homo = reshape([tmp.LLH_homo],[30,30]);

figure(99)
subplot(2,2,1)
imagesc(eBIC_complexity)
title('eBIC complexity')
subplot(2,2,2)
imagesc(LLH_full)
title('LLH full')
subplot(2,2,3)
imagesc(LLH_hetero)
title('LLH hetero')
subplot(2,2,4)
imagesc(LLH_homo)
title('LLH homo')
figure(98)
subplot(2,2,1)
imagesc(eBIC_complexity)
title('eBIC complexity')
subplot(2,2,2)
imagesc(LLH_full+eBIC_complexity)
title('eBIC full')
subplot(2,2,3)
imagesc(LLH_hetero+eBIC_complexity)
title('eBIC hetero')
subplot(2,2,4)
imagesc(LLH_homo+eBIC_complexity)
title('eBIC homo')


%%
tmp = [M.TDC.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];