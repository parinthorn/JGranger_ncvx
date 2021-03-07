clear
clc
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
[~,E.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)'));
[~,E.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));

E.TDC=E.TDC';
E.ADHD=E.ADHD';
score_diff = E.TDC-E.ADHD;

[C,I] = sort(score_diff(:),'descend');
fprintf('D2K\n')
for ii=1:5
    selected_index = I(ii);
    [i,j] = ind2sub([116,116],selected_index);
    fprintf([AAL_116.name_full{j},'->',AAL_116.name_full{i},'\n'])
end


%%
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
[~,E.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)'));
[~,E.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));

E.TDC=E.TDC';
E.ADHD=E.ADHD';
score_diff = E.TDC-E.ADHD;

[C,I] = sort(score_diff(:),'descend');
fprintf('S2K\n')
for ii=1:5
    selected_index = I(ii);
    [i,j] = ind2sub([116,116],selected_index);
    fprintf([AAL_116.name_full{j},'->',AAL_116.name_full{i},'\n'])
end
%%
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered_LLHcorrection')
GC.total(:,:,1) = mean(M.TDC.model(M.TDC.index.eBIC).GC,3);
GC.total(:,:,2) = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3);
[~,E.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)'));
[~,E.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));

E.TDC=E.TDC';
E.ADHD=E.ADHD';
score_diff = E.TDC-E.ADHD;

[C,I] = sort(score_diff(:),'descend');
fprintf('C18K\n')
for ii=1:5
    selected_index = I(ii);
    [i,j] = ind2sub([116,116],selected_index);
    fprintf([AAL_116.name_full{j},'->',AAL_116.name_full{i},'\n'])
end