
name_list = {'df','L','bic','aicc'};
h.df = zeros(30,30);
h.bic = zeros(30,30);
h.aicc = zeros(30,30);
h.L = zeros(30,30);
for ii=1:30
    for jj=1:30
        for kk = 1:length(name_list)
        h.(name_list{kk})(ii,jj) = M.model(ii,jj).stat.model_selection_score.(name_list{kk});%length(find(M.model(ii,jj).A(:)));
        end
        h.df_lasso(ii,jj) = length(find(M.model(ii,jj).A(:)));
    end
end
%%
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';
load([resultpath,'formulation_S_ALL_RESULT_fixdf'],'ALL_RESULT')
load([resultpath,'formulation_S_result_fixdf'],'R')
% summary(ii,jj).F1 = zeros(30,30);
% summary(ii,jj).MCC = zeros(30,30);
% summary(ii,jj).ACC = zeros(30,30);
for ii=1:2
    h = zeros(30,30);
    figure;
    for jj=1:27
        h(R.index(ii,jj).bic) = h(R.index(ii,jj).bic)+1;
%         F1(ii,jj,:,:) = ALL_RESULT(ii,jj).model_acc(R.index(ii,jj).bic).F1;
    end
    imagesc(h)
end