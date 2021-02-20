clear
clc
inpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';
load([inpath,'formulation_S_summary_model_selection_score.mat'],'DBG');

dd = size(DBG,1);
realz = size(DBG,2);
name_list = {'L','df','bic','aic','aicc'};
for ii=1:dd
    for jj=1:realz
        for kk=1:length(name_list)
            tmp = [DBG(ii,jj).model_selection_score];
            tmp = reshape([tmp.(name_list{kk})],30,30);
            summary.(name_list{kk}){ii,jj} =tmp;
            if kk>2
                [~,index.(name_list{kk})(ii,jj)] = min(tmp(:));
            end
            
        end
    end
end
% save([inpath,'formulation_S_index'],'index')
%% plot index
ii=1;
h.bic = zeros(30,30);
h.aic = zeros(30,30);
h.aicc = zeros(30,30);
figure;
for kk=3:5
    for jj=1:realz
        h.(name_list{kk})(index.(name_list{kk})(ii,jj))=h.(name_list{kk})(index.(name_list{kk})(ii,jj))+1;
    end
    
    subplot(1,3,kk-2)
    imagesc(h.(name_list{kk}))
end
% Conclusion is, BIC should be the best amoung BIC, AIC, AICc


