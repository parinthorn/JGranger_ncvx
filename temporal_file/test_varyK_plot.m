clear
clf;close all
type_acc = {'total','common','differential'};
% load('C:/Users/CU_EE_LAB408/Desktop/tmp/cvx_varyK_formulation_D_result_ALL_RESULT.mat')

% load('C:/Users/CU_EE_LAB408/Desktop/tmp/cvx_varyK_formulation_D_result_acc.mat')
% resultpath = 'D:/parinthorn_thesis/formulation_D_result/';
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
load([resultpath,'cvx_varyK_formulation_D_result_ALL_RESULT_50.mat'])
load([resultpath,'cvx_varyK_formulation_D_result_acc_50.mat'])
for tt=1:length(type_acc)
    figure(tt)
K_list = [5,15,25,35,50];
acc_list = {'F1','ACC','MCC'};
plot_seq = reshape(1:15,[3,5]);
plot_seq = plot_seq(:);
cnt = 0;
for ss=1:length(acc_list)
for kk=1:length(K_list)
    cnt=cnt+1;
    val = zeros(30,30);
    list_ff = 1:size(ALL_RESULT,3);
%     list_ff=2:2;
    val_bic = 0;
    for ff=list_ff
    tmp = [ALL_RESULT(1,kk,ff).model_acc.(type_acc{tt})];tmp = [tmp.(acc_list{ss})];val = val +reshape(tmp,30,30)/length(list_ff);
    tmp2 = (R.(type_acc{tt}).(acc_list{ss}));
    if size(tmp2,1)~=1
        tmp2 = squeeze(tmp2(1,:,:));
    end
    val_bic = val_bic + tmp2(kk,ff)/length(list_ff);
    end
    subplot(length(acc_list),length(K_list),plot_seq(cnt))
    imagesc(val);
    
    for ff=list_ff
        bic_index = ALL_RESULT(1,kk,ff).index.bic;
        [col_bic(ff),row_bic(ff)] = ind2sub([30,30],bic_index);
        
        aicc_index = ALL_RESULT(1,kk,ff).index.aicc;
        [col_aicc(ff),row_aicc(ff)] = ind2sub([30,30],aicc_index);
    end

    
    
    title(sprintf('K=%d,bestcase=%.3f,bic=%.3f',K_list(kk),max(max(val)),val_bic))
    axis('square')
    colormap((1-gray).^0.4)
    caxis([0,1])
    set(gca,'xticklabel',[],'yticklabel',[])
    hold on
    scatter(row_bic,col_bic,'r','filled')
    scatter(row_aicc,col_aicc,'b','filled')
    
    hold off
    if kk==1
        ylabel(acc_list{ss})
        legend('BIC','AICc')
    end
end
% sgtitle(acc_list{ss})
end
sgtitle(type_acc{tt})
end
%%
clf;close all
K_list = [5,15,25,35,50];
figure;
for kk=1:length(K_list)
    FPR = zeros(30,30);
    TPR = zeros(30,30);
    for ff=1:5
    tmp = [ALL_RESULT(1,kk,ff).model_acc.total];tmp = [tmp.FPR];FPR = FPR+reshape(tmp,30,30)/5;
    tmp = [ALL_RESULT(1,kk,ff).model_acc.total];tmp = [tmp.TPR];TPR = TPR+reshape(tmp,30,30)/5;
    end
    subplot(3,length(K_list),kk)
    plot(FPR,TPR)
    title(sprintf('K=%d',K_list(kk)))
    axis([0 1 0 1])
    axis('square')
%     caxis([0,1])
    if kk==1
        ylabel('total')
    end
end

for kk=1:length(K_list)
    FPR = zeros(30,30);
    TPR = zeros(30,30);
    for ff=1:5
    tmp = [ALL_RESULT(1,kk,ff).model_acc.common];tmp = [tmp.FPR];FPR = FPR+reshape(tmp,30,30)/5;
    tmp = [ALL_RESULT(1,kk,ff).model_acc.common];tmp = [tmp.TPR];TPR = TPR+reshape(tmp,30,30)/5;
    end
    subplot(3,length(K_list),kk+5)
    plot(FPR,TPR)
    title(sprintf('K=%d',K_list(kk)))
%     axis([0 0.2 0.9 1])
axis([0 1 0 1])
    axis('square')
%     caxis([0,1])
    if kk==1
        ylabel('common')
    end
end

for kk=1:length(K_list)
    if kk==1
    end
    FPR = zeros(30,30);
    TPR = zeros(30,30);
    for ff=1:5
    tmp = [ALL_RESULT(1,kk,ff).model_acc.differential];tmp = [tmp.FPR];FPR = FPR+reshape(tmp,30,30)/5;
    tmp = [ALL_RESULT(1,kk,ff).model_acc.differential];tmp = [tmp.TPR];TPR = TPR+reshape(tmp,30,30)/5;
    end
    subplot(3,length(K_list),kk+10)
    plot(FPR,TPR)
    title(sprintf('K=%d',K_list(kk)))
%     axis([0 0.2 0.9 1])
axis([0 1 0 1])
    axis('square')
%     caxis([0,1])
    if kk==1
        ylabel('differential')
    end
end