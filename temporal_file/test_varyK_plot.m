clf;close all
type_acc = {'total','common','differential'};
load('C:/Users/CU_EE_LAB408/Desktop/tmp/varyK_formulation_D_result_ALL_RESULT.mat')

load('C:/Users/CU_EE_LAB408/Desktop/tmp/varyK_formulation_D_resultall.mat')
VK = ALL_RESULT(1,:,2);
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
%     list_ff = 1:size(ALL_RESULT,3);
    list_ff=5;
    val_bic = 0;
    for ff=list_ff
    tmp = [ALL_RESULT(1,kk,ff).model_acc.(type_acc{tt})];tmp = [tmp.(acc_list{ss})];val = val +reshape(tmp,30,30)/length(list_ff);
    tmp2 = squeeze(R.(type_acc{tt}).(acc_list{ss}));
    val_bic = val_bic + tmp2(kk,ff)/length(list_ff);
    end
    subplot(length(acc_list),length(K_list),plot_seq(cnt))
    imagesc(val);
    title(sprintf('K=%d,bestcase=%.3f,bic=%.3f',K_list(kk),max(max(val)),val_bic))
    axis('square')
    colormap((1-gray).^0.4)
    caxis([0,1])
    set(gca,'xticklabel',[],'yticklabel',[])

    if kk==1
        ylabel(acc_list{ss})
    end
end
% sgtitle(acc_list{ss})
end
sgtitle(type_acc{tt})
end
%%
clf;close all
K_list = [5,15,25,35,50];
for kk=1:length(K_list)
    FPR = zeros(30,30);
    TPR = zeros(30,30);
    for ff=1:5
    tmp = [ALL_RESULT(1,kk,ff).model_acc.total];tmp = [tmp.FPR];FPR = FPR+reshape(tmp,30,30)/5;
    tmp = [ALL_RESULT(1,kk,ff).model_acc.total];tmp = [tmp.TPR];TPR = TPR+reshape(tmp,30,30)/5;
    end
    subplot(1,length(K_list),kk)
    plot(FPR,TPR)
    title(sprintf('K=%d',K_list(kk)))
    axis([0 0.2 0.9 1])
    axis('square')
%     caxis([0,1])
end