function s = plot_group_GC(GC)
% size(GC) = [n,n,K]
[n,~,K] = size(GC);
commonNZ = ones(n,n)-eye(n);
diag_ind = 1:n+1:n^2;
for kk=1:K
    tmp = GC(:,:,kk);
    tmp(diag_ind) = 0;
    commonNZ = commonNZ & (tmp~=0);
end
s = figure;
for kk =1:K
        tmp = GC(:,:,kk);
    tmp(diag_ind) = 0;
    subplot(1,K,kk)
    spy(tmp,'r')
    hold on
    spy(commonNZ,'k')
    hold off
    axis('square')
    set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
%     subplot(2,K,i+K)
%     imagesc((GC(i)))
%     axis('square')
%     colormap ((1-gray));
%     set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
end
hold off
end