function s = plot_group_GC(GC)
% size(GC) = [n,n,K]
[n,~,K] = size(GC);
commonNZ = ones(n,n);
for kk=1:K
    commonNZ = commonNZ & (GC(:,:,kk)~=0);
end
s = figure;
for kk =1:K
    subplot(1,K,kk)
    spy(GC(:,:,kk),'r')
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