n=4;
K=3;
CGN = zeros(n,n,K);
figurepath = './plotting/figures/';
h = [4,1;4,2;1,3;2,4];


common_ind = comm_net(h(:,1),h(:,2),n,K);
rng(15648);
CGN(common_ind(:))=rand(size(common_ind))/2+0.25;
% comm_net(h(:,1),h(:,2),n,K)
f1=figure(1);
tt = tiledlayout(1,K, 'Padding','compact', 'TileSpacing','compact');
for kk=1:K
    nexttile;
    imagesc(CGN(:,:,kk))
    set(gca,'xticklabel',[],'yticklabel',[],'FontSize',28)
    title(sprintf('model #%d',kk),'FontSize',40)
    axis('square')
    
%     grid on
    if kk==K
        colorbar;
    end
    colormap(flipud(hot))
end
pp = get(0, 'Screensize');
% pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
saveas(gcf,[figurepath,'example_CGN'])
print([figurepath,'example_CGN'],'-painters','-depsc','-r300')
%%

DGN = CGN;
DGN(2,1,1)=0.5;
DGN(3,2,2)=0.3;
DGN(2,3,3)=0.15;

f2=figure(2);
tt = tiledlayout(1,K, 'Padding','compact', 'TileSpacing','compact');
for kk=1:K
    nexttile;
    imagesc(DGN(:,:,kk))
    
    axis('square')
    set(gca,'xticklabel',[],'yticklabel',[],'FontSize',28)
    title(sprintf('model #%d',kk),'FontSize',40)
%     grid on
    if kk==K
        colorbar;
    end
    colormap(flipud(hot))
end
pp = get(0, 'Screensize');
% pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
saveas(gcf,[figurepath,'example_DGN'])
print([figurepath,'example_DGN'],'-painters','-depsc','-r300')
%%

FGN= CGN;
FGN(:,:,2)=FGN(:,:,1);
FGN(:,:,3)=FGN(:,:,2);
FGN(2,1,1)=0.5;
FGN(3,2,2)=0.3;
FGN(2,3,3)=0.15;
f3=figure(3);
tt = tiledlayout(1,K, 'Padding','compact', 'TileSpacing','compact');
for kk=1:K
    nexttile;
    imagesc(FGN(:,:,kk))
    
    axis('square')
    set(gca,'xticklabel',[],'yticklabel',[],'FontSize',28)
    title(sprintf('model #%d',kk),'FontSize',40)
%     grid on
    if kk==K
        colorbar;
    end
    colormap(flipud(hot))
end
pp = get(0, 'Screensize');
% pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
saveas(gcf,[figurepath,'example_FGN'])
print([figurepath,'example_FGN'],'-painters','-depsc','-r300')

function cind = comm_net(i,j,n,K)
ind = sub2ind([n n],i,j);
cind = [];
for a =1:length(ind)
    cind = [cind ;ind(a):n^2:n^2*K];
end

end
