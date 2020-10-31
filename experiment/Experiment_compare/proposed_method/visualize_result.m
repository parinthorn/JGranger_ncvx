close all
clc
cd = 3;

realz = 7;
type = 2; %D
mname = {'1','5'};
dd=1;
true_model = E{type,cd,dd,realz};
plot_group_GC(true_model.GC);
set(figure(1),'Position',[962 819 720 176])
sgtitle(['true diff density:',mname{dd}])


ind = model_acc(dd,realz).M.index.bic;
estimated_model = model_acc(dd,realz).M.model(ind);
plot_group_GC(estimated_model.GC);
set(figure(2),'Position',[963 585 719 148])
sgtitle(['est diff density:',mname{dd}])
model_acc(dd,realz).result(ind).differential

dd = 2;
true_model = E{type,cd,dd,realz};
plot_group_GC(true_model.GC);
WSIZE = 100;
set(figure(3),'Position',[962 819-500 720 176])
sgtitle(['true diff density:',mname{dd}])


ind = model_acc(dd,realz).M.index.bic;
estimated_model = model_acc(dd,realz).M.model(ind);
plot_group_GC(estimated_model.GC);
set(figure(4),'Position',[963 585-500 719 148])
sgtitle(['est diff density:',mname{dd}])
model_acc(dd,realz).result(ind).differential
%%

close all
clc
GridSize = 30;
cd = 3;
realz = 3;
type = 2; %D
mname = {'1','5'};

figure;
subplot(1,2,1)
dd=1;
tmp = [model_acc(dd,realz).result.total];
F1 = reshape([tmp.F1],[GridSize GridSize]);
max(F1,[],'all')
imagesc(F1);
subplot(1,2,2)
dd = 2;
tmp = [model_acc(dd,realz).result.total];
F1 = reshape([tmp.F1],[GridSize GridSize]);
imagesc(F1)
max(F1,[],'all')
