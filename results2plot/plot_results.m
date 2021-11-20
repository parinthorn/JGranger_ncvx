%% exp_CGN
clear
clc
clf
close all
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
table_head = {'CGN','cvxCGN','Song','Greg'};
table_head_show = {'CGN','cvx-CGN','Song17C','Greg15'};
row_name = {'F1','FPR','TPR','ACC','MCC', 'bias'};
load([resource_path,'CGN_result'])
result.CGN = R;
load([resource_path,'CGN_CVX_result'])
result.cvxCGN = R;
% load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\adaptive_formulation_C_JSS_result.mat'])
load([resource_path,'CGN_JSS_result'])
result.Song = R;
load([resource_path,'magda_result'])
result.Greg = R;
M = zeros(5,4,2);
STD = zeros(5,4,2);
MEDIAN = zeros(5,4,2);
summary = zeros(5,4,2,100);
for jj=1:length(table_head)
    for ii=1:length(row_name)
        if ii==6
            result.(table_head{jj}).common.bias = result.(table_head{jj}).bias;
        end
        
        for dd=1:2
            if  jj~=4
                M(ii,jj,dd) = mean(result.(table_head{jj}).common.(row_name{ii})(dd,:));
                STD(ii,jj,dd) = std(result.(table_head{jj}).common.(row_name{ii})(dd,:));
                MEDIAN(ii,jj,dd) = median(result.(table_head{jj}).common.(row_name{ii})(dd,:));
                summary(ii,jj,dd,:) = result.(table_head{jj}).common.(row_name{ii})(dd,:);
            else
                M(ii,jj,dd) = mean(result.(table_head{jj}).(row_name{ii})(dd,:));
                STD(ii,jj,dd) = std(result.(table_head{jj}).(row_name{ii})(dd,:));
                MEDIAN(ii,jj,dd) = median(result.(table_head{jj}).(row_name{ii})(dd,:));
                summary(ii,jj,dd,:) = result.(table_head{jj}).(row_name{ii})(dd,:);
            end
        end
    end
end
O(:,:,1) = -M(:,:,1)+M(:,1,1);
O(:,:,2) = -M(:,:,2)+M(:,1,2);
O(2,:,:) = -O(2,:,:);

P(:,:,1) = -MEDIAN(:,:,1)+MEDIAN(:,1,1);
P(:,:,2) = -MEDIAN(:,:,2)+MEDIAN(:,1,2);
P(2,:,:) = -P(2,:,:);
% table
% printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name)
printtable_withtoprow([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name,{'Common density: 10%','Common density: 20%'})

hh = tiledlayout(2,2);
hh.TileSpacing = 'none';
hh.Padding = 'none';
density_name = {'Common density: 10%','Common density: 20%'};
for ii=1:2
    
    
    
    
    for dd=1:2 % this should be
        clear ARR
        ax(dd) = nexttile;
        for jj=1:4
            ARR(:,jj) = summary(ii,jj,dd,:);
        end
        boxplot(100*ARR)
        
        if ii==2
            a = get(get(gca,'children'),'children');   % Get the handles of all the objects
            t = get(a,'tag');   % List the names of all the objects
            box1 = a(9:12);   % The 7th object is the first box
            set(box1, 'Color', [0 0.5 0]);   % Set the color of the first box to green
        end
        
        if ii==1
            title(density_name{dd})
            set(gca,'xticklabel',[])
            
        else
            set(gca,'xticklabel',table_head_show)
            ylim([0 55])
        end
        grid on
        if dd==1
            ylabel([row_name{ii},' (%)'])
        else
            
            set(gca,'yticklabel',[])
            linkaxes([ax(1) ax(2)],'y')
        end
        set(findobj(gca,'type','line'),'linew',3)
    end
    
    
end
set(findall(gcf,'-property','FontSize'),'FontSize',28)

hh.TileSpacing = 'none';
hh.Padding = 'none';
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
%     saveas(gcf,[figurepath,'exp_CGN'])
% print([figurepath,'exp_CGN'],'-painters','-depsc','-r300')
% print([figurepath,'exp_CGN_sample'],'-painters','-dpng','-r300')
%% exp_DGN_A
clear
clc
clf
close all
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
table_head = {'DGN','cvxDGN','Song','Skrip'};
table_head_show = {'DGN','cvx-DGN','Song17D','Skrip19b'};
row_name = {'F1','FPR','TPR','ACC','MCC','bias'};
load([resource_path,'DGN_result_K5'])
result.DGN = R;
load([resource_path,'DGN_CVX_result_K5'])
result.cvxDGN = R;
load([resource_path,'DGN_JSS_result_K5'])
result.Song = R;
tmp = load([resource_path,'ResultSkripD_K5']);
R = tmp.R;
clear tmp
result.Skrip = R;
M = zeros(5,4,2);
MEDIAN = zeros(5,4,2);
STD = zeros(5,4,2);
summary = zeros(5,4,2,100);
for jj=1:length(table_head)
    for ii=1:length(row_name)
        if ii==6
            result.(table_head{jj}).total.bias = result.(table_head{jj}).bias;
        end
        for dd=1:2
            M(ii,jj,dd) = mean(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            STD(ii,jj,dd) = std(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            MEDIAN(ii,jj,dd) = median(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            summary(ii,jj,dd,:) = result.(table_head{jj}).total.(row_name{ii})(dd,:);
            
        end
    end
end
O(:,:,1) = -M(:,:,1)+M(:,1,1);
O(:,:,2) = -M(:,:,2)+M(:,1,2);
O(2,:,:) = -O(2,:,:);

P(:,:,1) = -MEDIAN(:,:,1)+MEDIAN(:,1,1);
P(:,:,2) = -MEDIAN(:,:,2)+MEDIAN(:,1,2);
P(2,:,:) = -P(2,:,:);
% table
% printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name)
printtable_withtoprow([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name,{'Diff. density: 1%', 'Diff. density: 5%'})
hh = tiledlayout(2,2);
hh.TileSpacing = 'none';
hh.Padding = 'none';
density_name = {'Differential density: 1%','Differential density: 5%'};
for ii=1:2
    
    
    
    
    for dd=1:2 % this should be
        clear ARR
        ax(dd) = nexttile;
        for jj=1:4
            ARR(:,jj) = summary(ii,jj,dd,:);
        end
        boxplot(100*ARR)
        if ii==2
            a = get(get(gca,'children'),'children');   % Get the handles of all the objects
            t = get(a,'tag');   % List the names of all the objects
            box1 = a(9:12);   % The 7th object is the first box
            set(box1, 'Color', [0 0.5 0]);   % Set the color of the first box to green
        end
        
        if ii==1
            title(density_name{dd})
            set(gca,'xticklabel',[])
            
        else
            set(gca,'xticklabel',table_head_show)
            %                 ylim([0 55])
        end
        grid on
        if dd==1
            ylabel([row_name{ii},' (%)'])
        else
            
            set(gca,'yticklabel',[])
            linkaxes([ax(1) ax(2)],'y')
        end
        set(findobj(gca,'type','line'),'linew',3)
    end
    
    
end
set(findall(gcf,'-property','FontSize'),'FontSize',28)
pp = get(0, 'Screensize');
% pp(2) = pp(4);
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
set(findall(gcf,'-property','FontSize'),'FontSize',28)
% saveas(gcf,[figurepath,'exp_DGN_A'])
% print([figurepath,'exp_DGN_A'],'-painters','-depsc','-r300')
%% exp_DGN_B_xxx
clear
clc
clf
close all
toggle = 'total';
dd=2;
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
table_head = {'DGN','cvxDGN','Song','Skrip'};
table_head_show = {'DGN','cvx-DGN','Song17D','Skrip19b'};
row_name = {'F1','FPR','TPR','ACC','MCC', 'bias'};
K_list = {'K5','K50'};

load([resource_path,'DGN_result_K5'])
result.DGN.K5 = R;
load([resource_path,'DGN_CVX_result_K5'])
result.cvxDGN.K5 = R;
load([resource_path,'DGN_JSS_result_K5'])
result.Song.K5 = R;
tmp=load([resource_path,'ResultSkripD_K5']);
R = tmp.R;
clear tmp
result.Skrip.K5 = R;


load([resource_path,'DGN_result_K50'])
result.DGN.K50 = R;
load([resource_path,'DGN_CVX_result_K50'])
result.cvxDGN.K50 = R;
load([resource_path,'DGN_JSS_result_K50'])
result.Song.K50 = R;
tmp = load([resource_path,'ResultSkripD_K50']);
R = tmp.R;
clear tmp
result.Skrip.K50 = R;
M = zeros(6,4,2);
STD = zeros(6,4,2);
MEDIAN = zeros(6,4,2);
summary = zeros(6,4,2,100);

for jj=1:length(table_head)
    
    for ii=1:length(row_name)
        
        for kk=1:length(K_list)
            if ii==6
                result.(table_head{jj}).(K_list{kk}).(toggle).bias = result.(table_head{jj}).(K_list{kk}).bias;
            end
            dd=2;
            M(ii,jj,kk) = mean(result.(table_head{jj}).(K_list{kk}).(toggle).(row_name{ii})(dd,:));
            STD(ii,jj,kk) = std(result.(table_head{jj}).(K_list{kk}).(toggle).(row_name{ii})(dd,:));
            MEDIAN(ii,jj,kk) = median(result.(table_head{jj}).(K_list{kk}).(toggle).(row_name{ii})(dd,:));
            summary(ii,jj,kk,:) = result.(table_head{jj}).(K_list{kk}).(toggle).(row_name{ii})(dd,:);
            
        end
    end
end
O(:,:,1) = -M(:,:,1)+M(:,1,1);
O(:,:,2) = -M(:,:,2)+M(:,1,2);
O(2,:,:) = -O(2,:,:);

P(:,:,1) = -MEDIAN(:,:,1)+MEDIAN(:,1,1);
P(:,:,2) = -MEDIAN(:,:,2)+MEDIAN(:,1,2);
P(2,:,:) = -P(2,:,:);
H = table([O(1:2,2:end,1);O(1:2,2:end,2);P(1:2,2:end,1);P(1:2,2:end,2)]);

M = permute(M,[1,3,2]);
STD = permute(STD,[1,3,2]);

Mt = reshape(M,[6,8]);
STDt = reshape(STD,[6,8]);
M2 = Mt(:,[1,5,2,6,3,7,4,8]);
STD2 =  STDt(:,[1,5,2,6,3,7,4,8]);
% table
% tmp = {table_head_show{:};table_head_show{:}};
% tmp = {tmp{:}};
tmp = {'$5$','$50$','$5$','$50$','$5$','$50$','$5$','$50$'};
% printtable(Mt*100,STDt*100,tmp,row_name)
printtable_withtoprow(Mt*100,STDt*100,tmp,row_name,table_head_show,'$K$')
hh = tiledlayout(2,4);
hh.TileSpacing = 'none';
hh.Padding = 'none';
density = {'1%','5%'};
for ii=1:2
    
    for jj=1:4
        ax(jj) = nexttile;
        clear ARR
        for kk=1:2
            ARR(:,kk) = summary(ii,jj,kk,:);
        end
        boxplot(100*ARR)
        if ii==2
            a = get(get(gca,'children'),'children');   % Get the handles of all the objects
            t = get(a,'tag');   % List the names of all the objects
            box1 = a(5:6);   % The 7th object is the first box
            set(box1, 'Color', [0 0.5 0]);   % Set the color of the first box to green
        end
        grid on
        if jj==1
            ylabel(([row_name{ii},' (%)']))
        else
            set(gca,'yticklabel',[])
        end
        if ii==1
            set(gca,'xticklabel',[])
            title(table_head_show{jj})
        else
            
            set(gca,'xticklabel',{'5','50'})
        end
        set(findobj(gca,'type','line'),'linew',3)
    end
    linkaxes([ax(1) ax(2) ax(3) ax(4)],'y')
end
xlabel(hh, '$K$','Interpreter','Latex','FontSize',28)
set(findall(gcf,'-property','FontSize'),'FontSize',28)
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
% saveas(gcf,[figurepath,'exp_DGN_B_',toggle])
% print([figurepath,'exp_DGN_B_',toggle],'-painters','-depsc','-r300')
%% exp_FGN

clear
clc
clf
close all
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
table_head = {'FGN','cvxFGN','Song','Skrip'};
table_head_show = {'FGN','cvx-FGN','Song15','Skrip19a'};
row_name = {'F1','FPR','TPR','ACC','MCC', 'bias'};
load([resource_path,'FGN_result_K5'])
result.FGN = R;
load([resource_path,'FGN_CVX_result_K5'])
result.cvxFGN = R;
load([resource_path,'FGN_JSS_result_K5'])
result.Song = R;
load([resource_path,'ResultSkripS_K5'])
result.Skrip = R;
M = zeros(5,4,2);
MEDIAN = zeros(5,4,2);
STD = zeros(5,4,2);
summary = zeros(5,4,2,100);
for jj=1:length(table_head)
    for ii=1:length(row_name)
        if ii==6
            result.(table_head{jj}).total.bias = result.(table_head{jj}).bias;
        end
        for dd=1:2
            M(ii,jj,dd) = mean(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            STD(ii,jj,dd) = std(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            MEDIAN(ii,jj,dd) = median(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            summary(ii,jj,dd,:) = result.(table_head{jj}).total.(row_name{ii})(dd,:);
            
        end
    end
end
O(:,:,1) = -M(:,:,1)+M(:,1,1);
O(:,:,2) = -M(:,:,2)+M(:,1,2);
O(2,:,:) = -O(2,:,:);

P(:,:,1) = -MEDIAN(:,:,1)+MEDIAN(:,1,1);
P(:,:,2) = -MEDIAN(:,:,2)+MEDIAN(:,1,2);
P(2,:,:) = -P(2,:,:);
% table
printtable_withtoprow([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name,{'Diff. density: 1%', 'Diff. density: 5%'})

hh = tiledlayout(2,2);
hh.TileSpacing = 'none';
hh.Padding = 'none';
density_name = {'Differential density: 1%','Differential density: 5%'};
for ii=1:2
    
    
    
    
    for dd=1:2 % this should be
        clear ARR
        ax(dd) = nexttile;
        for jj=1:4
            ARR(:,jj) = summary(ii,jj,dd,:);
        end
        boxplot(100*ARR)
        if ii==2
            a = get(get(gca,'children'),'children');   % Get the handles of all the objects
            t = get(a,'tag');   % List the names of all the objects
            box1 = a(9:12);   % The 7th object is the first box
            set(box1, 'Color', [0 0.5 0]);   % Set the color of the first box to green
        end
        
        if ii==1
            title(density_name{dd})
            set(gca,'xticklabel',[])
            
        else
            set(gca,'xticklabel',table_head_show)
            %                 ylim([0 55])
        end
        grid on
        if dd==1
            ylabel([row_name{ii},' (%)'])
        else
            
            set(gca,'yticklabel',[])
            linkaxes([ax(1) ax(2)],'y')
        end
        set(findobj(gca,'type','line'),'linew',3)
    end
    
    
end
set(findall(gcf,'-property','FontSize'),'FontSize',28)
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
% saveas(gcf,[figurepath,'exp_FGN'])
% print([figurepath,'exp_FGN'],'-painters','-depsc','-r300')
%% exp_cvxcompare
clear
clc
clf
close all
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
table_head = {'C','D','F'};

row_name = {'F1','FPR','TPR','ACC','MCC', 'bias'};

load([resource_path,'T150_CGN_result_K5'])
result.C.ncvx = R;
load([resource_path,'T150_CGN_CVX_result_K5'])
result.C.cvx = R;

load([resource_path,'T150_DGN_result_K5'])
result.D.ncvx = R;
load([resource_path,'T150_DGN_CVX_result_K5'])
result.D.cvx = R;

load([resource_path,'T150_FGN_result_K5'])
result.F.ncvx = R;
load([resource_path,'T150_FGN_CVX_result_K5'])
result.F.cvx = R;

M = zeros(3,2,5); % [C,D,F] x [ncvx, cvx] x [F1, FPR, TPR, ACC , MCC]
STD = zeros(3,2,5);
MEDIAN = zeros(3,2,5);
summary = zeros(3,2,5,100);
type_acc = {'ncvx','cvx'};
for jj=1:length(table_head)
    if jj==1
        toggle= 'common';
    else
        toggle = 'total';
    end
    for ii=1:length(row_name)
        
        for t=1:2
            if ii==6
                result.(table_head{jj}).(type_acc{t}).(toggle).bias = result.(table_head{jj}).(type_acc{t}).bias;
            end
            dd=2;
            
            M(jj,t,ii) = mean(result.(table_head{jj}).(type_acc{t}).(toggle).(row_name{ii})(dd,:));
            STD(jj,t,ii) = std(result.(table_head{jj}).(type_acc{t}).(toggle).(row_name{ii})(dd,:));
            MEDIAN(jj,t,ii) = median(result.(table_head{jj}).(type_acc{t}).(toggle).(row_name{ii})(dd,:));
            summary(jj,t,ii,:) = result.(table_head{jj}).(type_acc{t}).(toggle).(row_name{ii})(dd,:);
        end
    end
end
M = permute(M,[3,2,1]);
STD = permute(STD,[3,2,1]);
MEDIAN = permute(MEDIAN,[3,2,1]);

O(:,:,1) = -M(:,:,1)+M(:,1,1);
O(:,:,2) = -M(:,:,2)+M(:,1,2);
O(:,:,3) = -M(:,:,2)+M(:,1,3);
O(2,:,:) = -O(2,:,:);

P(:,:,1) = -MEDIAN(:,:,1)+MEDIAN(:,1,1);
P(:,:,2) = -MEDIAN(:,:,2)+MEDIAN(:,1,2);
P(:,:,3) = -MEDIAN(:,:,3)+MEDIAN(:,1,3);
P(2,:,:) = -P(2,:,:);

% table
table_head_show = {'CGN','cvx-CGN','DGN','cvx-DGN','FGN','cvx-FGN'};
tmp_M =100*[M(:,1,1) M(:,2,1) M(:,1,2) M(:,2,2) M(:,1,3) M(:,2,3)];
tmp_STD =100*[STD(:,1,1) STD(:,2,1) STD(:,1,2) STD(:,2,2) STD(:,1,3) STD(:,2,3)];
printtable(tmp_M,tmp_STD,table_head_show,row_name)

hh = tiledlayout(2,3);
hh.TileSpacing = 'none';
hh.Padding = 'none';
density = {'1%','5%'};
formulation_name = {'common part','total part','total part'};
cvx_list = {'CGN','cvx-CGN';'DGN','cvx-DGN';'FGN','cvx-FGN'};
dd=2;
% [C,D,F] x [ncvx, cvx] x [F1, FPR, TPR, ACC , MCC]
for ii=1:2
    for jj=1:3
        ax(jj) = nexttile;
        clear ARR
        for tt=1:2
            ARR(:,tt) = summary(jj,tt,ii,:) ;
        end
        boxplot(100*ARR)
        if ii==2
            a = get(get(gca,'children'),'children');   % Get the handles of all the objects
            t = get(a,'tag');   % List the names of all the objects
            box1 = a(5:6);   % The 7th object is the first box
            set(box1, 'Color', [0 0.5 0]);   % Set the color of the first box to green
        end
        
        
        if ii==1
            title(formulation_name{jj})
            set(gca,'xticklabel',[])
            
        else
            set(gca,'xticklabel',{cvx_list{jj,:}})
        end
        grid on
        
        if jj==1
            ylabel([row_name{ii},' (%)'])
        else
            
            set(gca,'yticklabel',[])
            
        end
        set(findobj(gca,'type','line'),'linew',3)
        
        
    end
    linkaxes([ax(1) ax(2) ax(3)],'y')
    
end
set(findall(gcf,'-property','FontSize'),'FontSize',28)
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
% saveas(gcf,[figurepath,'exp_cvxcompare'])
% print([figurepath,'exp_cvxcompare'],'-painters','-depsc','-r300')
%% exp_CGN_ROC

clear
clc
clf
close all
type_list = {'total','common','differential'};
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
method_path = {'CGN_ALL_RESULT.mat','CGN_CVX_ALL_RESULT.mat'};
for mm=1:2
    load([resource_path,method_path{mm}])
    ii =2;sample = 6;
    sample_acc = ALL_RESULT(1,sample);
    metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmpval(1,:)=[metrics.F1];max1val = {max(tmpval(1,:))};
    
    sample_acc = ALL_RESULT(2,sample);
    metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmpval(2,:)=[metrics.F1];max2val = {max(tmpval(2,:))};
    % tt=tiledlayout(1,2,'padding','none','tilespacing','none');
    
    tmp_avg = zeros(2,30);
    max1_avg = 0;
    max2_avg = 0;
    for jj=1:100
        sample_acc = ALL_RESULT(1,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(1,:)=[metrics.FPR];max1 = {max(tmp(1,:))};
        
        sample_acc = ALL_RESULT(2,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(2,:)=[metrics.FPR];max2 = {max(tmp(2,:))};
        
        tmp_avg = tmp_avg+tmp/100;
        max1_avg = max1_avg+max1{1}/100;
        max2_avg = max2_avg+max2{1}/100;
    end
    FPR_array{mm} = tmp_avg;
    
    tmp_avg = zeros(2,30);
    max1_avg = 0;
    max2_avg = 0;
    for jj=1:100
        sample_acc = ALL_RESULT(1,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(1,:)=[metrics.TPR];max1 = {max(tmp(1,:))};
        
        sample_acc = ALL_RESULT(2,jj);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; tmp(2,:)=[metrics.TPR];max2 = {max(tmp(2,:))};
        
        tmp_avg = tmp_avg+tmp/100;
        max1_avg = max1_avg+max1{1}/100;
        max2_avg = max2_avg+max2{1}/100;
    end
    TPR_array{mm} = tmp_avg;
end
plot(100*[FPR_array{1};FPR_array{2}]',100*[TPR_array{1};TPR_array{2}]','Linewidth',4)
legend('CGN on common density 10%','CGN on common density 20%','cvx-CGN on common density 10%','cvx-CGN on common density 20%','location','southeast')
% u = gca;
% color_list = {[0,0,1],[0,1,0],[0,0,0],[0,0,0]};
% for kk=1:4
% u.Children(kk).Color = color_list{kk};
% end
% legend
axis([0 30 70 100])
axis('square')
% legend('Common density 10%','Common density 20%','location','northeast')
% title(sprintf('10%%, F1 bestcase: %2.1f, 20%%,  F1 bestcase: %2.1f',100*max1_avg{1},100*max2_avg{1}))
% ylabel('F1 score')
% xlabel('denser $\leftarrow\lambda\rightarrow$ sparser','Interpreter','latex')
xlabel('FPR (%)')
ylabel('TPR (%)')
grid on
set(gca,'FontSize',28)
pp = get(0, 'Screensize');
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
% saveas(gcf,[figurepath,'exp_CGN_ROC'])
% print([figurepath,'exp_CGN_ROC'],'-painters','-depsc','-r300')
%% exp_DGN_supp K5 (cvx)

clear
clc
clf
close all
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
load([resource_path,'DGN_CVX_ALL_RESULT_K5.mat'])
sample = 74;
tt=tiledlayout(1,2);
type_list = {'common','differential'};
type_list_show = {'common part','differential part'};
score_type = 'F1';
text_label = {'Diff. density 1%','Diff. density 5%'};
for dd=2:2
    sample_acc = ALL_RESULT(dd,sample);
    
    
    
    for ii=1:length(type_list)
        nexttile;
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; metrics=[metrics.F1];
        bestcase = max(metrics);
        metrics = reshape(metrics,[30,30]);
        imagesc(100*metrics)
        %     grid on
        %     text_show = sprintf([type_list_show{ii},' best ',score_type,' score: %2.1f'],100*bestcase);
        text_show = type_list_show{ii};
        title(text_show)
        axis('square')
        colormap((1-gray).^0.4)
        caxis([0,100])
        set(gca,'xticklabel',[],'yticklabel',[])
        if ii==2
            %         rr= colorbar('location','eastoutside');
            %         set(get(rr,'label'),'string','F1 (%)');
        end
        if ii==1
            ylabel('$\leftarrow\lambda_{1}$','Interpreter','latex')
        end
        xlabel('$\lambda_{2}\rightarrow$','Interpreter','latex')
        
        set(gca,'FontSize',36)
        
        %     set(gca,)
    end
    
end
rr=colorbar('FontSize',36);
rr.Title.String='F1 (%)';
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'exp_DGN_supp_K5'])
% print([figurepath,'exp_DGN_supp_K5'],'-painters','-depsc','-r300')
%% exp_DGN_supp K50 (cvx)

clear
clc
clf
close all
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
load([resource_path,'DGN_CVX_ALL_RESULT_K50.mat'])
sample = 74;
tt=tiledlayout(1,2);
type_list = {'common','differential'};
type_list_show = {'common part','differential part'};
score_type = 'F1';
text_label = {'Diff. density 1%','Diff. density 5%'};
for dd=2:2
    sample_acc = ALL_RESULT(dd,sample);
    
    
    
    for ii=1:length(type_list)
        nexttile;
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; metrics=[metrics.F1];
        bestcase = max(metrics);
        metrics = reshape(metrics,[30,30]);
        imagesc(100*metrics)
        %     grid on
        %     text_show = sprintf([type_list_show{ii},' best ',score_type,' score: %2.1f'],100*bestcase);
        text_show = type_list_show{ii};
        title(text_show)
        axis('square')
        colormap((1-gray).^0.4)
        caxis([0,100])
        set(gca,'xticklabel',[],'yticklabel',[])
        if ii==2
            %         rr= colorbar('location','eastoutside');
            %         set(get(rr,'label'),'string','F1 (%)');
        end
        if ii==1
            ylabel('$\leftarrow\lambda_{1}$','Interpreter','latex')
        end
        xlabel('$\lambda_{2}\rightarrow$','Interpreter','latex')
        
        set(gca,'FontSize',36)
        
        %     set(gca,)
    end
    
end
rr=colorbar('FontSize',36);
rr.Title.String='F1 (%)';
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'exp_DGN_supp_K50'])
% print([figurepath,'exp_DGN_supp_K50'],'-painters','-depsc','-r300')

%% exp_FGN_supp (cvx)

clear
clc
clf
close all
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
load([resource_path,'FGN_CVX_ALL_RESULT_K5.mat'])
load([resource_path,'FGN_CVX_result_K5.mat'])
sample = [1:100];%72
selected_index = [1:100];
tt=tiledlayout(1,2,'padding','compact','tilespacing','none');
type_list = {'total','common','differential'};
type_list_show = {'cvx-FGN, total bestcase'};
score_type = 'F1';
text_label = {'Differential  density 1%','Differential density 5%'};
for dd=1:2
    ii =1;
    GridF1_avg = zeros(30,30);
    GridF1 = zeros(30,30,length(sample));
    I_best = zeros(length(sample),1);
    J_best = zeros(length(sample),1);
    I_score = zeros(length(sample),1);
    J_score = zeros(length(sample),1);
    F1_avg = 0;
    nexttile;
    for sample_ind=1:length(sample)
        sample_acc = ALL_RESULT(dd,sample(sample_ind));
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; metrics=[metrics.(score_type)];
        
        [bestcase,best_index] = max(metrics);
        score_index = R.index(dd,sample_ind).eBIC;
        
        [I_best(sample_ind),J_best(sample_ind)] = ind2sub([30,30],best_index);
        [I_score(sample_ind),J_score(sample_ind)] = ind2sub([30,30],score_index);
        
        metrics = reshape(metrics,[30,30]);
        GridF1(:,:,sample_ind) = metrics;
        GridF1_avg=GridF1_avg+metrics/length(sample);
        F1_avg=F1_avg+bestcase/length(sample);
    end
    distance{dd} = [I_best J_best]-[I_score J_score];
    bestval(dd) = F1_avg;
    imagesc(100*mean(GridF1(:,:,selected_index),3))
    hold on
    scatter(J_best(selected_index),I_best(selected_index),300,'sr','filled')
    scatter(J_score(selected_index),I_score(selected_index),150,'ob','filled')
    %     plot([J_best(selected_index) J_score(selected_index)]',[I_best(selected_index) I_score(selected_index)]','-k','LineWidth',4)
    hold off
    if dd==2
        legend('maximum F1 index','eBIC selected index','FontSize',28)
    end
    %     grid on
    text_show = sprintf([type_list_show{1},' ',score_type,' score: %2.1f'],100*F1_avg);
    %     title(text_show)
    axis('square')
    colormap((1-gray).^0.4)
    caxis([0,100])
    set(gca,'xticklabel',[],'yticklabel',[])
    if dd==1
        %                 colorbar('Location','westoutside')
        xlabel('$\lambda_{2}\rightarrow$','Interpreter','latex','FontSize',38)
        
        set(gca,'xaxisLocation','top')
        ylabel('$\leftarrow \lambda_{1}$','Interpreter','latex','FontSize',38)
        
    else
        xlabel('$\lambda_{2}\rightarrow$','Interpreter','latex','FontSize',38)
        
        set(gca,'xaxisLocation','top')
    end
    
    if ii==1
        title([text_label{dd},sprintf('\n avg best F1: %.2f %%',100*F1_avg)],'FontSize',38)
        
    end
end

rr=colorbar('FontSize',18);
set(get(rr,'label'),'string','F1 (%)');
% rr.Position(4)=rr.Position(4)-0.005;
% rr.Position(2)=rr.Position(2)+0.005/2;
pp = get(0, 'Screensize');
pp(3) = pp(3)*0.8;
set(gcf, 'Position', pp);
% saveas(gcf,[figurepath,'exp_FGN_supp_cvx'])
% print([figurepath,'exp_FGN_supp_cvx'],'-painters','-depsc','-r300')

%% exp_FGN_supp (non-cvx)

clear
clc
clf
close all
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
load([resource_path,'FGN_ALL_RESULT_K5.mat'])
load([resource_path,'FGN_result_K5.mat'])
sample = [1:100];%72
selected_index = [1:1:100];
tt=tiledlayout(1,2,'padding','compact','tilespacing','none');
type_list = {'total','common','differential'};
type_list_show = {'FGN, total bestcase'};
score_type = 'F1';
text_label = {'Differential  density 1%','Differential density 5%'};
for dd=1:2
    ii =1;
    GridF1_avg = zeros(30,30);
    GridF1 = zeros(30,30,length(sample));
    I_best = zeros(length(sample),1);
    J_best = zeros(length(sample),1);
    I_score = zeros(length(sample),1);
    J_score = zeros(length(sample),1);
    F1_avg = 0;
    nexttile;
    for sample_ind=1:length(sample)
        sample_acc = ALL_RESULT(dd,sample(sample_ind));
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; metrics=[metrics.F1];
        
        [bestcase,best_index] = max(metrics);
        score_index = R.index(dd,sample_ind).eBIC;
        
        [I_best(sample_ind),J_best(sample_ind)] = ind2sub([30,30],best_index);
        [I_score(sample_ind),J_score(sample_ind)] = ind2sub([30,30],score_index);
        
        metrics = reshape(metrics,[30,30]);
        GridF1(:,:,sample_ind) = metrics;
        GridF1_avg=GridF1_avg+metrics/length(sample);
        F1_avg=F1_avg+bestcase/length(sample);
    end
    distance{dd} = [I_best J_best]-[I_score J_score];
    bestval(dd) = F1_avg;
    imagesc(100*mean(GridF1(:,:,selected_index),3))
    hold on
    scatter(J_best(selected_index),I_best(selected_index),300,'sr','filled')
    scatter(J_score(selected_index),I_score(selected_index),150,'ob','filled')
    %     plot([J_best(selected_index) J_score(selected_index)]',[I_best(selected_index) I_score(selected_index)]','-k','LineWidth',4)
    hold off
    if dd==2
        legend('maximum F1 index','eBIC selected index','FontSize',28)
    end
    %     grid on
    text_show = sprintf([type_list_show{1},' ',score_type,' score: %2.1f'],100*F1_avg);
    %     title(text_show)
    axis('square')
    colormap((1-gray).^0.4)
    caxis([0,100])
    set(gca,'xticklabel',[],'yticklabel',[])
    if dd==1
        %                 colorbar('Location','westoutside')
        xlabel('$\lambda_{2}\rightarrow$','Interpreter','latex','FontSize',38)
        
        set(gca,'xaxisLocation','top')
        ylabel('$\leftarrow \lambda_{1}$','Interpreter','latex','FontSize',38)
        
    else
        xlabel('$\lambda_{2}\rightarrow$','Interpreter','latex','FontSize',38)
        
        set(gca,'xaxisLocation','top')
    end
    
    if ii==1
        title([text_label{dd},sprintf('\n avg best F1: %.2f %%',100*F1_avg)],'FontSize',38)
        
    end
end

rr=colorbar('FontSize',18);
set(get(rr,'label'),'string','F1 (%)');
pp = get(0, 'Screensize');
pp(3) = pp(3)*0.8;
set(gcf, 'Position', pp);
% saveas(gcf,[figurepath,'exp_FGN_supp_noncvx'])
% print([figurepath,'exp_FGN_supp_noncvx'],'-painters','-depsc','-r300')
%%
clear
clf
clc
close all
figurepath = './results2plot/figures/';
resource_path = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_S_result\';
load([resource_path,'LLHcorrected_result_adaptive_cvx_formulationS_1percent_lag1_K5_5.mat'])
Mcvx = M;
load([resource_path,'LLHcorrected_result_adaptive_formulationS_1percent_lag1_K5_5.mat'])
Mncvx = M;
width = 100/2000;
figure(2)
tt = tiledlayout(1,2,'padding','compact','tilespacing','none');
nexttile;
tmp = [Mncvx.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];tmp = reshape([tmp.df],30,30)/2000;
histogram(tmp,'BinWidth', width,'Normalization','probability')
title('Differential density: 1%')
% title('non-convex')
% nexttile;
hold on
tmp = [Mcvx.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];tmp = reshape([tmp.df],30,30)/2000;
histogram(tmp,'BinWidth',  width,'Normalization','probability')
% title('convex')
hold off
colormap((gray).^(0.4))

% legend('FGN','cvx-FGN')
xlabel('estimated model density','FontSize',36)
ylabel('Probability','FontSize',36)
ylim([0 0.4])
grid on
set(gca,'FontSize',36)
load([resource_path,'LLHcorrected_result_adaptive_cvx_formulationS_5percent_lag1_K5_5.mat'])
Mcvx = M;
load([resource_path,'LLHcorrected_result_adaptive_formulationS_5percent_lag1_K5_5.mat'])
Mncvx = M;
nexttile;
tmp = [Mncvx.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];tmp = reshape([tmp.df],30,30)/2000;
histogram(tmp,'BinWidth',  width,'Normalization','probability')

title('Differential density: 5%')
hold on
% nexttile;
tmp = [Mcvx.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];tmp = reshape([tmp.df],30,30)/2000;
histogram(tmp, 'BinWidth', width,'Normalization','probability')
hold off
colormap((gray).^(0.4))
legend('FGN','cvx-FGN')
xlabel('estimated model density','FontSize',36)
%     ylabel('probability','FontSize',36)
ylim([0 0.4])
grid on

set(gca,'FontSize',36,'yticklabel',[])
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'exp_FGN_densityhistogram'])
% print([figurepath,'exp_FGN_densityhistogram'],'-depsc','-r300')

%% supplementary material
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
n=5;p=10;K=3;
T= 80;
N = T-p;
Y = randn(n,N,K);
H = randn(n*p,N,K);

[yc,gc] = vectorize_VAR(Y,H,[n,p,K,N]);

IND_DIAG = 1:n+1:n^2; % indices of diagonal elements
P1 = speye(n^2);
P1(IND_DIAG,:) = [];
P = sparse(kron(P1,speye(p*K)));
Dtmp = diffmat(n,p,K);
D = sparse(Dtmp*P);

imagesc(gc ~=0)
pbaspect([size(gc,2) size(gc,1) 1])
colormap(1-gray)
grid on
set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'supplementary_Gmatrix'])
% print([figurepath,'supplementary_Gmatrix'],'-depsc','-r300')

imagesc(P ~=0)
pbaspect([size(P,2) size(P,1) 1])
colormap(1-gray)
grid on
set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'supplementary_Gmatrix'])
% print([figurepath,'supplementary_Pmatrix'],'-depsc','-r300')

imagesc(D ~=0)
pbaspect([size(D,2) size(D,1)  1])
colormap(1-gray)
grid on
set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'supplementary_Gmatrix'])
% print([figurepath,'supplementary_Dmatrix'],'-depsc','-r300')

imagesc(gc'*gc ~=0)
pbaspect([1 1 1])
colormap(1-gray)
grid on
set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'supplementary_GtGmatrix'])
% print([figurepath,'supplementary_GtGmatrix'],'-depsc','-r300')

imagesc(P'*P ~=0)
pbaspect([1 1 1])
colormap(1-gray)
grid on
set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'supplementary_PtPmatrix'])
% print([figurepath,'supplementary_PtPmatrix'],'-depsc','-r300')

imagesc(D'*D ~=0)
pbaspect([1 1 1])
colormap(1-gray)
grid on
set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
set(gcf,'WindowState','fullscreen')
% saveas(gcf,[figurepath,'supplementary_DtDmatrix'])
% print([figurepath,'supplementary_DtDmatrix'],'-depsc','-r300')
%%
n=5;p=10;K=5;
T= 80;
IND_DIAG = 1:n+1:n^2; % indices of diagonal elements
P1 = speye(n^2);
P1(IND_DIAG,:) = [];
P = sparse(kron(P1,speye(p*K)));
Dtmp = diffmat(n,p,K);
D = sparse(Dtmp*P);

%%

figurepath = './results2plot/figures/';
resource_path = './results2plot/';
Connectivity_List = [25,31; ...
    26,32; ...
    27,40];
Adj = zeros(116,116);
Adj(25,31) = -1;
Adj(26,32) = 1;
Adj(27,40) = -1;
% T = table(Adj);
writematrix(Adj,[figurepath,'VisualizeReal.txt'],'Delimiter','\t')

%%
figurepath = './results2plot/figures/';
resource_path = './results2plot/';
Connectivity_List = [25,31; ...
    26,32; ...
    27,40];
Adj = zeros(116,116);
% Adj(25,31) = -1;
% Adj(26,32) = 1;
% Adj(27,40) = -1;

Adj(26,25) = -1;
Adj(6,10) = 1;
Adj(26,28) = 1;
Adj(9,25) = 1;
Adj(26,9) = 1;

% T = table(Adj);
writematrix(Adj,[figurepath,'VisualizeReal_ORB.txt'],'Delimiter','\t')

%%

%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm SECTION %%%%%%%%%%%%%%%%%%%%%%%%%


%% Algorithm performance

clear
clc
inpath = './experiment/model_parameters/';
% outpath = './results2plot/';
% mkdir(outpath)
type = 2; %D type
cd = 3;
T = 100;
p_true = 1;
p_est = 2;
K = 5;
load([inpath,'model_K',int2str(K),'_p',int2str(p_true)]) % struct E
[~,~,dd,m] = size(E);
realz = m;
mname = {'1','5'};
model = E{type,cd,1,1};
y = sim_VAR(model.A,T,1,model.seed,0);
n = size(y,1);
parameter.varorder = p_est;
parameter.formulation = 'fgn'; % cgn, dgn, fgn
parameter.penalty_weight = 'LS'; % LS, uniform
parameter.GridSize = 30;
parameter.data_concat = 0;
parameter.noisecov = 'full'; % 'full', 'diag', 'identity'

eff_T = T-p_est;
H = zeros(n*p_est,eff_T,K);
Y = zeros(n,eff_T,K);
disp('Generating H, Y matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p_est);
end
disp('vectorizing model')
[b,G] = vectorize_VAR(Y,H,[n,p_est,K,eff_T]);
xLS = G\b;
[P,~] = offdiagJSS(n,p_est,K);
Dtmp = diffmat(n,p_est,K);
D = Dtmp*P;
L1 = sparse(P);

if strcmp(parameter.formulation,'fgn')
    MODEL_DIM = [n, p_est,K,p_est,p_est];
    L2 = sparse(D);
else
    MODEL_DIM = [n, p_est,K,p_est,p_est*K];
    L2 = sparse(P);
end






[Lambdacrit_1.cvx,Lambdacrit_2.cvx] = gen_critical_lambdas(G,b, xLS,MODEL_DIM,'cvx',parameter.penalty_weight,parameter.formulation);
[Lambdacrit_1.ncvx,Lambdacrit_2.ncvx] = gen_critical_lambdas(G,b, xLS,MODEL_DIM,'ncvx',parameter.penalty_weight,parameter.formulation);
Lambda = logspace(-6,0,parameter.GridSize);
ii = 21;
jj = 23;
a1.cvx = Lambdacrit_1.cvx*Lambda(ii);
a2.cvx = Lambdacrit_2.cvx*Lambda(jj);
a1.ncvx = Lambdacrit_1.ncvx*Lambda(ii);
a2.ncvx = Lambdacrit_2.ncvx*Lambda(jj);

if strcmp(parameter.formulation,'cgn')
    a1.cvx = a1.cvx*0;
    a1.ncvx = a1.ncvx*0;
end


ALG_PARAMETER.cvx = gen_alg_params('cvx', parameter.formulation);
ALG_PARAMETER.cvx.PRINT_RESULT = 0;
ALG_PARAMETER.cvx.dim = MODEL_DIM;
ALG_PARAMETER.cvx.x0 = xLS;
[x.cvx,~,~,history.cvx] = spectral_ADMM(G, b,a1.cvx,a2.cvx,L1,L2,ALG_PARAMETER.cvx);
ALG_PARAMETER.ncvx = gen_alg_params('ncvx', parameter.formulation);
ALG_PARAMETER.ncvx.PRINT_RESULT = 0;
ALG_PARAMETER.ncvx.dim = MODEL_DIM;
ALG_PARAMETER.ncvx.x0 = xLS;

[x.ncvx,~,~,history.ncvx] = adaptive_ADMM(G, b,a1.ncvx,a2.ncvx,L1,L2,ALG_PARAMETER.ncvx);
%%
clf
figurepath = './results2plot/figures/';
formulation = parameter.formulation;

% tt=tiledlayout(2,1);
figure(1);
% tmp = diff(history.cvx.objval);

semilogy(history.cvx.reldiff_norm,'linewidth',2)
grid on
xlabel('iterations $(k)$','Interpreter','latex')
h=ylabel('$\frac{\Vert x_{k+1}-x_{k} \Vert_{2}}{\Vert x_{k} \Vert_{2}}$','Interpreter','latex');
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',36)
set(h, 'FontSize', 60)
set(gcf,'WindowState','fullscreen')
ylim([10^-7,2*sqrt(2)])
print([figurepath,'algperf_xreldiff_cvx_',formulation],'-depsc','-r300')


figure(2);
semilogy(abs(((history.cvx.objval)-history.cvx.objval(end))/history.cvx.objval(end)),'linewidth',2)
yyaxis right
plot(history.cvx.rho,'linewidth',2)
ylabel('$\rho$','Interpreter','latex')
yyaxis left
grid on
xlabel('iterations $(k)$','Interpreter','latex')
h=ylabel('$\frac{F(x_{k})-p^{*}}{p^{*}}$','Interpreter','latex');
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',36)
set(h, 'FontSize', 60)
set(gcf,'WindowState','fullscreen')
ylim([10^-11,2*sqrt(2)])
print([figurepath,'algperf_objreldiff_cvx_',formulation],'-depsc','-r300')


figure(3);
% tmp = diff(history.ncvx.objval);
semilogy(history.ncvx.reldiff_norm,'linewidth',2)
grid on
xlabel('iterations $(k)$','Interpreter','latex')
h=ylabel('$\frac{\Vert x_{k+1}-x_{k} \Vert_{2}}{\Vert x_{k} \Vert_{2}}$','Interpreter','latex');
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',36)
set(h, 'FontSize', 60)
set(gcf,'WindowState','fullscreen')
ylim([10^-7,2*sqrt(2)])
print([figurepath,'algperf_xreldiff_ncvx_',formulation],'-depsc','-r300')


figure(4);
semilogy(abs(((history.ncvx.objval)-history.ncvx.objval(end))/history.ncvx.objval(end)),'linewidth',2)
yyaxis right
plot(history.ncvx.rho,'linewidth',2)
ylabel('$\rho$','Interpreter','latex')
yyaxis left
grid on
xlabel('iterations $(k)$','Interpreter','latex')
h=ylabel('$\frac{F(x_{k})-p^{*}}{p^{*}}$','Interpreter','latex');
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',36)
set(h, 'FontSize', 60)
set(gcf,'WindowState','fullscreen')
ylim([10^-8,2*sqrt(2)])
print([figurepath,'algperf_objreldiff_ncvx_',formulation],'-depsc','-r300')

%% convex CGN comparison
PARAMETER.proj_ind = P*(1:1:n^2*p_est*K)';
PARAMETER.p=2;
PARAMETER.q=1;
PARAMETER.dim(1)=n;
PARAMETER.dim(2)=p_est;
PARAMETER.dim(3)=K;
PARAMETER.delta=0.7;
PARAMETER.eta=0.5;
PARAMETER.rho=0.5;
PARAMETER.IS_LINESEARCH = 1;
[x.cvx_nmAPG,history.cvx_nmAPG]= nmAPG_BB_cgn(G,b,a2.cvx,PARAMETER,xLS);
PARAMETER.q=0.5;
[x.ncvx_nmAPG,history.ncvx_nmAPG]= nmAPG_BB_cgn(G,b,a2.ncvx,PARAMETER,xLS);

PARAMETER.IS_LINESEARCH = 0;
PARAMETER.q=1;
[x.cvx_nmAPG_noBB,history.cvx_nmAPG_noBB]= nmAPG_BB_cgn(G,b,a2.cvx,PARAMETER,xLS);
PARAMETER.q=0.5;
[x.ncvx_nmAPG_noBB,history.ncvx_nmAPG_noBB]= nmAPG_BB_cgn(G,b,a2.ncvx,PARAMETER,xLS);

ALG_PARAMETER.cvx = gen_alg_params('cvx', 'cgn');
% ALG_PARAMETER.cvx.epscor = 0.01;
ALG_PARAMETER.cvx.PRINT_RESULT = 0;
ALG_PARAMETER.cvx.dim = MODEL_DIM;
ALG_PARAMETER.cvx.x0 = xLS;
ALG_PARAMETER.ncvx = gen_alg_params('ncvx', 'cgn');
ALG_PARAMETER.ncvx.PRINT_RESULT = 0;
ALG_PARAMETER.ncvx.dim = MODEL_DIM;
ALG_PARAMETER.ncvx.x0 = xLS;

[x.cvx,~,~,history.cvx] = spectral_ADMM(G, b,a1.cvx*0,a2.cvx,L1,L2,ALG_PARAMETER.cvx);

[x.ncvx,~,~,history.ncvx] = adaptive_ADMM(G, b,a1.ncvx*0,a2.ncvx,L1,L2,ALG_PARAMETER.ncvx);

ALG_PARAMETER.static = ALG_PARAMETER.cvx;
ALG_PARAMETER.static.IS_ADAPTIVE = 0;
ALG_PARAMETER.static.rho_init = 600;
[x.cvx_static,~,~,history.cvx_static] = spectral_ADMM(G, b,a1.cvx*0,a2.cvx,L1,L2,ALG_PARAMETER.static);
%%
figure(1)
tt=tiledlayout(2,1);
nexttile;
semilogy(history.cvx.reldiff_norm,'linewidth',3)
hold on
semilogy(history.cvx_nmAPG.reldiff_norm,'linewidth',3)
semilogy(history.cvx_nmAPG_noBB.reldiff_norm,'linewidth',3)
hold off
xlabel('iterations')
ylabel('$\frac{\Vert x^{+}-x \Vert_{2}}{\Vert x \Vert_{2}}$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Spectral ADMM', 'nmAPG with BB', 'nmAPG no BB','location','eastoutside')
nexttile;
semilogy(((history.cvx.objval)-history.cvx.objval(end))/history.cvx.objval(end),'linewidth',3)
grid on
hold on
semilogy(((history.cvx_nmAPG.objval)-history.cvx_nmAPG.objval(end))/history.cvx_nmAPG.objval(end),'linewidth',3)
semilogy(((history.cvx_nmAPG_noBB.objval)-history.cvx_nmAPG_noBB.objval(end))/history.cvx_nmAPG_noBB.objval(end),'linewidth',3)
grid on
hold off
xlabel('iterations')
ylabel('$\frac{f(x^{+})-p^{*}}{p^{*}}$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Spectral ADMM', 'nmAPG with BB', 'nmAPG no BB','location','eastoutside')

figure(2)
tt=tiledlayout(2,1);
nexttile;
semilogy(history.ncvx.reldiff_norm,'linewidth',3)
grid on
hold on
semilogy(history.ncvx_nmAPG.reldiff_norm,'linewidth',3)
semilogy(history.ncvx_nmAPG_noBB.reldiff_norm,'linewidth',3)
grid on
hold off
xlabel('iterations')
ylabel('$\frac{\Vert x^{+}-x \Vert_{2}}{\Vert x \Vert_{2}}$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Adaptive ADMM', 'nmAPG with BB', 'nmAPG no BB','location','eastoutside')
nexttile;
semilogy(((history.ncvx.objval)-history.ncvx.objval(end))/history.ncvx.objval(end),'linewidth',3)
grid on
hold on
semilogy(((history.ncvx_nmAPG.objval)-history.ncvx_nmAPG.objval(end))/history.ncvx_nmAPG.objval(end),'linewidth',3)
semilogy(((history.ncvx_nmAPG_noBB.objval)-history.ncvx_nmAPG_noBB.objval(end))/history.ncvx_nmAPG_noBB.objval(end),'linewidth',3)
grid on
hold off
xlabel('iterations')
ylabel('$\frac{f(x^{+})-p^{*}}{p^{*}}$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Adaptive ADMM', 'nmAPG with BB', 'nmAPG no BB','location','eastoutside')
%% time usage nonconvex

figure(3)
% plot([0,history.ncvx.t],history.ncvx.objval,'linewidth',3)
% hold on
% plot([0;history.ncvx_nmAPG.t],history.ncvx_nmAPG.objval,'linewidth',3)
% plot([0;history.ncvx_nmAPG_noBB.t],history.ncvx_nmAPG_noBB.objval,'linewidth',3)
% hold off
plot(history.ncvx.objval,'linewidth',3)
hold on
plot(history.ncvx_nmAPG.objval,'linewidth',3)
plot(history.ncvx_nmAPG_noBB.objval,'linewidth',3)
hold off
xlabel('iteration')
ylabel('objective')
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
grid on
legend('Adaptive ADMM', 'nmAPG with BB', 'nmAPG no BB')
print([figurepath,'algperf_CGN'],'-dsvg','-r300')
%% time usage convex

figure(4)
% h1=plot([0,history.cvx.t],history.cvx.objval,'linewidth',3);
% hold on
% h2=plot([0,history.cvx_static.t],history.cvx_static.objval,'linewidth',3);
% plot([0;history.cvx_nmAPG.t],history.cvx_nmAPG.objval,'linewidth',3)
% plot([0;history.cvx_nmAPG_noBB.t],history.cvx_nmAPG_noBB.objval,'linewidth',3)
% hold off

h1=plot(history.cvx.objval,'linewidth',3);
hold on
h2=plot(history.cvx_static.objval,'linewidth',3);
h3=plot(history.cvx_nmAPG.objval,'linewidth',3);
h4=plot(history.cvx_nmAPG_noBB.objval,'linewidth',3);
hold off

xlabel('iteration')
ylabel('objective')
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',24)
legend('Spectral ADMM', 'fixed ADMM', 'nmAPG with BB', 'nmAPG no BB')
uistack(h4,'top')
uistack(h3,'top')
uistack(h2,'top')
uistack(h1,'top')

print([figurepath,'algperf_cvxCGN'],'-dsvg','-r300')
%% vectorization visualization

n = 3;p=5;K=3;
T = 10;
y = randn(n,T,K);
x = ones(n^2,1);
x_tmp = zeros(p,K);
for kk=1:K
    x_tmp(:,kk) = kk;
end
x = kron(x,x_tmp(:));
b=(1:K)';
b = kron(b,ones(n*T,1));

N = T-p;
H = zeros(n*p,N,K);
Y = zeros(n,N,K);
disp('Generating H, Y matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
end

G = [];
%blkH = zeros(Num,n*p,K);
BLKH = zeros(N,n*p*K,K);
E1 = zeros(K,1); E1(1) = 1;
PLOT_G = zeros(n*N, p*n^2*K,K);
for k=1:K
    H0 = zeros(n,N,p);
    for j=1:p
        H0(:,:,j) = H((j-1)*n+1:j*n,:,k);
    end
    ind = linindex(n,N,p,'col');
    vecH = H0(ind);
    TMP = zeros(n*p*K*N,1);
    [~,IND] = blocksub(TMP,p,p*(K-1));
    TMP(IND) = vecH;
    TMP = reshape(TMP,n*p*K,N);
    TMP = TMP';
    BLKH(:,:,k) = sparse([zeros(N,(k-1)*p) TMP(:,1:end-(k-1)*p)]);
    BLKG = kron(speye(n),BLKH(:,:,k)); % not efficient when n is large
    % size of BLKG is nN x pn^2K
    G = [G;BLKG];
    PLOT_G(:,:,k) = BLKG;
end
G = sparse(G);

tt = tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
color = {[1,0,0],[0,0.5,0],[0,0,1]};
for k=1:K
    nexttile(3*k-2);
    hAxes = gca;
    imagesc(hAxes, ones(n*T,1))
    colormap( hAxes , [1 1 1;color{k}] )
    pbaspect([6 length(b) 1])
    set(gca,'xticklabel',[],'yticklabel',[])
    ylabel(sprintf('$v^{(%d)}$',k),'Interpreter','latex')
    %     axis('square')
    grid on
    grid minor
    set(gca,'FontSize',30)
    if k==K
        xlabel('$b$','Interpreter','latex')
    end
    
    nexttile(3*k-1);
    hAxes = gca;
    imagesc(hAxes, PLOT_G(:,:,k)~=0)
    colormap( hAxes , [1 1 1;color{k}] )
    set(gca,'xticklabel',[],'yticklabel',[])
    ylabel(sprintf('$G^{(%d)}$',k),'Interpreter','latex')
    %     axis('square')
    grid on
    grid minor
    set(gca,'FontSize',30)
end
xlabel('$G$','Interpreter','latex')
nexttile(3,[3,1]);
hAxes = gca;
imagesc(hAxes, x)
colormap( hAxes , [1,0,0;0,0.5,0;0,0,1] )
pbaspect([5 length(x) 1])
set(gca,'xticklabel',[],'yticklabel',[])
xlabel('$x$','Interpreter','latex')
set(gca,'FontSize',30)
grid on
grid minor

%% REVIEWER RESPONSE RV. 1


%% plot ground-truth
clear
clc
model_parameter_path = './experiment/model_parameters/';
n = 20; % time-series dimension
p = 2;  % ground-truth VAR order
K = 3; % number of models
realization = 1; % number of model's realizations
common_density = [0.1]; % for p=1, common GC network density [We actually used only 0.1, 0.2]
differential_density = [0.05]; % differential GC network density
model = {'common','differential','similar'}; % type of similarity
mname = {'C','D','S'};
cnt = 0;
E=cell(length(model),length(common_density),length(differential_density),realization);
for m=1:length(model)
    for d=1:length(common_density)
        opts.common_density = common_density(d);
        for diff_d =1:length(differential_density)
            opts.differential_density = differential_density(diff_d);
            opts.type = model{m};
            
            for b=1:realization %number of [C,S,D] VAR model generated
                if strcmp(mname{m},'D')
                    E{m,d,diff_d,b} = gen_multi_VAR([n,p,K],opts,E{1,d,diff_d,b}.A); % look for C type model [code 1] to generate D type
                else
                    E{m,d,diff_d,b} = gen_multi_VAR([n,p,K],opts);
                end
            end
        end
    end
end

sample_model = 1;
tt = tiledlayout(K,3);
tt.TileSpacing = 'None';
tt.Padding = 'None';
mname = {'common','differential','fused'};
figurepath = './results2plot/figures/';
for ii=1:3
    for kk =1:K
        
        model = E{ii,1,1,sample_model};
        GC = model.GC(:,:,1:K);
        A = model.A;
        A = squeeze(sqrt(sum(A.^2,3)));
        
        
        
        [n,~,K] = size(A);
        commonNZ = ones(n,n)-eye(n);
        diag_ind = 1:n+1:n^2;
        for zz=1:K
            tmp = A(:,:,zz);
            commonNZ = commonNZ & (tmp~=0);
        end
        s = [];
        
        tmp = A(:,:,kk);
        nexttile;
        
        imagesc(tmp)
        set(gca,'box','on', 'BoxStyle','full')
        
        
        hold on
        spy(tmp.*(1-commonNZ),30,'r')
        hold off
        
        set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
        if ii==1
            title(sprintf('model #%d', kk))
            
        end
        if kk==1
            ylabel(mname{ii})
        end
    end
    hold off
    colormap((1-gray))
    
    
end
pp = get(0, 'Screensize');
pp(3) = pp(3)*1/3;
pp(4) = pp(4)*1/3;
pp = pp+200;
set(gcf, 'Position', pp);
set(findall(gcf,'-property','FontSize'),'FontSize',28)
% print([figurepath,'reviewer_response_sampleGC'],'-painters','-depsc','-r300')
% print([figurepath,'reviewer_response_sampleGC'],'-painters','-dpng','-r300')

%% fMRI connectogram visualization [D2K, F2K, C18K]
clear
clc
clf
close all
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('.\experiment\ADHD\AAL_116.mat')


filename_list = {'estim_C18K', 'estim_D2K', 'estim_F2K'};
for file_index = 1:length(filename_list)
    clf
    close all
    file_name  = filename_list{file_index};
    load([inpath,file_name])
    Right_Left_permutation = [116,115,114,113,112,111,110,109,107,105,103,101,99,97,95,93,91,89,87,85,83,81,79,77,75,73,71,69,67,65,63,61,59,57,55,53,51,49,47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108];
    
    if strcmp(file_name, 'estim_C18K')
        unstable_TDC = [5, 6, 11, 14, 17];
        unstable_ADHD = [5];
        GC_TDC = M.TDC.model(M.TDC.index.eBIC).GC;
        GC_TDC(:,:,unstable_TDC) = [];
        GC_ADHD = M.TDC.model(M.TDC.index.eBIC).GC;
        GC_ADHD(:,:,unstable_ADHD) = [];
        
        
        GC.total(:,:,1) = mean(GC_TDC,3);
        GC.total(:,:,2) = mean(GC_ADHD,3);
    else
        GC.total = M.model(M.index.eBIC).GC;
    end
    
    GC_TDC = GC.total(:,:,1).*(1-eye(116));
    GC_ADHD = GC.total(:,:,2).*(1-eye(116));
    selected_index_TDC = find(GC_TDC);
    selected_index_ADHD = find(GC_ADHD);
    [~,I_TDC] = sort(abs(GC_TDC(selected_index_TDC)),'descend');
    [~,I_ADHD] = sort(abs(GC_ADHD(selected_index_ADHD)),'descend');
    
    GC_TDC_toplot = zeros(116,116);
    GC_ADHD_toplot = zeros(116,116);
    
    first_n_rank = 50;
    
    for ii=1:first_n_rank
        % for ii=1:length(selected_index_TDC)
        [ee,cc] = ind2sub([116,116],selected_index_TDC(I_TDC(ii)));
        
        
        GC_TDC_toplot(ee, cc) = GC.total(ee, cc, 1);
        
        
    end
    
    for ii=1:first_n_rank
        % for ii=1:length(selected_index_ADHD)
        
        [ee,cc] = ind2sub([116,116],selected_index_ADHD(I_ADHD(ii)));
        GC_ADHD_toplot(ee, cc) = GC.total(ee, cc, 2);
        
        
    end
    ROI_color = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1];
    
    Right_Left_permutation = Right_Left_permutation([89:116, 1:88]);
    ROI_color = ROI_color([89:116, 1:88]);
    myColorMap = zeros(116,3);
    
    myLabel = AAL_116.name(Right_Left_permutation);
    
    ROI_name = {'Posterior Fossa', ...
        'Temporal lobe', ...
        'Central Structures', ...
        'Parietal Lobe', ...
        'Occipital Lobe', ...
        'Insula and Cingulate Gyri', ...
        'Frontal Lobe'};
    
    for jj = 1:116
        newStr = strrep(myLabel{jj},'_',' ');
        myLabel{jj} = newStr;
        scale = 0.8;
        switch ROI_color(jj)
            case 1
                xx = [1 0 0].*scale;
            case 2
                xx = [0 1 0].*scale;
            case 3
                xx = [0 0 1].*scale;
            case 4
                xx = [1 1 0].*scale;
            case 5
                xx = [1 0 1].*scale;
            case 6
                xx = [0 1 1].*scale;
            case 7
                xx = [0.5 0.6 0.8].*scale;
        end
        myColorMap(jj, :) = xx;
        
    end
    figurepath = './results2plot/figures/';
    figure(1)
    % nexttile;
    % myColorMap = lines(length(GC_TDC_toplot));
    h1=circularGraph(GC_TDC_toplot(Right_Left_permutation, Right_Left_permutation)','Colormap',myColorMap,'Label',myLabel);
    pp = get(0, 'Screensize');
    pp(3) = pp(3)*1;
    set(gcf, 'Position', pp);
    delete(h1.ShowButton)
    delete(h1.HideButton)
    exportgraphics(gcf,[figurepath,'reviewer_response_',file_name,'_circular_TDC.png'],'ContentType','vector')
    exportgraphics(gcf,[figurepath,'reviewer_response_',file_name,'_circular_TDC.eps'],'ContentType','vector')
    
    % nexttile;
    figure(2)
    h2=circularGraph(GC_ADHD_toplot(Right_Left_permutation, Right_Left_permutation)','Colormap',myColorMap,'Label',myLabel);
    pp = get(0, 'Screensize');
    pp(3) = pp(3)*1;
    delete(h2.ShowButton)
    delete(h2.HideButton)
    set(gcf, 'Position', pp);
    exportgraphics(gcf,[figurepath,'reviewer_response_',file_name,'_circular_ADHD.png'])
    exportgraphics(gcf,[figurepath,'reviewer_response_',file_name,'_circular_ADHD.eps'])
end

figure(3)
hold on
for jj=1:7
    switch jj
        case 1
            xx = [1 0 0].*scale;
        case 2
            xx = [0 1 0].*scale;
        case 3
            xx = [0 0 1].*scale;
        case 4
            xx = [1 1 0].*scale;
        case 5
            xx = [1 0 1].*scale;
        case 6
            xx = [0 1 1].*scale;
        case 7
            xx = [0.5 0.6 0.8].*scale;
    end
    plot(NaN,NaN,'.','Color',xx, 'MarkerSize',69);
end
hold off
lh = legend(ROI_name);
lh.Position(1) = 0.5 - lh.Position(3)/2;
lh.Position(2) = 0.518 - lh.Position(4)/2;
set(gca,'FontSize', 32, 'xticklabel',[], 'yticklabel',[])

exportgraphics(gcf,[figurepath,'reviewer_response_fMRI_legend.png'],'ContentType','vector')
exportgraphics(gcf,[figurepath,'reviewer_response_fMRI_legend.eps'],'ContentType','vector')
%% fMRI connectogram visualization [D2K bootstrap]

clear
clc
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('.\experiment\ADHD\AAL_116.mat')

bootstrap_path = dir('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_D2K_seed_*_bootstrap_*_DIAGCOV.mat');
GC.total = zeros(116,116, 2);
GC.repetition = zeros(116, 116, 2);

bootstrap_sample = length(bootstrap_path);

for xx=1:bootstrap_sample
    fprintf("bootstrap sample No. %d\n", xx)
    
    load([inpath,bootstrap_path(xx).name])
    
    GC.total = GC.total +  (M.model(M.index.eBIC).GC)/bootstrap_sample;  % to filter with mean
    GC.repetition = GC.repetition + (M.model(M.index.eBIC).GC~=0)/bootstrap_sample;  % to filter with occurence
end


clf
close all
GC_TDC = GC.total(:,:,1).*(1-eye(116)).*(GC.repetition(:, :, 1)>0.95);
GC_ADHD = GC.total(:,:,2).*(1-eye(116)).*(GC.repetition(:, :, 2)>0.95);
GC_TDC_toplot = GC_TDC;
GC_ADHD_toplot = GC_ADHD;

ROI_color = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1];
Right_Left_permutation = [116,115,114,113,112,111,110,109,107,105,103,101,99,97,95,93,91,89,87,85,83,81,79,77,75,73,71,69,67,65,63,61,59,57,55,53,51,49,47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108];
Right_Left_permutation = Right_Left_permutation([89:116, 1:88]);
ROI_color = ROI_color([89:116, 1:88]);
myColorMap = zeros(116,3);

myLabel = AAL_116.name(Right_Left_permutation);

ROI_name = {'Posterior Fossa', ...
    'Temporal lobe', ...
    'Central Structures', ...
    'Parietal Lobe', ...
    'Occipital Lobe', ...
    'Insula and Cingulate Gyri', ...
    'Frontal Lobe'};

for jj = 1:116
    newStr = strrep(myLabel{jj},'_',' ');
    myLabel{jj} = newStr;
    scale = 0.8;
    switch ROI_color(jj)
        case 1
            xx = [1 0 0].*scale;
        case 2
            xx = [0 1 0].*scale;
        case 3
            xx = [0 0 1].*scale;
        case 4
            xx = [1 1 0].*scale;
        case 5
            xx = [1 0 1].*scale;
        case 6
            xx = [0 1 1].*scale;
        case 7
            xx = [0.5 0.6 0.8].*scale;
    end
    myColorMap(jj, :) = xx;
    
end
figurepath = './results2plot/figures/';
figure(1)
% nexttile;
% myColorMap = lines(length(GC_TDC_toplot));
h1=circularGraph(GC_TDC_toplot(Right_Left_permutation, Right_Left_permutation)','Colormap',myColorMap,'Label',myLabel);
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
delete(h1.ShowButton)
delete(h1.HideButton)
exportgraphics(gcf,[figurepath,'reviewer_response_estim_D2Kbootstrap_circular_TDC_first50.png'],'ContentType','vector')
exportgraphics(gcf,[figurepath,'reviewer_response_estim_D2Kbootstrap_circular_TDC_first50.eps'],'ContentType','vector')

figure(2)
h2=circularGraph(GC_ADHD_toplot(Right_Left_permutation, Right_Left_permutation)','Colormap',myColorMap,'Label',myLabel);
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
delete(h2.ShowButton)
delete(h2.HideButton)
exportgraphics(gcf,[figurepath,'reviewer_response_estim_D2Kbootstrap_circular_ADHD_first50.png'],'ContentType','vector')
exportgraphics(gcf,[figurepath,'reviewer_response_estim_D2Kbootstrap_circular_ADHD_first50.eps'],'ContentType','vector')
%% case n=100

clear
clc
performance_path = './results2plot/';

load([performance_path,'reviewer_response_CVX_CGN_result'],'R')
load([performance_path,'reviewer_response_CVX_CGN_ALL_RESULT'],'ALL_RESULT')
R_cvx = R;
ALL_RESULT_CVX = ALL_RESULT;

load([performance_path,'reviewer_response_CGN_result'],'R')
load([performance_path,'reviewer_response_CGN_ALL_RESULT'],'ALL_RESULT')
R_ncvx=R;

ALL_RESULT_NCVX = ALL_RESULT;

group_index = [ones(1, 100) 2*ones(1, 100)];

group_index = cell(200, 1);
group_index(1:100) = {'CGN'};
group_index(101:200) = {'cvx-CGN'};

performance_name = {'F1', 'MCC', 'TPR', 'FPR', 'ACC'};

ARRcvx_to_plot = zeros(100, 6);
ARRncvx_to_plot = zeros(100, 6);

tt = tiledlayout(1,6);

for ii = 1:length(performance_name)
    tmp = R_cvx.common.(performance_name{ii})';
    DATA = [R_ncvx.common.(performance_name{ii})'; tmp];
    
    
    ax(ii) = nexttile;
    
    
    
    boxplot(100*DATA, group_index)
    title(performance_name{ii})
    
    if ii == 1
        ylabel('performance index')
    end
    if ii == 4
        ylim([0 1.2])
    end
    if ii == 3
        ylim([90 100])
    end
    if ii == 5
        ylim([95 100])
    end
    set(findobj(gca,'type','line'),'linew',3)
    ytickformat('percentage')
    xtickangle(15)
    %     a = get(gca,'XTickLabel');
    %     set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
end
linkaxes([ax(1) ax(2)])
% linkaxes([ax(3) ax(5)])
% for ii = 1:length(performance_name)
%     ax(ii).YTickLabel = strcat(ax(ii).YTickLabel,'%');
%     ax(ii).XAxis.FontSize = 10;
% end

DATA = [R_ncvx.bias'; R_cvx.bias'];


axb = nexttile;
boxplot(100*DATA, group_index);
set(findobj(gca,'type','line'),'linew',3)
axb.YTickLabel = strcat(axb.YTickLabel,'%');
title('Estimation error')
xtickangle(15)

pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
pp(4) = pp(4)*75;
pp(2) = pp(2) + 100;
set(gcf, 'Position', pp);
set(findall(gcf,'-property','FontSize'),'FontSize',28)
for ii = 1:length(performance_name)
    %     ax(ii).YTickLabel = strcat(ax(ii).YTickLabel,'%');
    ax(ii).XAxis.FontSize = 25;
end
axb.XAxis.FontSize = 25;
figurepath = './results2plot/figures/';
exportgraphics(gcf,[figurepath,'n100_eBIC_results.eps'],'ContentType','vector')
exportgraphics(gcf,[figurepath,'n100_eBIC_results.png'])


%% plot ROC
clc
close all
clear
performance_path = './results2plot/';

load([performance_path,'reviewer_response_CVX_CGN_result'],'R')
load([performance_path,'reviewer_response_CVX_CGN_ALL_RESULT'],'ALL_RESULT')
R_cvx = R;
ALL_RESULT_CVX = ALL_RESULT;

load([performance_path,'reviewer_response_CGN_result'],'R')
load([performance_path,'reviewer_response_CGN_ALL_RESULT'],'ALL_RESULT')
R_ncvx=R;

ALL_RESULT_NCVX = ALL_RESULT;

figure(1)
% cvx-CGN

FPR_array = zeros(1,30);
TPR_array = zeros(1,30);
for ii=1:length(ALL_RESULT_CVX)
    tmp = [ALL_RESULT_CVX(ii).model_acc.common]; FPR = [tmp.FPR]; TPR=[tmp.TPR];
    
    FPR_array = FPR_array+FPR/100;
    TPR_array = TPR_array+TPR/100;
end
plot(FPR_array, TPR_array, '.-b', 'LineWidth', 4)
hold on

FPR_array = zeros(1,30);
TPR_array = zeros(1,30);
for ii=1:length(ALL_RESULT_NCVX)
    if isempty(ALL_RESULT_NCVX(ii).model_acc)
        continue
    end
    tmp = [ALL_RESULT_NCVX(ii).model_acc.common]; FPR = [tmp.FPR]; TPR=[tmp.TPR];
    
    FPR_array = FPR_array+FPR/100;
    TPR_array = TPR_array+TPR/100;
end
plot(FPR_array, TPR_array, '--r', 'LineWidth', 4)
hold off

legend('cvx-CGN', 'CGN')
xlabel('FPR')
ylabel('TPR')
xlim([0 1])
ylim([0 1])
axis('square')


pp = get(0, 'Screensize');
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
set(findall(gcf,'-property','FontSize'),'FontSize',28)
figurepath = './results2plot/figures/';
exportgraphics(gcf,[figurepath,'n100_ROC_results.eps'],'ContentType','vector')
exportgraphics(gcf,[figurepath,'n100_ROC_results.png'])