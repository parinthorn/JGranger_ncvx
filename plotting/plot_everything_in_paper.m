%% 3.2: Common network comparison (CommonGrangerNet, convex CommonGrangerNet, JSS, Magda Gregorova)
clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'CGN','cvxCGN','Song','Greg'};
table_head_show = {'CGN','cvx-CGN','Songsiri, 17','Gregorova, 15'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'adaptive_formulation_C_result'])
result.CGN = R;
load([resource_path,'adaptive_formulation_C_cvx_result_K5'])
result.cvxCGN = R;
load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\adaptive_formulation_C_JSS_result.mat'])
result.Song = R;
load([resource_path,'magda_result'])
result.Greg = R;
M = zeros(5,4,2);
STD = zeros(5,4,2);
summary = zeros(5,4,2,100);
for jj=1:length(table_head)
    for ii=1:length(row_name)
        for dd=1:2
            if  jj~=4
                M(ii,jj,dd) = mean(result.(table_head{jj}).common.(row_name{ii})(dd,:));
                STD(ii,jj,dd) = std(result.(table_head{jj}).common.(row_name{ii})(dd,:));
                summary(ii,jj,dd,:) = result.(table_head{jj}).common.(row_name{ii})(dd,:);
            else
                M(ii,jj,dd) = mean(result.(table_head{jj}).(row_name{ii})(dd,:));
                STD(ii,jj,dd) = std(result.(table_head{jj}).(row_name{ii})(dd,:));
                summary(ii,jj,dd,:) = result.(table_head{jj}).(row_name{ii})(dd,:);
            end
        end
    end
end

% table
printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name)

hh = tiledlayout(2,2);
hh.TileSpacing = 'compact';
hh.Padding = 'compact';
density = {'10%','20%'};

for ii=1:2
    for dd=1:2
        ARR = zeros(100,4);
        for jj=1:length(table_head)
            ARR(:,jj) = summary(ii,jj,dd,:);
        end
        nexttile;
        boxplot(100*ARR);
        
        set(findobj(gca,'type','line'),'linew',3)
        set(gca,'xticklabel',table_head_show)
        grid on
        if ii==1
            title(['common density: ', density{dd}])
            
        end
        if dd==1
            ylabel([row_name{ii},'(%)'])
        end
    end
    
    
end
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf,[figurepath,'3_2'])
%     print([figurepath,'3_2'],'-depsc')
%% 3.3 A

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'CGN','cvxCGN','Song','Skrip'};
table_head_show = {'CGN','cvx-CGN','Songsiri, 17','Skripikov, 19b'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'adaptive_formulation_D_result_K5'])
result.CGN = R;
load([resource_path,'adaptive_formulation_D_cvx_result_K5'])
result.cvxCGN = R;
load([resource_path,'adaptive_formulation_D_JSS_result_K5'])
result.Song = R;
load([resource_path,'ResultSkripD_K5'])
result.Skrip = R;
M = zeros(5,4,2);
STD = zeros(5,4,2);
summary = zeros(5,4,2,100);
for jj=1:length(table_head)
    for ii=1:length(row_name)
        for dd=1:2
            M(ii,jj,dd) = mean(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            STD(ii,jj,dd) = std(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            summary(ii,jj,dd,:) = result.(table_head{jj}).total.(row_name{ii})(dd,:);
            
        end
    end
end

% table
printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name)

tt = tiledlayout(2,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
density = {'1%','5%'};
for ii=1:2
    for dd=1:2
        ARR = zeros(100,4);
        for jj=1:length(table_head)
            ARR(:,jj) = summary(ii,jj,dd,:);
        end
        nexttile;
        boxplot(100*ARR);
        set(findobj(gca,'type','line'),'linew',3)
        set(gca,'xticklabel',table_head_show)
        grid on
        if ii==1
            title(['differential density: ', density{dd}])
        end
        if dd==1
            ylabel([row_name{ii},'(%)'])
        end
    end
end
set(gcf, 'Position', get(0, 'Screensize'));

set(findall(gcf,'-property','FontSize'),'FontSize',24)

%     saveas(gcf,[figurepath,'3_3_A'])
%     print([figurepath,'3_3_A'],'-depsc')
%% 3.3 B
clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'CGN','cvxCGN','Song','Skrip'};
table_head_show = {'CGN','cvx-CGN','Songsiri, 17','Skripikov, 19b'};
row_name = {'F1','FPR','TPR','ACC','MCC'};

load([resource_path,'adaptive_formulation_D_result_K5'])
result.CGN = R;
load([resource_path,'adaptive_formulation_D_cvx_result_K5'])
result.cvxCGN = R;
load([resource_path,'adaptive_formulation_D_JSS_result_K5'])
result.Song = R;
load([resource_path,'ResultSkripD_K5'])
result.Skrip = R;


load([resource_path,'adaptive_formulation_D_result_K50'])
result.CGN = R;
load([resource_path,'adaptive_formulation_D_cvx_result_K50'])
result.cvxCGN = R;
load([resource_path,'adaptive_formulation_D_JSS_result_K50'])
result.Song = R;
load([resource_path,'ResultSkripD_K50'])
result.Skrip = R;
M = zeros(5,4);
STD = zeros(5,4);
summary = zeros(5,4,100);
for jj=1:length(table_head)
    for ii=1:length(row_name)
        dd=2;
            M(ii,jj) = mean(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            STD(ii,jj) = std(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            summary(ii,jj,:) = result.(table_head{jj}).total.(row_name{ii})(dd,:);
            
        
    end
end

% table
printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name)

tt = tiledlayout(2,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
density = {'1%','5%'};
for ii=1:2
    for dd=1:2
        ARR = zeros(100,4);
        for jj=1:length(table_head_show)
            ARR(:,jj) = summary(ii,jj,dd,:);
        end
        nexttile;
        boxplot(100*ARR);
        set(findobj(gca,'type','line'),'linew',3)
        set(gca,'xticklabel',table_head_show)
        grid on
        if ii==1
            title(['differential density: ', density{dd}])
        end
        if dd==1
            ylabel([row_name{ii},'(%)'])
        end
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf,[figurepath,'3_3_B'])
%     print([figurepath,'3_3_B'],'-depsc')
%% 3.4

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'FGN','cvxFGN','Song','Skrip'};
table_head_show = {'FGN','cvx-FGN','Songsiri, 17','Skripnikov, 19a'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'adaptive_formulation_S_result_K5'])
result.FGN = R;
load([resource_path,'adaptive_formulation_S_cvx_result_K5'])
result.cvxFGN = R;
load([resource_path,'adaptive_formulation_S_JSS_result_K5'])
result.Song = R;
load([resource_path,'skripS_result'])
result.Skrip = R;
M = zeros(5,4,2);
STD = zeros(5,4,2);
summary = zeros(5,4,2,100);
for jj=1:length(table_head)
    for ii=1:length(row_name)
        for dd=1:2
            M(ii,jj,dd) = mean(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            STD(ii,jj,dd) = std(result.(table_head{jj}).total.(row_name{ii})(dd,:));
            summary(ii,jj,dd,:) = result.(table_head{jj}).total.(row_name{ii})(dd,:);
            
        end
    end
end

% table
printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head_show{:},table_head_show{:}},row_name)

tt = tiledlayout(2,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
density = {'1%','5%'};
for ii=1:2
    for dd=1:2
        ARR = zeros(100,4);
        for jj=1:length(table_head)
            ARR(:,jj) = summary(ii,jj,dd,:);
        end
        nexttile;
        boxplot(100*ARR);
        set(findobj(gca,'type','line'),'linew',3)
        set(gca,'xticklabel',table_head_show)
        grid on
        if ii==1
            title(['differential density: ', density{dd}])
        end
        if dd==1
            ylabel([row_name{ii},'(%)'])
        end
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf,[figurepath,'3_4'])
%     print([figurepath,'3_4'],'-depsc')
%% 3.5 CGN

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'CGN','cvxCGN'};
table_head_show = {'CGN','cvx-CGN'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'adaptive_formulation_CT150_result_K5'])
result.CGN = R;
load([resource_path,'adaptive_formulation_CT150_cvx_result_K5'])
result.cvxCGN = R;
M = zeros(3,5,2);
STD = zeros(3,5,2);
summary = zeros(3,5,2,100);
type_acc = {'total','common','differential'};
for jj=1:length(table_head)
    for ii=1:length(row_name)
        for t=1:3
            dd=2;
            M(t,ii,jj) = mean(result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:));
            STD(t,ii,jj) = std(result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:));
            summary(t,ii,jj,:) = result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:);
        end
    end
end

% table
printtable([squeeze(M(1,:,:)*100) squeeze(M(2,:,:)*100) squeeze(M(3,:,:)*100)],[squeeze(STD(1,:,:)*100) squeeze(STD(2,:,:)*100) squeeze(STD(3,:,:)*100)],{table_head_show{:},table_head_show{:},table_head_show{:}},row_name)

tt = tiledlayout(2,1);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
density = {'1%','5%'};
    dd=2;
for ii=1:2
    nexttile;
    for jj=1:length(table_head)
        for tt=1:3
            ARR(:,tt) = summary(tt,ii,jj,:);
        end
        data_to_plot{1,jj} = 100*ARR;
    end
    
    a=boxplotGroup(data_to_plot, 'PrimaryLabels', {'CGN', 'cvx-CGN'}, ...
    'SecondaryLabels',{'total' 'common' 'differential'}, 'GroupLabelType', 'Vertical');
grid on
ylabel([row_name{ii},'(%)'])
        set(findobj(gca,'type','line'),'linew',3)
end
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf,[figurepath,'3_5_CGN'])
%     print([figurepath,'3_5_CGN'],'-depsc')
%% 3.5 DGN

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'DGN','cvxDGN'};
table_head_show = {'DGN','cvx-DGN'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'adaptive_formulation_DT150_result_K5'])
result.DGN = R;
load([resource_path,'adaptive_formulation_DT150_cvx_result_K5'])
result.cvxDGN = R;
M = zeros(3,5,2);
STD = zeros(3,5,2);
summary = zeros(3,5,2,100);
type_acc = {'total','common','differential'};
for jj=1:length(table_head)
    for ii=1:length(row_name)
        for t=1:3
            dd=2;
            M(t,ii,jj) = mean(result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:));
            STD(t,ii,jj) = std(result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:));
            summary(t,ii,jj,:) = result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:);
        end
    end
end

% table
printtable([squeeze(M(1,:,:)*100) squeeze(M(2,:,:)*100) squeeze(M(3,:,:)*100)],[squeeze(STD(1,:,:)*100) squeeze(STD(2,:,:)*100) squeeze(STD(3,:,:)*100)],{table_head_show{:},table_head_show{:},table_head_show{:}},row_name)

hh = tiledlayout(2,1);
hh.TileSpacing = 'compact';
hh.Padding = 'compact';
density = {'1%','5%'};
    dd=2;
for ii=1:2
    nexttile;
    
    for jj=1:length(table_head)
        for tt=1:3
            ARR(:,tt) = summary(tt,ii,jj,:);
        end
        data_to_plot{1,jj} = 100*ARR;
    end
    
    a=boxplotGroup(data_to_plot, 'PrimaryLabels', {'CGN', 'cvx-CGN'}, ...
    'SecondaryLabels',{'total' 'common' 'differential'}, 'GroupLabelType', 'Vertical');
grid on
ylabel([row_name{ii},'(%)'])

        set(findobj(gca,'type','line'),'linew',3)

%         set(gca,'xticklabel',table_head,'fontsize',20)
end
set(findall(gcf,'-property','FontSize'),'FontSize',24)

set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf,[figurepath,'3_5_DGN'])
%     print([figurepath,'3_5_DGN'],'-depsc')
%% 3.5 FGN

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'FGN','cvxFGN'};
table_head_show = {'FGN','cvx-FGN'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'adaptive_formulation_ST150_result_K5'])
result.FGN = R;
load([resource_path,'adaptive_formulation_ST150_cvx_result_K5'])
result.cvxFGN = R;
M = zeros(3,5,2);
STD = zeros(3,5,2);
summary = zeros(3,5,2,100);
type_acc = {'total','common','differential'};
for jj=1:length(table_head)
    for ii=1:length(row_name)
        for t=1:3
            dd=2;
            M(t,ii,jj) = mean(result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:));
            STD(t,ii,jj) = std(result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:));
            summary(t,ii,jj,:) = result.(table_head{jj}).(type_acc{t}).(row_name{ii})(dd,:);
        end
    end
end

% table
printtable([squeeze(M(1,:,:)*100) squeeze(M(2,:,:)*100) squeeze(M(3,:,:)*100)],[squeeze(STD(1,:,:)*100) squeeze(STD(2,:,:)*100) squeeze(STD(3,:,:)*100)],{table_head_show{:},table_head_show{:},table_head_show{:}},row_name)

tt = tiledlayout(2,1);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
density = {'1%','5%'};
    dd=2;
for ii=1:2
    nexttile;
    
    for jj=1:length(table_head)
        for tt=1:3
            ARR(:,tt) = summary(tt,ii,jj,:);
        end
        data_to_plot{1,jj} = 100*ARR;
    end
    
    a=boxplotGroup(data_to_plot, 'PrimaryLabels', {'CGN', 'cvx-CGN'}, ...
    'SecondaryLabels',{'total' 'common' 'differential'}, 'GroupLabelType', 'Vertical');
grid on
ylabel([row_name{ii},'(%)'])
        set(findobj(gca,'type','line'),'linew',3)
end
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf,[figurepath,'3_5_FGN'])
%     print([figurepath,'3_5_FGN'],'-depsc')