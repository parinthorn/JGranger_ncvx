%% This script is for printing result

%% 3.2: Common network comparison (CommonGrangerNet, convex CommonGrangerNet, JSS, Magda Gregorova)
clear
clc
clf
close all
resource_path = './experiment/result_to_plot/';
table_head = {'CGN','cvxCGN','Song','Greg'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'formulation_C_eBICresult'])
result.CGN = R;
load([resource_path,'formulation_C_cvx_eBICresult'])
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
printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head{:},table_head{:}},row_name)

tt = tiledlayout(2,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
density = {'10%','20%'};

for ii=1:2
    for dd=1:2
        ARR = zeros(100,4);
        for jj=1:length(table_head)
            ARR(:,jj) = summary(ii,jj,dd,:);
        end
        nexttile;
        boxplot(ARR);
        set(findobj(gca,'type','line'),'linew',2)
        set(gca,'xticklabel',table_head,'fontsize',20)
        grid on
        if ii==1
            title(['common density: ', density{dd}])
        
        end
        if dd==1
            ylabel(row_name{ii})
        end
    end
    
    
end
set(gcf, 'Position', get(0, 'Screensize'));

%% 3.3 A

clear
clc
clf
close all
resource_path = './experiment/result_to_plot/';
table_head = {'CGN','cvxCGN','Song','Skrip'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'adaptive_formulation_D_result_K5'])
result.CGN = R;
load([resource_path,'formulation_D_cvx_eBICresult'])
result.cvxCGN = R;
load([resource_path,'adaptive_formulation_D_JSS_result_K5'])
result.Song = R;
load([resource_path,'skripD_result'])
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
printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head{:},table_head{:}},row_name)

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
        boxplot(ARR);
        set(findobj(gca,'type','line'),'linew',2)
        set(gca,'xticklabel',table_head,'fontsize',20)
        grid on
        if ii==1
            title(['common density: ', density{dd}])
        end
        if dd==1
            ylabel(row_name{ii})
        end
    end    
end
set(gcf, 'Position', get(0, 'Screensize'));




%% 3.3 B



%% 3.4

clear
clc
clf
close all
resource_path = './experiment/result_to_plot/';
table_head = {'CGN','cvxCGN','Song','Skrip'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'formulation_S_eBICresult'])
result.CGN = R;
load([resource_path,'formulation_S_cvx_eBICresult'])
result.cvxCGN = R;
load([resource_path,'adaptive_formulation_S_JSS_result_K5'])
result.Song = R;
load([resource_path,'skripD_result'])
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
printtable([M(:,:,1)*100 M(:,:,2)*100],[STD(:,:,1)*100 STD(:,:,2)*100],{table_head{:},table_head{:}},row_name)

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
        boxplot(ARR);
        set(findobj(gca,'type','line'),'linew',2)
        set(gca,'xticklabel',table_head,'fontsize',20)
        grid on
        if ii==1
            title(['common density: ', density{dd}])
        end
        if dd==1
            ylabel(row_name{ii})
        end
    end    
end
set(gcf, 'Position', get(0, 'Screensize'));

%% 3.5