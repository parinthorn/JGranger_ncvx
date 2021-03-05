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
set(findall(gcf,'-property','FontSize'),'FontSize',28)
set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,[figurepath,'3_5_CGN'])
    print([figurepath,'3_5_CGN'],'-depsc')
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
    
    a=boxplotGroup(data_to_plot, 'PrimaryLabels', {'DGN', 'cvx-DGN'}, ...
        'SecondaryLabels',{'total' 'common' 'differential'}, 'GroupLabelType', 'Vertical');
    grid on
    ylabel([row_name{ii},'(%)'])
    
    set(findobj(gca,'type','line'),'linew',3)
    
    %         set(gca,'xticklabel',table_head,'fontsize',20)
end
set(findall(gcf,'-property','FontSize'),'FontSize',28)

set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,[figurepath,'3_5_DGN'])
    print([figurepath,'3_5_DGN'],'-depsc')
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
    
    a=boxplotGroup(data_to_plot, 'PrimaryLabels', {'FGN', 'cvx-FGN'}, ...
        'SecondaryLabels',{'total' 'common' 'differential'}, 'GroupLabelType', 'Vertical');
    grid on
    ylabel([row_name{ii},'(%)'])
    set(findobj(gca,'type','line'),'linew',3)
end
set(findall(gcf,'-property','FontSize'),'FontSize',28)
set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,[figurepath,'3_5_FGN'])
    print([figurepath,'3_5_FGN'],'-depsc')