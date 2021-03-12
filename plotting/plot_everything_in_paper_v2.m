%% exp_CGN
clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'CGN','cvxCGN','Song','Greg'};
table_head_show = {'CGN','cvx-CGN','Song17C','Greg15'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'LLHcorrected_adaptive_formulation_C_result'])
result.CGN = R;
load([resource_path,'LLHcorrected_adaptive_formulation_C_cvx_result'])
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
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
    saveas(gcf,[figurepath,'exp_CGN'])
    print([figurepath,'exp_CGN'],'-depsc')
%% exp_DGN_A

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'DGN','cvxDGN','Song','Skrip'};
table_head_show = {'DGN','cvx-DGN','Song17D','Skrip19b'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'LLHcorrected_adaptive_formulation_D_result_K5'])
result.DGN = R;
load([resource_path,'LLHcorrected_adaptive_formulation_D_cvx_result_K5'])
result.cvxDGN = R;
load([resource_path,'LLHcorrected_adaptive_formulation_D_JSS_result_K5'])
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
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);

set(findall(gcf,'-property','FontSize'),'FontSize',28)
% 
    saveas(gcf,[figurepath,'exp_DGN_A'])
    print([figurepath,'exp_DGN_A'],'-depsc')
%% exp_DGN_B_xxx
clear
clc
clf
close all
toggle = 'total';
dd=2;
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'DGN','cvxDGN','Song','Skrip'};
table_head_show = {'DGN','cvx-DGN','Song17D','Skrip19b'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
K_list = {'K5','K50'};

load([resource_path,'LLHcorrected_adaptive_formulation_D_result_K5'])
result.DGN.K5 = R;
load([resource_path,'LLHcorrected_adaptive_formulation_D_cvx_result_K5'])
result.cvxDGN.K5 = R;
load([resource_path,'LLHcorrected_adaptive_formulation_D_JSS_result_K5'])
result.Song.K5 = R;
load([resource_path,'ResultSkripD_K5'])
result.Skrip.K5 = R;


load([resource_path,'LLHcorrected_adaptive_formulation_D_result_K50'])
result.DGN.K50 = R;
load([resource_path,'LLHcorrected_adaptive_formulation_D_cvx_result_K50'])
result.cvxDGN.K50 = R;
load([resource_path,'LLHcorrected_adaptive_formulation_D_JSS_result_K50'])
result.Song.K50 = R;
load([resource_path,'ResultSkripD_K50'])
result.Skrip.K50 = R;
M = zeros(5,4,2);
STD = zeros(5,4,2);
summary = zeros(5,4,2,100);

for jj=1:length(table_head)
    
    for ii=1:length(row_name)
        for kk=1:length(K_list)
            dd=2;
            M(ii,jj,kk) = mean(result.(table_head{jj}).(K_list{kk}).(toggle).(row_name{ii})(dd,:));
            STD(ii,jj,kk) = std(result.(table_head{jj}).(K_list{kk}).(toggle).(row_name{ii})(dd,:));
            summary(ii,jj,kk,:) = result.(table_head{jj}).(K_list{kk}).(toggle).(row_name{ii})(dd,:);
            
        end
    end
end
M = permute(M,[1,3,2]);
STD = permute(STD,[1,3,2]);

Mt = reshape(M,[5,8]);
STDt = reshape(STD,[5,8]);
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
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
    saveas(gcf,[figurepath,'exp_DGN_B_',toggle])
    print([figurepath,'exp_DGN_B_',toggle],'-depsc')
%% exp_FGN

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'FGN','cvxFGN','Song','Skrip'};
table_head_show = {'FGN','cvx-FGN','Song17F','Skrip19a'};
row_name = {'F1','FPR','TPR','ACC','MCC'};
load([resource_path,'LLHcorrected_adaptive_formulation_S_result_K5'])
result.FGN = R;
load([resource_path,'LLHcorrected_adaptive_formulation_S_cvx_result_K5'])
result.cvxFGN = R;
load([resource_path,'LLHcorrected_adaptive_formulation_S_JSS_result_K5'])
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
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
    saveas(gcf,[figurepath,'exp_FGN'])
    print([figurepath,'exp_FGN'],'-depsc')
%% exp_cvxcompare
clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
table_head = {'C','D','F'};

row_name = {'F1','FPR','TPR','ACC','MCC'};

load([resource_path,'LLHcorrected_adaptive_formulation_CT150_result_K5'])
result.C.ncvx = R;
load([resource_path,'LLHcorrected_adaptive_formulation_CT150_cvx_result_K5'])
result.C.cvx = R;

load([resource_path,'LLHcorrected_adaptive_formulation_DT150_result_K5'])
result.D.ncvx = R;
load([resource_path,'LLHcorrected_adaptive_formulation_DT150_cvx_result_K5'])
result.D.cvx = R;

load([resource_path,'LLHcorrected_adaptive_formulation_ST150_result_K5'])
result.F.ncvx = R;
load([resource_path,'LLHcorrected_adaptive_formulation_ST150_cvx_result_K5'])
result.F.cvx = R;

M = zeros(3,2,5); % [C,D,F] x [ncvx, cvx] x [F1, FPR, TPR, ACC , MCC]
STD = zeros(3,2,5);
summary = zeros(3,2,5,100);
type_acc = {'ncvx','cvx'};
for jj=1:length(table_head)
    for ii=1:length(row_name)
        for t=1:2
            dd=2;
            if jj==1
                toggle= 'common';
            else
                toggle = 'total';
            end
            M(jj,t,ii) = mean(result.(table_head{jj}).(type_acc{t}).(toggle).(row_name{ii})(dd,:));
            STD(jj,t,ii) = std(result.(table_head{jj}).(type_acc{t}).(toggle).(row_name{ii})(dd,:));
            summary(jj,t,ii,:) = result.(table_head{jj}).(type_acc{t}).(toggle).(row_name{ii})(dd,:);
        end
    end
end
M = permute(M,[3,2,1]);
STD = permute(STD,[3,2,1]);
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
pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
    saveas(gcf,[figurepath,'exp_cvxcompare'])
    print([figurepath,'exp_cvxcompare'],'-depsc')
%% exp_CGN_ROC

clear
clc
clf
close all
type_list = {'total','common','differential'};
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
method_path = {'LLHcorrected_adaptive_formulation_C_ALL_RESULT.mat','LLHcorrected_adaptive_formulation_C_cvx_ALL_RESULT.mat'};
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
    saveas(gcf,[figurepath,'exp_CGN_ROC'])
    print([figurepath,'exp_CGN_ROC'],'-depsc')
%% exp_DGN_supp

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
load([resource_path,'LLHcorrected_adaptive_formulation_D_cvx_ALL_RESULT_K50.mat'])
sample = 74;
tt=tiledlayout(1,2,'padding','compact','tilespacing','none');
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
    
    set(gca,'FontSize',28,'xaxisLocation','top')
    
%     set(gca,)
end

end
rr=colorbar('FontSize',18);
set(get(rr,'label'),'string','F1 (%)');
pp = get(0, 'Screensize');
% pp(3) = pp(3)*0.75;
set(gcf, 'Position', pp);
    saveas(gcf,[figurepath,'exp_DGN_supp'])
    print([figurepath,'exp_DGN_supp'],'-depsc')
    
%% exp_FGN_supp

clear
clc
clf
close all
figurepath = './plotting/figures/';
resource_path = './experiment/result_to_plot/';
load([resource_path,'LLHcorrected_adaptive_formulation_S_cvx_ALL_RESULT_K5.mat'])
sample = 74;
tt=tiledlayout(1,2,'padding','compact','tilespacing','none');
type_list = {'total','common','differential'};
type_list_show = {'cvx-FGN, total bestcase'};
score_type = 'F1';
text_label = {'Differential  density 1%','Differential density 5%'};
for dd=1:2




    % for ii=2:2
    ii =1;
    GridF1_avg = zeros(30,30);
    F1_avg = 0;
    nexttile;
    for sample=1:100
        sample_acc = ALL_RESULT(dd,sample);
        metrics = [sample_acc.model_acc]; metrics= [metrics.(type_list{ii})]; metrics=[metrics.F1];
        
        bestcase = max(metrics);
        metrics = reshape(metrics,[30,30]);
        GridF1_avg=GridF1_avg+metrics/100;
        F1_avg=F1_avg+bestcase/100;
    end

    imagesc(100*GridF1_avg)
%     grid on
    text_show = sprintf([type_list_show{1},' ',score_type,' score: %2.1f'],100*F1_avg);
%     title(text_show)
    axis('square')
    colormap((1-gray).^0.4)
    caxis([0,100])
    set(gca,'xticklabel',[],'yticklabel',[])
    if dd==1
%                 colorbar('Location','westoutside')
        xlabel('$\lambda_{2}\rightarrow$','Interpreter','latex')
        
        set(gca,'xaxisLocation','top')
        ylabel('$\leftarrow \lambda_{1}$','Interpreter','latex')
        
    else
        xlabel('$\lambda_{2}\rightarrow$','Interpreter','latex')
        
        set(gca,'xaxisLocation','top')
    end
    
    if ii==1
        title(text_label{dd})
        
    end
    
    set(gca,'FontSize',28)
% end

end

rr=colorbar('FontSize',18);
set(get(rr,'label'),'string','F1 (%)');
% rr.Position(4)=rr.Position(4)-0.005;
% rr.Position(2)=rr.Position(2)+0.005/2;
set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,[figurepath,'exp_FGN_supp'])
    print([figurepath,'exp_FGN_supp'],'-depsc')