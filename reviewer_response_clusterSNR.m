%% SNR results clustering

% The setting will be SNR 0-6, 6-10, 10++ [low, mid, high]

%% 
% setting 2 diff density, 3 formulation [only convex]

clear
clc
inpath = './experiment/model_parameters/';
performance_path = './results2plot/';

file_list = {'CGN_CVX_result.mat', 'DGN_CVX_result_K5.mat', 'FGN_CVX_result_K5.mat'};
         
load([inpath,'model_K5_p1.mat']) % load ground-truth to compute snr
model_type = [2, 2, 3];
         
K=5;
p=1;
n=20;
sigmasq=1;



tt = tiledlayout(2, 3);
tt.Padding = 'compact';
tt.TileSpacing = 'compact';
tile_index = [1 2 3;4 5 6];
for f_type = 1:3
    
    
    realization = 100;
    snr_label = zeros(realization, 1);
    file_name = file_list{f_type};
    load([performance_path,file_name])
    
    
    
    for density_type = 1:2
        
        for jj=1:realization
            fprintf('(%d)\n',jj)
            
            if f_type == 1
                ind1 = density_type + 2;
                ind2 = 2;
                toggle = 'common';
            else
                ind1 = 3;
                ind2 = density_type;
                toggle = 'total';
            end
                
            
            
            GTmodel = E{model_type(f_type), ind1, ind2,jj};
            trSy_total = 0;
            trSe_total = 0;
            for kk=1:K
                Ak = GTmodel.A(:,:,:,kk); % the matrices in a cell array.
                [~, trSy, trSe] = compute_ar_SNR(Ak, sigmasq, [n, p]);
                trSy_total = trSy_total+trSy;
                trSe_total = trSe_total+trSe;
                
            end
            snr_label(jj) = 10*log10(trSy_total/trSe_total);
            fprintf('snr: %.2f\n', snr_label(jj))
        end
        
        
        score.F1 = R.(toggle).F1(density_type,:);
        group_index = zeros(1, realization);
        for ii=1:realization
            if snr_label(ii)<5
                group_index(ii) = 1;
            elseif snr_label(ii)<10
                group_index(ii) = 2;
            else
                group_index(ii) = 3;
            end
        end
        
        nexttile(tile_index(density_type, f_type));
        boxplot(score.F1, group_index)
        set(gca, 'xticklabel', {'snr<5','5<=snr<10','snr>=10'})
        ylabel('F1 score')
%         title(name_list{density_type, f_type})
        
        
    end
end
% set(findall(gcf,'-property','FontSize'),'FontSize',28)
pp = get(0, 'Screensize');
pp(3) = pp(3)*0.9;
pp(4) = pp(4)*0.7;
set(gcf, 'Position', pp);



%% cvx vs ncvx
clear
clc
inpath = './experiment/model_parameters/';
performance_path = './results2plot/';
file_list = {'T150_CGN_CVX_result_K5.mat', 'T150_DGN_CVX_result_K5.mat', 'T150_FGN_CVX_result_K5.mat'; ...
             'T150_CGN_result_K5.mat', 'T150_DGN_result_K5.mat', 'T150_FGN_result_K5.mat'};
         
name_list = {'cvx-CGN', 'cvx-DGN', 'cvx-FGN'; ...
                'CGN', 'DGN', 'FGN'};

K=5;
p = 3;
n = 20;
sigmasq = 1;
load([inpath,'compare_convex_model_K',int2str(K),'_p',int2str(p)]) % load ground-truth to compute snr
model_type = [2, 2, 3];

tt = tiledlayout(2, 3);
tt.Padding = 'compact';
tt.TileSpacing = 'compact';
tile_index = [1 2 3;4 5 6];
for f_type = 1:3
    
    
    realization = 100;
    snr_label = zeros(realization, 1);
    for jj=1:realization
        fprintf('(%d)\n',jj)
        GTmodel = E{model_type(f_type),jj};
        trSy_total = 0;
        trSe_total = 0;
        for kk=1:K
            Ak = GTmodel.A(:,:,:,kk); % the matrices in a cell array.
            [~, trSy, trSe] = compute_ar_SNR(Ak, sigmasq, [n, p]);
            trSy_total = trSy_total+trSy;
            trSe_total = trSe_total+trSe;
            
        end
        snr_label(jj) = 10*log10(trSy_total/trSe_total);
        fprintf('snr: %.2f\n', snr_label(jj))
    end
    
    for reg_type = 1:2
        
        file_name = file_list{reg_type, f_type};
        load([performance_path,file_name])
        score.bias = R.bias(2,:);
        group_index = zeros(1, realization);
        for ii=1:realization
            if snr_label(ii)<10
                group_index(ii) = 1;
            elseif snr_label(ii)<15
                group_index(ii) = 2;
            elseif snr_label(ii)<20
                group_index(ii) = 3;
            elseif snr_label(ii)<25
                group_index(ii) = 4;
            else
                group_index(ii) = 5;
            end
        end
        
        mean1 = mean(score.bias(group_index==1));
        mean2 = mean(score.bias(group_index==2));
        mean3 = mean(score.bias(group_index==3));
        mean4 = mean(score.bias(group_index==4));
        mean5 = mean(score.bias(group_index==5));
        
        siz1 = sum(group_index==1);
        siz2 = sum(group_index==2);
        siz3 = sum(group_index==3);
        siz4 = sum(group_index==4);
        siz5 = sum(group_index==5);
        
        nexttile(tile_index(reg_type, f_type));
        yyaxis left
        
        boxplot(score.bias, group_index)
        ylabel('bias')
        hold on
        
        b1 = plot([mean1 mean2 mean3 mean4 mean5], '-r', 'Linewidth', 2);
        yyaxis right
        ylabel('sample')
        b2 = bar([siz1 siz2 siz3 siz4 siz5]);
        legend('mean', 'sample size')
        hold off
        set(gca,'children',flipud(get(gca,'children')))
        b2.FaceAlpha = 0.5;

        
%         uistack(b1,'back')
        
        
        set(gca, 'xticklabel', {'snr<10','10<=snr<15','15<=snr<20','20<=snr<25','snr>=25'})
        
        title(name_list{reg_type, f_type})
        
        
        
    end
end
% set(findall(gcf,'-property','FontSize'),'FontSize',28)

pp = get(0, 'Screensize');
pp(3) = pp(3)*0.9;
pp(4) = pp(4)*0.7;
set(gcf, 'Position', pp);
figurepath = './results2plot/figures/';
exportgraphics(gcf,[figurepath,'reviewer_response_clusterSNR.png'])
exportgraphics(gcf,[figurepath,'reviewer_response_clusterSNR.eps'])






% tt = tiledlayout(1,1);


% group_index = (snr_label>20);


%%

ii=2;
tt = tiledlayout(1,1);
nexttile;
boxplot(R.total.F1(ii,:), ((snr_label(:)>20) +(snr_label(:)<15)));
%% FGN

%%