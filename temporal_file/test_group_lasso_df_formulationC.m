%% This script iterates through model and evaluate new model criteria score
% FORMULATION C
clear
clc
load('D:\JGranger_ncvx\data_compare\model_K5_p1.mat') % gt model


realization=3;
dd = 2;
cd = [3];
mname = {'10','20'};
n=20;p=1;K=5;
T = 100;
for cd_index=1:length(cd)
    for realization=1:1
        GTmodel = E{2,cd(cd_index),dd,realization};
        load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\formulationC_',mname{cd_index},'percent_',int2str(realization),'.mat'])
        GridSize = M.GridSize;
        for ii=1:1
            for jj=1:GridSize
                model = M.model(ii,jj);
                
                idx = efficient_vect([n,p,K]);
                x_reg = model.A_reg(idx);
                x_cls = model.A(idx);
                tmp_reg = reshape(x_reg,[p*K,n^2]);
                tmp_ls = reshape(x_cls,[p*K,n^2]);
                tmp_ls =sqrt(sum(tmp_ls.^2,1));
                tmp_reg = sqrt(sum(tmp_reg.^2,1));
                df(ii,jj) = length(find(tmp_reg))+(p*K-1)*sum(tmp_reg./tmp_ls,'omitnan');
                df_lasso(ii,jj) = length(find(x_reg));
                
                
                bic(ii,jj) = log(T-p)*df(ii,jj)+model.stat.model_selection_score.bic-log(T-p)*df_lasso(ii,jj);
                bic_lasso(ii,jj) = model.stat.model_selection_score.bic;
            end
            [~,M.index.bic_new] = min(bic);
%             save(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_C_magda\new_bic_formulationC_',mname{cd},'percent_',int2str(realization),'.mat'])

        end
        
        
        
    end
end

% GridSize = M.GridSize;
% y = sim_VAR(GTmodel.A,T,1,GTmodel.seed,0);
% H = zeros(n*p,T-p,K);
% Y = zeros(n,T-p,K);
% disp('Generating H matrix')
% for kk=1:K
%     [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
% end
%%

%
clf;close all;
figure;
subplot(1,2,1);imagesc(df);subplot(1,2,2);imagesc(df_lasso)
figure;
subplot(1,2,1);imagesc(bic);subplot(1,2,2);imagesc(bic_lasso)
plot_group_GC(GTmodel.GC(:,:,1:5))
sgtitle('ground-truth model')

plot_group_GC(M.model(M.index.bic).GC(:,:,1:5))
sgtitle('Degrees of freedom: Lasso')

[t1,t2] = min(bic(:));
plot_group_GC(M.model(t2).GC(:,:,1:5))
sgtitle('Degrees of freedom: Group-lasso block-sized pK')