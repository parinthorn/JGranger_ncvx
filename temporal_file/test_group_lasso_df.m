%% This script iterates through model and evaluate new model criteria score


clear
clc
load('.\data_compare\model_K50_p1.mat') % gt model
realization=1;
GTmodel = E{2,3,1,realization};
[n,p,K] = feval(@(x) x{:}, num2cell(GTmodel.dim));
T = 100;

y_sim = sim_VAR(GTmodel.A,T,1,GTmodel.seed,0);
H = zeros(n*p,T-p,K);
Y = zeros(n,T-p,K);
disp('Generating H matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y_sim(:,:,kk),p);
end
[yc,gc] = vectorize_VAR(Y,H,[n,p,K,T]);


n=20;p=1;K=50;
T = 100;
load(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\result_formulationD_1percent_lag1_K50_',int2str(realization),'.mat'])
GridSize = M.GridSize;
% y = sim_VAR(GTmodel.A,T,1,GTmodel.seed,0);
% H = zeros(n*p,T-p,K);
% Y = zeros(n,T-p,K);
% disp('Generating H matrix')
% for kk=1:K
%     [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
% end


%%
for ii=1:GridSize
    for jj=1:GridSize
        model = M.model(ii,jj);
        BLOCK_SIZE = p*K;
        idx = efficient_vect([n,p,K]);
        x_reg = model.A_reg(idx);
        x_cls = model.A(idx);
        for rm_ind=1:n^2*p*K
            gc_rm = gc; gc_rm(:,rm_ind) = [];
            x_rm =x_reg; x_rm(rm_ind) = [];
            y_rm = yc-gc_rm*x_rm;
            x_star(rm_ind,1) = gc(:,rm_ind)'*y_rm/sum(gc(:,rm_ind).^2);
        end
        
        tmp_reg = reshape(x_reg,[BLOCK_SIZE,n^2*K/BLOCK_SIZE]);
        tmp_ls = reshape(x_cls,[BLOCK_SIZE,n^2*K/BLOCK_SIZE]);
        tmp_ls =sqrt(sum(tmp_ls.^2,1));
        tmp_reg = sqrt(sum(tmp_reg.^2,1));
        df(ii,jj) = length(find(tmp_reg))+(BLOCK_SIZE-1)*sum(tmp_reg./tmp_ls,'omitnan');
        df_lasso(ii,jj) = length(find(x_reg));
        bic(ii,jj) = log(T-p)*df(ii,jj)+model.stat.model_selection_score.L;
        bic_lasso(ii,jj) = model.stat.model_selection_score.bic;
        df_JSS(ii,jj) = sum(x_reg./x_star);
        %         disp(bic_lasso(ii,jj)-log(T-p)*df_lasso(ii,jj)+log(T-p)*df(ii,jj)-bic(ii,jj))
    end
end
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