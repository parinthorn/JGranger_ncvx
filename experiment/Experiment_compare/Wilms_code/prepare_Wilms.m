%% Load data
clear
clc
D = readtable('C:\Users\CU_EE_LAB408\Downloads\Code\wilms_data.csv');
D = table2array(D);
D = D(:,2:end);
n=15;K=15;T=76;
Y = reshape(D,[T,n,K]);
y = permute(Y,[2,1,3]);
GridSize = 30;
%% Vectorizing model
[n,Num,K] = size(y);

p=2;
[P,~] = offdiagJSS(n,p,K);
Dmat = diffmat(n,p,K);
D = sparse(Dmat*P);

Lambda = logspace(-3,0,GridSize);
tmp(K)= struct();
disp('Generating H matrix')
for kk=1:K
    [tmp(kk).H,tmp(kk).Y] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc,~] = veccoefmatgroup(reshape([tmp.Y],[n,Num-p,K]),reshape([tmp.H],[n*p,Num-p,K]),[],[n,p,K,Num-p]);
disp('calculating Lambda max')
Lmax = lambdamax_grouplasso(gc,yc,[n ,p ,K]);
Lambda = [0 Lambda*Lmax];
xLS = gc\yc;
a1 = Lambda(1);
a2 = Lambda(25);
pp=2;qq=0.5;
x_S_cvx = zeros(n^2*p*K,GridSize,GridSize);
x_D_cvx = zeros(n^2*p*K,GridSize,GridSize);
x_S = zeros(n^2*p*K,GridSize,GridSize);
x_D = zeros(n^2*p*K,GridSize,GridSize);
A_S_cvx = zeros(n,n,p,K,GridSize,GridSize);
A_D_cvx = zeros(n,n,p,K,GridSize,GridSize);
A_S = zeros(n,n,p,K,GridSize,GridSize);
A_D = zeros(n,n,p,K,GridSize,GridSize);
for ii=1:30
    for jj=1:30
        disp('Intitialize formulation S')
        [x_S_cvx(:,ii,jj), ~, ~,history_S_cvx] = spectral_ADMM_heuristic(gc, yc, P,D,a1, a2,pp,1, [n,p,K,p,p], 1,0.1,2,xLS);
        [x_S(:,ii,jj), ~,~, history_S] = spectral_ADMM_heuristic(gc, yc, P,D,a1, a2,pp,qq, [n,p,K,p,p], 1,0.1,100,x_S_cvx(:,ii,jj));
        disp('Intitialize formulation D')
        [x_D_cvx(:,ii,jj), ~, ~,history_D_cvx] = spectral_ADMM_heuristic(gc, yc, P,P,a1, a2,pp,1, [n,p,K,p,p*K], 1,0.1,2,xLS);
        [x_D(:,ii,jj), ~,L2x_D, history_D] = spectral_ADMM_heuristic(gc, yc, P,P,a1, a2,pp,qq, [n,p,K,p,p*K], 1,0.1,100,x_D_cvx(:,ii,jj));
        
        
        A_S(:,:,:,:,ii,jj) = devect(x_S(:,ii,jj),n,p,K);
        A_D(:,:,:,:,ii,jj) = devect(x_D(:,ii,jj),n,p,K);
        A_S_cvx(:,:,:,:,ii,jj) = devect(x_S_cvx(:,ii,jj),n,p,K);
        A_D_cvx(:,:,:,:,ii,jj) = devect(x_D_cvx(:,ii,jj),n,p,K);
    end
end

