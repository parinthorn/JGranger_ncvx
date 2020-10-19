%% Generate data
clear all; close all;clc
% addpath('../lasso/');
% addpath('../lasso/admm');

n = 10; p = 3; N = 100; K = 4;
density = 0.5;
B = zeros(n,n*p,K);
H = zeros(n*p,N,K);
Y = zeros(n,N,K);

for k=1:K,
    [ind_nz,phis,y] = gen_sparseAR(n,p,1,density,N);
    [tmp,Hk,Yk] = gen_timeseries(phis,N+p);

    
    B(:,:,k) = reshape(phis,n,n*p); 
    H(:,:,k) = Hk;
    Y(:,:,k) = Yk;
end

[b,G,D,x] = veccoefmatgroup(Y,H,B,[n p K N]);

nn = n^2*p*K; nd = n^2*p*(K-1);

close all;
%% Solve fused lasso problem
a1 = 20; a2 = 50; 
rho = 100; % ADMM parameter

[x,history] = fusedlasso_ADMM(G, b, a1, a2, [n p K], rho,D);

%% Detect the zero pattern
X = reshape(x,p,K,n^2);
figure;
for k=1:K,
%     TMP1 = squeeze(all(abs(X(:,k,:)) == 0)); 
    TMP2 = squeeze(any(abs(X(:,k,:)) > 0));
%     IND_Z = find(TMP1);
    IND_NZ{k} = find(TMP2);
    subplot(2,K,k);plot_spy(IND_NZ{k},n);title(['# ',num2str(k)]);
end

IND_NZ_COMMON = intersect(IND_NZ{1},IND_NZ{2},'rows');
if K > 2,
    for k=3:K,
        IND_NZ_COMMON = intersect(IND_NZ_COMMON,IND_NZ{k},'rows');
    end
end
subplot(2,K,K+1);plot_spy(IND_NZ_COMMON,n);title('common pattern');

%% group lasso
figure;
for k=1:K,
    [bk,Gk] = veccoefmat(Y(:,:,k),H(:,:,k),B(:,:,k));
    [zk] = group_lasso(Gk, bk, a1, p*ones(n^2,1), rho, 1);
    Z = reshape(zk,p,n^2);
    TMP = any(abs(Z) > 0,1); TMP = TMP';
    IND_NZ_GL{k} = find(TMP);
    
    subplot(2,K,k);plot_spy(IND_NZ_GL{k},n);title(['# ',num2str(k)]);
end
IND_NZ_GL_COMMON = intersect(IND_NZ_GL{1},IND_NZ_GL{2},'rows');
if K > 2,
    for k=3:K,
        IND_NZ_GL_COMMON = intersect(IND_NZ_GL_COMMON,IND_NZ_GL{k},'rows');
    end
end
subplot(2,K,K+1);plot_spy(IND_NZ_GL_COMMON,n);title('common pattern');
