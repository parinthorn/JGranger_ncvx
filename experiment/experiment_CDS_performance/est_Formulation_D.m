function M = est_Formulation_D(y,A_true,P)
%% This program estimates Ground truth with their original model
% Input are
%       y : 3D array [n,Num,K] which is dimension of timeseries,
%       timepoints, # models respectively
%       A : true VAR parameters [n,n,p,K]
%       P : Off-diagonal projection matrix
%       f : model diff density 1 -> 5%, 2-> 20%
% Output is
%       M : structure containing
%           M.stat
%           M.x_true
%           M.x_est
[n,Num,K] = size(y);
p=2;
GridSize = 30;
Lambda = logspace(-3,0,GridSize);
tmp(K)= struct();
disp('Generating H matrix')
for kk=1:K
    [tmp(kk).H,tmp(kk).Y] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc,~] = veccoefmatgroup(reshape([tmp.Y],[n,Num-p,K]),reshape([tmp.H],[n*p,Num-p,K]),A_true,[n,p,K,Num-p]);
disp('calculating Lambda max')
Lmax = lambdamax_grouplasso(gc,yc,[n ,p ,K]);
Lambda = Lambda*Lmax;
xLS = gc\yc;
M.stat.Lambda_max = Lmax;
M.stat.consistency = zeros(2,K,GridSize);
M.stat.bias = zeros(GridSize);
M.flag = zeros(GridSize);
t1 = tic;
% M.Lipschitz = eigs(gc'*gc,1);
for ii=1:GridSize
    a1 = Lambda(ii);
    tmpA = A_true(:);
    parfor jj=1:GridSize
        fprintf('Grid : (%d,%d)/(%d, %d) \n',ii,jj,GridSize,GridSize)
        [x0, ~, ~] = spectral_ADMM(gc, yc, P,P,a1, Lambda(jj),2,1, [n,p,K,p,p*K], 1,0.6,2,xLS);
        [x_est, ~,~, history] = spectral_ADMM(gc, yc, P,P,a1, Lambda(jj),2,0.5, [n,p,K,p,p*K], 1,0.1,100,x0);
        A_est = devect(full(x_est),n,p,K);
        A_cvx = devect(full(x0),n,p,K);
        A_Lpq(:,:,:,:,1,jj) = A_est;
        A_L21(:,:,:,:,1,jj) = A_cvx;
        nz_ind{1,jj} = find(A_est);
        flag(1,jj) = history.flag;
        if flag(1,jj) ==-1
            fprintf('max iteration exceed at grid (%d,%d)\n',ii,jj)
        end
        bias(1,jj) = norm(A_est(:)-tmpA)/norm(tmpA);
    end
    M.A_Lpq(:,:,:,:,ii,:) = A_Lpq;
    M.A_cvx(:,:,:,:,ii,:) = A_L21;
    M.nz_ind{ii,:} = nz_ind;
    M.flag(ii,:) = flag;
    M.stat.bias(ii,:) = bias;
end
M.time = toc(t1);
end