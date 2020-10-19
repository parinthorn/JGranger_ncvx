function M = formulation_D(y,P,varargin)
%% This program estimates Ground truth with their original model
% Input are
%             y : 3D array [n,Num,K] which is dimension of timeseries,
%       timepoints, # models respectively
%             P : Off-diagonal projection matrix
% (optional)  p : lag-order
%      GridSize : dimension of a regularization grid
% Output is
%       M : structure containing
%           M.x_est
[n,Num,K] = size(y);
len_varargin = length(varargin);
if isempty(varargin)
  p=2;
  GridSize = 30;
elseif len_varargin==1
  p = varargin{1};
  GridSize = 30;
elseif len_varargin ==2
  p = varargin{1};
  GridSize = varargin{2};
else
  error('must be atmost 5 input')
end
Lambda = logspace(-3,0,GridSize);
tmp(K)= struct();
disp('Generating H matrix')
for kk=1:K
    [tmp(kk).H,tmp(kk).Y] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc,~] = vectorize_VAR(reshape([tmp.Y],[n,Num-p,K]),reshape([tmp.H],[n*p,Num-p,K]));
disp('calculating Lambda max')
Lmax = lambdamax_grouplasso(gc,yc,[n ,p ,K]);
Lambda = Lambda*Lmax;
xLS = gc\yc;
if Num>n*p
  init_cvx = 0;
else
  init_cvx = 1;
end
M.lambda_crit = Lmax;
M.lambda_range = [Lambda(1) Lambda(end)];
M.GridSize = GridSize;
M.flag = zeros(GridSize);
t1 = tic;
for ii=1:GridSize
    a1 = Lambda(ii);
    A_reg = zeros(n,n,p,K,1,GridSize);
    A = zeros(n,n,p,K,1,GridSize);
    parfor jj=1:GridSize
        fprintf('Grid : (%d,%d)/(%d, %d) \n',ii,jj,GridSize,GridSize)
        if init_cvx
          [x0, ~, ~] = spectral_ADMM(gc, yc, P,P,a1, Lambda(jj),2,1, [n,p,K,p,p*K], 1,0.6,2,xLS);
        else
          x0 = xLS;
        end
        [x_est, ~,~, history] = spectral_ADMM(gc, yc, P,P,a1, Lambda(jj),2,0.5, [n,p,K,p,p*K], 1,0.1,100,x0);
        A_reg_tmp = devect(full(x_est),n,p,K); % convert to (n,n,p,K) format
        A_reg(:,:,:,:,1,jj) = A_reg_tmp; % this is for arranging result into parfor format
        [x_cls,ls_flag(1,jj)] = constrained_LS_D(gc,yc,find(x_est));
        A_cls =devect(full(x_cls),n,p,K);
        A(:,:,:,:,1,jj) = A_cls;
        score(1,jj) = model_selection(yc,A_cls)
        nz_ind = cell(K,1);
        diag_ind=1:n+1:n^2;
        for kk=1:K
          tmp_nz_ind{kk} = setdiff(squeeze(A_reg_tmp(:,:,1,:)),diag_ind)
        end
        nz_ind(1,jj) = tmp_nz_ind;
        [ind_common(1,jj),ind_differential(1,jj)] = split_common_diff(tmp_nz_ind,[n,p,K]); % find common and differential off-diagonal nonzero index
        flag(1,jj) = history.flag;
        if flag(1,jj) ==-1
            fprintf('max iteration exceed at grid (%d,%d)\n',ii,jj)
        end
    end
    M.stat.score(ii,:) = score;
    M.A_reg(:,:,:,:,ii,:) = A_reg;
    M.A(:,:,:,:,ii,:) = A;
    M.nz_ind(ii,:) = nz_ind;
    M.flag(ii,:) = flag;
end
M.time = toc(t1);
end
