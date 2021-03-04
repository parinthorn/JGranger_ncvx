function M = test_cvxformulation_C(y,varargin)
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
[n,T,K] = size(y);
len_varargin = length(varargin);
% toggle = 'static';
% toggle = 'adaptive_P';
% toggle = 'adaptive_L';
if isempty(varargin)
  p=1;
  GridSize = 30;
  toggle = 'adaptive_L';
elseif len_varargin==1
  p = varargin{1};
  GridSize = 30;
  toggle = 'adaptive_L';
elseif len_varargin ==2
  p = varargin{1};
  GridSize = varargin{2};
  toggle = 'adaptive_L';
elseif len_varargin ==3
  p = varargin{1};
  GridSize = varargin{2};
  toggle = varargin{3};
else
  error('must be atmost 4 input')
end
% Lambda = logspace(-6,0,GridSize);
% Lambda = logspace(-3.5,0,GridSize);
H = zeros(n*p,T-p,K);
Y = zeros(n,T-p,K);
disp('Generating H matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc] = vectorize_VAR(Y,H,[n,p,K,T]);
xLS = gc\yc;
% xLS = x_cheat;
if T>n*p
  init_cvx = 0;
else
  init_cvx = 1;
end
ALG_PARAMETER.PRINT_RESULT=0;
ALG_PARAMETER.IS_ADAPTIVE =1;
ALG_PARAMETER.dim = [n,p,K,p,p*K];
ALG_PARAMETER.rho_init = 1;
ALG_PARAMETER.epscor = 0.1;
ALG_PARAMETER.Ts = 50;
ALG_PARAMETER.is_chol = 1;
ALG_PARAMETER.multiplier = 2;
ALG_PARAMETER.toggle = 'formulationD';
ALG_PARAMETER.gamma = 1; % for adaptive case
ALG_PARAMETER.is_spectral = 1;
disp('calculating Lambda max')
qq=1; %non-convex case
[Lambda_1,Lambda_2,opt] = grid_generation(gc,yc,GridSize,ALG_PARAMETER,qq,toggle);
ALG_PARAMETER.L1 = opt.L1;
ALG_PARAMETER.L2 = opt.L2;
M.GridSize = GridSize;
M.flag = zeros(GridSize);
if isvector(Lambda_1)
    M.lambda2_crit = Lambda_2(end);
    M.lambda2_range = [Lambda_2(1) Lambda_2(end)];
    M.lambda1_crit = Lambda_1(end);
    M.lambda1_range = [Lambda_1(1) Lambda_1(end)];
else
    M.lambda1 = Lambda_1;
    M.lambda2 = Lambda_2;
end

t1 = tic;


    A_reg = zeros(n,n,p,K,GridSize);
    A = zeros(n,n,p,K,GridSize);
    ls_flag = zeros(1,GridSize);
    ind_common = cell(1,GridSize);
    ind_differential = cell(1,GridSize);
    flag = zeros(1,GridSize);
    ind = cell(1,GridSize);
parfor ii=1:GridSize
    a = Lambda_2(:,ii);
        fprintf('Grid : (%d)/(%d) \n',ii,GridSize)
        if init_cvx
            cvx_param = ALG_PARAMETER;
            cvx_param.Ts = 2;
          [x0, ~, ~] = spectral_ADMM_adaptive(gc, yc,0, a,2,1, cvx_param,xLS);
        else
          x0 = xLS;
        end
        [x_reg,~, ~, history] = spectral_ADMM_adaptive(gc, yc,0, a,2,1, ALG_PARAMETER,x0);
%         [x_reg, ~,~, history] = spectral_ADMM_C(gc, yc, a, ALG_PARAMETER,x0);
        A_reg_tmp = devect(full(x_reg),n,p,K); % convert to (n,n,p,K) format
        A_reg(:,:,:,:,ii) = A_reg_tmp; % this is for arranging result into parfor format
        [x_cls,ls_flag(1,ii)] = constrained_LS_D(gc,yc,find(x_reg));
        A_cls =devect(full(x_cls),n,p,K);
        A(:,:,:,:,ii) = A_cls;
        
        tmp_reg = reshape(x_reg,[p*K,n^2]);
        tmp_ls = reshape(x_cls,[p*K,n^2]);
        tmp_ls =sqrt(sum(tmp_ls.^2,1));
        tmp_reg = sqrt(sum(tmp_reg.^2,1));
        df = length(find(tmp_reg))+(p*K-1)*sum(tmp_reg./tmp_ls,'omitnan');
        score(1,ii) = model_selection_S(Y,H,A_cls,df,'LLH_full');
        tmp_ind = cell(1,K);
        diag_ind=1:n+1:n^2;
        for kk=1:K
          tmp_ind{kk} = setdiff(find(squeeze(A_reg_tmp(:,:,1,kk))),diag_ind);
%           ind_VAR{kk} = find(A);
        end
        ind(1,ii) = {tmp_ind};
        [ind_common(1,ii),ind_differential(1,ii)] = split_common_diff(tmp_ind,[n,p,K]); % find common and differential off-diagonal nonzero index
        flag(1,ii) = history.flag;
        if flag(1,ii) ==-1
            fprintf('max iteration exceed at grid (%d)\n',ii)
        end
end
% GIC_LIST = {'bic_lasso','bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
GIC_LIST = fieldnames(score);
for nn=1:length(GIC_LIST)
[~,M.index.(GIC_LIST{nn})] = min([score.(GIC_LIST{nn})]);
end

for ii=1:GridSize
    M.model(ii).stat.model_selection_score = score(ii);
    M.model(ii).A_reg = A_reg(:,:,1:p,:,ii);
    M.model(ii).A = A(:,:,1:p,:,ii);
    M.model(ii).GC = squeeze(sqrt(sum(M.model(ii).A.^2,3)));
    for kk=1:K
        M.model(ii).ind_VAR{kk} = find(M.model(ii).A(:,:,:,kk));
    end
    M.model(ii).ind = ind(ii);
    M.model(ii).ind_common = ind_common(ii);
    M.model(ii).ind_differential = ind_differential(ii);
    M.model(ii).flag = flag(ii);
    M.model(ii).ls_flag = ls_flag(ii);
  
end

M.time = toc(t1);
end
