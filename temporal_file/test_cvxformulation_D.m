function M = test_cvxformulation_D(y,varargin)
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
% global x_cheat
[n,T,K] = size(y);
len_varargin = length(varargin);
% toggle = 'static';
% toggle = 'adaptive_P';
toggle = 'adaptive_L';
if isempty(varargin)
  p=1;
  GridSize = 30;
elseif len_varargin==1
  p = varargin{1};
  GridSize = 30;
elseif len_varargin ==2
  p = varargin{1};
  GridSize = varargin{2};
else
  error('must be atmost 3 input')
end
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
ALG_PARAMETER.Ts = 2;
ALG_PARAMETER.is_chol = 1;
ALG_PARAMETER.multiplier = 2;
ALG_PARAMETER.toggle = 'formulationD';
ALG_PARAMETER.gamma = 1; % for adaptive case

disp('calculating Lambda max')
qq=1; %convex case
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
for ii=1:GridSize
    a1 = Lambda_1(:,ii);
    A_reg = zeros(n,n,p,K,1,GridSize);
    A = zeros(n,n,p,K,1,GridSize);
    ls_flag = zeros(1,GridSize);
    ind_common = cell(1,GridSize);
    ind_differential = cell(1,GridSize);
    flag = zeros(1,GridSize);
    ind = cell(1,GridSize);

   	parfor jj=1:GridSize
        fprintf('Grid : (%d,%d)/(%d, %d) \n',ii,jj,GridSize,GridSize)
        if init_cvx
            cvx_param = ALG_PARAMETER;
            cvx_param.Ts = 2;
          [x0, ~, ~] = spectral_ADMM_adaptive(gc, yc, a1, Lambda_2(:,jj),2,1, cvx_param);
        else
          x0 = xLS;
        end
        [x_reg, ~,~, history] = spectral_ADMM_adaptive(gc, yc, a1, Lambda_2(:,jj),2,1, ALG_PARAMETER,x0);
        A_reg_tmp = devect(full(x_reg),n,p,K); % convert to (n,n,p,K) format
        A_reg(:,:,:,:,1,jj) = A_reg_tmp; % this is for arranging result into parfor format
        [x_cls,ls_flag(1,jj)] = constrained_LS_D(gc,yc,find(x_reg));
        A_cls =devect(full(x_cls),n,p,K);
        A(:,:,:,:,1,jj) = A_cls;
        score(1,jj) = model_selection(Y,A_cls);
        tmp_ind = cell(1,K);
        diag_ind=1:n+1:n^2;
        for kk=1:K
          tmp_ind{kk} = setdiff(find(squeeze(A_reg_tmp(:,:,1,kk))),diag_ind);
        end
        ind(1,jj) = {tmp_ind};
        [ind_common(1,jj),ind_differential(1,jj)] = split_common_diff(tmp_ind,[n,p,K]); % find common and differential off-diagonal nonzero index
        flag(1,jj) = history.flag;
        if flag(1,jj) ==-1
            fprintf('max iteration exceed at grid (%d,%d)\n',ii,jj)
        end
    end
    tmp_struct.stat.model_selection_score(ii,:) = score;
    tmp_struct.A_reg(:,:,1:p,1:K,ii,:) = A_reg;
    tmp_struct.A(:,:,1:p,1:K,ii,:) = A;
    tmp_struct.ind(ii,:) = ind;
    tmp_struct.ind_common(ii,:) = ind_common;
    tmp_struct.ind_differential(ii,:) = ind_differential;
    tmp_struct.flag(ii,:) = flag;
    tmp_struct.ls_flag(ii,:) = ls_flag;
end
% GIC_LIST = {'bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
GIC_LIST = fieldnames(tmp_struct.stat.model_selection_score);
for nn=1:length(GIC_LIST)
    if strcmp(GIC_LIST{nn},'L')&&strcmp(GIC_LIST{nn},'df')&&strcmp(GIC_LIST{nn},'SSE')
        continue;
    end
[~,M.index.(GIC_LIST{nn})] = min([tmp_struct.stat.model_selection_score.(GIC_LIST{nn})]);
end
for ii=1:GridSize
  for jj=1:GridSize
    M.model(ii,jj).stat.model_selection_score = tmp_struct.stat.model_selection_score(ii,jj);
    M.model(ii,jj).A_reg = tmp_struct.A_reg(:,:,1:p,1:K,ii,jj);
    M.model(ii,jj).A = tmp_struct.A(:,:,1:p,1:K,ii,jj);
    M.model(ii,jj).GC = squeeze(sqrt(sum(M.model(ii,jj).A.^2,3)));
    for kk=1:K
        M.model(ii,jj).ind_VAR{kk} = find(M.model(ii,jj).A(:,:,:,kk));
    end
    M.model(ii,jj).ind = tmp_struct.ind(ii,jj);
    M.model(ii,jj).ind_common = tmp_struct.ind_common(ii,jj);
    M.model(ii,jj).ind_differential = tmp_struct.ind_differential(ii,jj);
    M.model(ii,jj).flag = tmp_struct.flag(ii,jj);
  end
end

M.time = toc(t1);
end
