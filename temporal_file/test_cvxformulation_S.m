function M = test_cvxformulation_S(y,P,D,varargin)
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
  error('must be atmost 5 input')
end
Lambda = logspace(-6,0,GridSize);
H = zeros(n*p,T-p,K);
Y = zeros(n,T-p,K);
disp('Generating H matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc] = vectorize_VAR(Y,H,[n,p,K,T]);
disp('calculating Lambda max')
Lmax = lambdamax_grouplasso(gc,yc,[n ,p ,K]);
Lambda = Lambda*Lmax;
xLS = gc\yc;
if T>n*p
  init_cvx = 0;
else
  init_cvx = 1;
end
M.lambda_crit = Lmax;
M.lambda_range = [Lambda(1) Lambda(end)];
M.GridSize = GridSize;
M.flag = zeros(GridSize);
ALG_PARAMETER.PRINT_RESULT=0;
ALG_PARAMETER.IS_ADAPTIVE =1;
ALG_PARAMETER.L1 = P;
ALG_PARAMETER.L2 = D;
ALG_PARAMETER.dim = [n,p,K,p,p*K];
ALG_PARAMETER.rho_init = 1;
ALG_PARAMETER.epscor = 0.1;
ALG_PARAMETER.Ts = 100;
ALG_PARAMETER.is_chol = 1;
ALG_PARAMETER.multiplier = 2;
ALG_PARAMETER.toggle = 'formulationS';

Ind = (1:1:(size(D,2)))';
Dplus=D;Dminus=D;
Dplus(D==-1) = 0;
Dminus(D==1) = 0;
Dminus = abs(Dminus);
Indplus_in = Dplus*Ind;
Indminus_in = abs(Dminus*Ind);

t1 = tic;
for ii=1:GridSize % test 12
    a1 = Lambda(ii);
    A_reg = zeros(n,n,p,K,1,GridSize);
    A = zeros(n,n,p,K,1,GridSize);
    ls_flag = zeros(1,GridSize);
    ind_common = cell(1,GridSize);
    ind_differential = cell(1,GridSize);
    flag = zeros(1,GridSize);
    ind = cell(1,GridSize);
    parfor jj=1:GridSize %test 18
        Indplus = Indplus_in;
        Indminus = Indminus_in;
        fprintf('Grid : (%d,%d)/(%d, %d) \n',ii,jj,GridSize,GridSize)
        if init_cvx
            cvx_param = ALG_PARAMETER;
            cvx_param.Ts = 2;
          [x0, ~, ~,~] = spectral_ADMM(gc, yc, a1, Lambda(jj),2,1, cvx_param);
        else
          x0 = xLS;
        end
        [x_reg, Px,Dx, history] = spectral_ADMM(gc, yc, a1, Lambda(jj),2,1, ALG_PARAMETER,x0);
        A_reg_tmp = devect(full(x_reg),n,p,K); % convert to (n,n,p,K) format
        A_reg(:,:,:,:,1,jj) = A_reg_tmp; % this is for arranging result into parfor format
        x_cls = constrained_LS_S(gc,yc,D,Dx,P,Px,'off');
        A_cls =devect(full(x_cls),n,p,K);
        fused_index=intersect( ...
            union(unique(Indplus(Dplus*x_cls~=0)), ...
                  unique(Indminus(Dminus*x_cls~=0))), ...
            union(Indplus(D*x_cls==0), ...
                  Indminus(D*x_cls==0)));
        tmp = (reshape(x_cls(fused_index),[p,length(x_cls(fused_index))/p]));
        df = length(find(x_cls))-length(unique(tmp(1,:)));
        A(:,:,:,:,1,jj) = A_cls;
%         error('please insert model_selection for formulationS')
        score(1,jj) = model_selection_S(Y,A_cls,df);
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
    tmp_struct.A_reg(:,:,1:p,:,ii,:) = A_reg;
    tmp_struct.A(:,:,1:p,:,ii,:) = A;
    tmp_struct.ind(ii,:) = ind;
    tmp_struct.ind_common(ii,:) = ind_common;
    tmp_struct.ind_differential(ii,:) = ind_differential;
    tmp_struct.flag(ii,:) = flag;
    tmp_struct.ls_flag(ii,:) = ls_flag;
end
[~,M.index.bic] = min([tmp_struct.stat.model_selection_score.bic]);
[~,M.index.aicc] = min([tmp_struct.stat.model_selection_score.aicc]);
for ii=1:GridSize
  for jj=1:GridSize
    M.model(ii,jj).stat.model_selection_score = tmp_struct.stat.model_selection_score(ii,jj);
    M.model(ii,jj).A_reg = tmp_struct.A_reg(:,:,1:p,:,ii,jj);
    M.model(ii,jj).A = tmp_struct.A(:,:,1:p,:,ii,jj);
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
