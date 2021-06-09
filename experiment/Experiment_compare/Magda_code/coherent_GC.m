function M = coherent_GC(y,P,varargin)
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

P1 = P;
% ind_diag = setdiff((1:n^2*p*K)',P1*(1:n^2*p*K)');
P2 = eye(n^2*p*K);
P2(P1*(1:n^2*p*K)',:) = [];

Lambda = logspace(-3,0,GridSize);
H = zeros(n*p,T-p,K);
Y = zeros(n,T-p,K);
disp('Generating H matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc] = vectorize_VAR(Y,H,[n,p,K,T-p]);
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
ALG_PARAMETER.PRINT_RESULT=1;
ALG_PARAMETER.IS_ADAPTIVE =1;
ALG_PARAMETER.L1 = P1;
ALG_PARAMETER.L2 = P2;
ALG_PARAMETER.dim = [n,p,K,p,p*K];
ALG_PARAMETER.rho_init = 1;
ALG_PARAMETER.epscor = 0.1;
ALG_PARAMETER.Ts = 2;
t1 = tic;


    A_reg = zeros(n,n,p,K,GridSize);
    A = zeros(n,n,p,K,GridSize);
    ls_flag = zeros(1,GridSize);
    ind_common = cell(1,GridSize);
    ind_differential = cell(1,GridSize);
    flag = zeros(1,GridSize);
    ind = cell(1,GridSize);
for ii=1:GridSize
    a = Lambda(ii);
        fprintf('Grid : (%d)/(%d) \n',ii,GridSize)
        if init_cvx
            cvx_param = ALG_PARAMETER;
            cvx_param.Ts = 2;
          [x0, ~, ~] = spectral_ADMM(gc, yc, P,P,0, Lambda(ii),2,1, cvx_param,xLS);
        else
          x0 = xLS;
        end
        [x_reg, ~,~, history] = spectral_ADMM_C(gc, yc, a, ALG_PARAMETER,x0);
        A_reg_tmp = devect(full(x_reg),n,p,K); % convert to (n,n,p,K) format
        A_reg(:,:,:,:,ii) = A_reg_tmp; % this is for arranging result into parfor format
        [x_cls,ls_flag(1,ii)] = jointvargc_constrainedLS_DGN(gc,yc,find(x_reg));
        A_cls =devect(full(x_cls),n,p,K);
        A(:,:,:,:,ii) = A_cls;
        score(1,ii) = model_selection(Y,H,A_cls,'full');
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
[~,M.index.bic] = min([score.bic]);
[~,M.index.aicc] = min([score.aicc]);

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
