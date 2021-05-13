function M = cvx_CGN(y,varargin)
%% This core function estimate the common part of Granger network
%  using formulation cvx-CGN
%  input: y, multiple multivariate time-series with dimension (n,T,K), n is
%  #channels, T is time-points, K is # of models to be estimated.
%      : p, VAR order (default, p=1)
%      : GridSize, the resolution of solution path (default, GridSize=30)
%      : weight_def, choice of weighting, weight_def = 'static': no weight
%                                            = 'adaptive_L' (default) : included weight, available when T>=np
%
% parameter.p : 1,2,3, (default = 1)
% parameter.gridsize:  integer (default = 20x20)
% parameter.penalty_weight: 'uniform' (uniform weight = 1)  'LS'  (reciprocal of LS estimate) (default = 'LS')
% parameter.noisecov: 'full', 'diag', 'identity'   (default = 'full')  This is the structure of covariance of residual error
% parameter.qnorm = 'cvx','ncvx'  (1,1/2) % toggle

% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce
[n,T,K] = size(y);
len_varargin = length(varargin);
if isempty(varargin)
    p=1;
    GridSize = 30;
    weight_def = 'adaptive_L';
    toggle = 'LLH_full';
elseif len_varargin==1
    p = varargin{1};
    GridSize = 30;
    weight_def = 'adaptive_L';
    toggle = 'LLH_full';
elseif len_varargin ==2
    p = varargin{1};
    GridSize = varargin{2};
    weight_def = 'adaptive_L';
    toggle = 'LLH_full';
elseif len_varargin ==3
    p = varargin{1};
    GridSize = varargin{2};
    weight_def = varargin{3};
    toggle = 'LLH_full';
elseif len_varargin ==4
    p = varargin{1};
    GridSize = varargin{2};
    weight_def = varargin{3};
    toggle = varargin{4};
else
    error('must be atmost 5 input')
end
eff_T = T-p;
H = zeros(n*p,eff_T,K);
Y = zeros(n,eff_T,K);
disp('Generating H matrix')
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc] = vectorize_VAR(Y,H,[n,p,K,eff_T]);

% produce L1 L2 from [n,p,K] CGN
switch parameter.formulation

    L1 = XXXX
    L2 = XXXXX
    parameter.weight = vector of two/one column (depend on formulation)
    
    
    xLS = gc\yc;

% use this one
lambda_range = grid_generattion (gc,yc,parameter.gridsize,parameter.weight, parameter.penalty_weight,L1,L2)


% create lambda range 
[Lambda_1,Lambda_2,opt] = grid_generation(gc,yc,GridSize,ALG_PARAMETER,qq,weight_def);

% Algorithm parameter
if isnull(ALG_PARAMETER)
    set default value of ALG_PARAMETER.XXX = good default
    

if eff_T>n*p
    init_cvx = 0;
else
    init_cvx = 1;
end
ALG_PARAMETER.PRINT_RESULT=0;
ALG_PARAMETER.IS_ADAPTIVE =1;
ALG_PARAMETER.dim = [n,p,K,p,p*K];
ALG_PARAMETER.rho_init = 1;
ALG_PARAMETER.epscor = 0.5;
ALG_PARAMETER.Ts = 2;
ALG_PARAMETER.is_chol = 1;
ALG_PARAMETER.multiplier = 2;
ALG_PARAMETER.toggle = 'formulationD';
ALG_PARAMETER.gamma = 1; % for adaptive case
ALG_PARAMETER.is_spectral = 1;
disp('calculating Lambda max')
qq=1; %convex case


ALG_PARAMETER.L1 = opt.L1; % remove
ALG_PARAMETER.L2 = opt.L2; % remove
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

% use this one    
adaptive_ADMM(gc,yc,L1,L2,a1,a2,ALG_PARAMETER) % call prox of q = 1/2, heuristic update rho
spectral_ADMM(gc,yc,L1,L2,a1,a2,ALG_PARAMETER) % call prox of q = 1, spectral ADMM
    
CHOOSE ALGORITHM
if non-convex
    run function adaptive_ADMM
    else 
    run function spectral_ADMM

for ii=1:GridSize
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
    score(1,ii) = model_selection_S(Y,H,A_cls,df,toggle);
    tmp_ind = cell(1,K);
    diag_ind=1:n+1:n^2;
    for kk=1:K
        tmp_ind{kk} = setdiff(find(squeeze(A_reg_tmp(:,:,1,kk))),diag_ind);
    end
    ind(1,ii) = {tmp_ind};
    [ind_common(1,ii),ind_differential(1,ii)] = split_common_diff(tmp_ind,[n,p,K]); % find common and differential off-diagonal nonzero index
    flag(1,ii) = history.flag;
    if flag(1,ii) ==-1
        fprintf('max iteration exceed at grid (%d)\n',ii)
    end
end
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
M.flag = reshape([M.model.flag],[1,GridSize]);
M.LLH_type = toggle;
M.time = toc(t1);
end
