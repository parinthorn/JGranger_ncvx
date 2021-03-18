function M = cvx_FGN(y,varargin)
%% This core function estimates differential parts of Granger network 
% with group sparse, group fused lasso from multiple multivariate
% time-series y using formulation cvx-FGN
% input: y, multiple multivariate time-series with dimension (n,T,K), n is
% #channels, T is time-points, K is # of models to be estimated.
%      : p, VAR order (default, p=1)
%      : GridSize, the resolution of solution grid(GridSize x GridSize) (default, GridSize=30)
%      : weight_def, choice of weighting, weight_def = 'static': no weight
%                                                    = 'adaptive_D' (default) : included weight, available when T>=np
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce.
[n,T,K] = size(y);
len_varargin = length(varargin);
if isempty(varargin)
  p=1;
  GridSize = 30;
  weight_def = 'adaptive_D';
  toggle = 'LLH_full';
  data_concat=0;
elseif len_varargin==1
  p = varargin{1};
  GridSize = 30;
  weight_def = 'adaptive_D';
  toggle = 'LLH_full';
  data_concat = 0;
elseif len_varargin ==2
  p = varargin{1};
  GridSize = varargin{2};
  weight_def = 'adaptive_D';
  toggle = 'LLH_full';
  data_concat = 0;
elseif len_varargin ==3
  p = varargin{1};
  GridSize = varargin{2};
  weight_def = varargin{3};
  toggle = 'LLH_full';
  data_concat=0;
elseif len_varargin ==4
  p = varargin{1};
  GridSize = varargin{2};
  weight_def = varargin{3};
  toggle = varargin{4};
  data_concat = 0;
elseif len_varargin ==5
    p = varargin{1};
  GridSize = varargin{2};
  weight_def = varargin{3};
  toggle = varargin{4};
  data_concat = varargin{5};  
else
    error('must be atmost 6 input')
end
if (data_concat)
    Ktmp = K/2;K=2; % divides to 2 groups and set K=2
    eff_T = Ktmp*(T-p);
    y1 = y(:,:,1:Ktmp);
    y2 = y(:,:,(Ktmp+1):end);
    H = zeros(n*p,eff_T,2);
    Y = zeros(n,eff_T,2);
    disp('Concatenating H, Y matrix')
    for kk=1:Ktmp
        [H(:,(T-p)*(kk-1)+1:(T-p)*(kk),1),Y(:,(T-p)*(kk-1)+1:(T-p)*(kk),1)] = H_gen(y1(:,:,kk),p);
        [H(:,(T-p)*(kk-1)+1:(T-p)*(kk),2),Y(:,(T-p)*(kk-1)+1:(T-p)*(kk),2)] = H_gen(y2(:,:,kk),p);
    end
    disp('vectorizing model')
    [yc,gc] = vectorize_VAR(Y,H,[n,p,2,eff_T]);
else
    eff_T = T-p;
    H = zeros(n*p,eff_T,K);
    Y = zeros(n,eff_T,K);
    disp('Generating H, Y matrix')
    for kk=1:K
        [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
    end
    disp('vectorizing model')
    [yc,gc] = vectorize_VAR(Y,H,[n,p,K,eff_T]);
end
xLS = gc\yc;
if eff_T>n*p
  init_cvx = 0;
else
  init_cvx = 1;
end
ALG_PARAMETER.PRINT_RESULT=0;
ALG_PARAMETER.IS_ADAPTIVE =1;
ALG_PARAMETER.dim = [n,p,K,p,p];
ALG_PARAMETER.rho_init = 1;
ALG_PARAMETER.epscor = 0.5;
ALG_PARAMETER.Ts = 2;
ALG_PARAMETER.is_chol = 1;
ALG_PARAMETER.multiplier = 2;
ALG_PARAMETER.toggle = 'formulationS';
ALG_PARAMETER.gamma = 1; % for adaptive case
ALG_PARAMETER.is_spectral = 1;

disp('setting up solution path')
qq=1; %convex case
[Lambda_1,Lambda_2,opt] = grid_generation(gc,yc,GridSize,ALG_PARAMETER,qq,weight_def);
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

P = opt.L1;
D = opt.L2;
Ind = (1:1:(size(D,2)))';
Dplus=D;Dminus=D;
Dplus(D==-1) = 0;
Dminus(D==1) = 0;
Dminus = abs(Dminus);
Indplus_in = Dplus*Ind;
Indminus_in = abs(Dminus*Ind);

t1 = tic;
for ii=1:GridSize % test 20
    a1 = Lambda_1(:,ii);
    A_reg = zeros(n,n,p,K,1,GridSize);
    A = zeros(n,n,p,K,1,GridSize);
    ls_flag = zeros(1,GridSize);
    ind_common = cell(1,GridSize);
    ind_differential = cell(1,GridSize);
    flag = zeros(1,GridSize);
    ind = cell(1,GridSize);
    parfor (jj=1:GridSize)
        Indplus = Indplus_in;
        Indminus = Indminus_in;
        fprintf('Grid : (%d,%d)/(%d, %d) \n',ii,jj,GridSize,GridSize)
        if init_cvx
            cvx_param = ALG_PARAMETER;
            cvx_param.Ts = 2;
          [x0, ~, ~,~] = spectral_ADMM_adaptive(gc, yc, a1, Lambda_2(:,jj),2,1, cvx_param);
        else
          x0 = xLS;
        end
        [x_reg, Px,Dx, history] = spectral_ADMM_adaptive(gc, yc, a1, Lambda_2(:,jj),2,qq, ALG_PARAMETER,x0);
        A_reg_tmp = devect(full(x_reg),n,p,K); % convert to (n,n,p,K) format
        A_reg(:,:,:,:,1,jj) = A_reg_tmp; % this is for arranging result into parfor format
        x_cls = constrained_LS_S(gc,yc,D,Dx,P,Px,'off');
        A_cls =devect(full(x_cls),n,p,K);
        fused_index=intersect( ...
            union(unique(Indplus(Dplus*x_cls~=0)), ...
                  unique(Indminus(Dminus*x_cls~=0))), ...
            union(Indplus(D*x_cls==0), ...
                  Indminus(D*x_cls==0))); % intuitively, this operation is to find indices of nonzero variables but with zero differences
        tmp =length(find(diff(x_cls(fused_index))==0));
        df = length(find(x_cls))-tmp;
        A(:,:,:,:,1,jj) = A_cls;
        score(1,jj) = model_selection_S(Y,H,A_cls,df,toggle);
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
% GIC_LIST = {'bic_lasso','bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
GIC_LIST = fieldnames(tmp_struct.stat.model_selection_score);
for nn=1:length(GIC_LIST)
[~,M.index.(GIC_LIST{nn})] = min([tmp_struct.stat.model_selection_score.(GIC_LIST{nn})]);
end
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
M.flag = reshape([M.model.flag],[GridSize,GridSize]);
M.LLH_type = toggle;
M.time = toc(t1);
end