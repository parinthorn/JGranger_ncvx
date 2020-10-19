function M = gen_multi_VAR(PARAMETER,opts,varargin)
% gen_VAR_param generates multiple sparse VAR parameters [A_1, ..., A_p]
% y_k(t) = A_k_1*y_k(t-1) + ... + A_k_p*y_k(t-p) + u(t)
% Where 'k' is k th model
% Input parameters are
% 'PARAMETER' = [n p k]
% 'opts' : a struct that contains field
%          '.diff_density : differential level on similar structure, all
%          similar model have same differential level
%          '.common_density' : nonzero common_density of VAR coefficient
%          '.type' : 'common', 'similar', 'differential' or manually
%          specify Kx2 vector where type(i,:) = [a,b] refers to
%           index (k=i) has relation type 'b' with index (k=j) (j<i)
%                    'b' = 2, 1 ,0 if no relation,  similar edges, common sparsity
%                    respectively
% The output are
% 'M' : M contains k struct with field
%        '.ind_nz' : nonzero indices
%        '.A' : VAR parameter with size [n,n,p,k]
%        '.info' : contains cells ind_nz
%% Function init
n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);
common_density = opts.common_density;
diff_density = opts.diff_density;
M = struct();
M.info.dim = [n,p,K];
if numel(varargin) > 0
    DMAT = varargin{1};
end
switch opts.type
  case 'common'
    tmp = [1 2;repmat([1 0],K-1,1)];
  case 'similar'
    K=K+1; % if all K-1 models are similar to the first one, the first one will not have any differential part
           % so we generate K+1 models and truncate the first model
    tmp = [1 2;repmat([1 1],K-1,1)];
  case 'differential' % must be generated after 'common' type
    tmp = [1 2;repmat([1 1],K-1,1)];
end

M.A = zeros(n,n,p,K); % allocate ground truth VAR parameter
M.info.common_density = opts.common_density;
M.info.type = opts.type;
M.info.diff_density = opts.diff_density;
if opts.diff_density ==0
  opts.diff_density = 0;
end

% our strategy to generate K VAR models that
% 1. C type: all K models share same pattern but can have different value of VAR coefficient
% 2. D type: add differential part to C type
% 3. S type: make the parameter in C type similar
for kk=1:K
    fprintf('Group %d\n',kk)
    if exist('DMAT','var')
        [M.info.ind_nz{kk},M.A(:,:,:,kk),M.info.spectral_radius{kk},M.info.seed(kk)] = gen_VAR(n,p,common_density,diff_density,1,DMAT(:,:,:,kk)); % look for generated C model
    else
        [M.info.ind_nz{kk},M.A(:,:,:,kk)] = gen_VAR(n,p,common_density,diff_density,tmp(kk,2),M.A(:,:,:,tmp(kk,1)));

    end
end
if strcmp(opts.type,'similar')
    M.A = M.A(:,:,:,2:end);
end
end
