function model = gen_multi_VAR(PARAMETER,opts,varargin)
% gen_VAR_param generates multiple sparse VAR parameters [A_1, ..., A_p]
% y_k(t) = A_k_1*y_k(t-1) + ... + A_k_p*y_k(t-p) + u(t)
% Where 'k' is k th model
% Input parameters are
% 'PARAMETER' = [n p k]
% 'opts' : a struct that contains field
%          '.differential_density : differential level on similar structure, all
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
differential_density = opts.differential_density;
model.dim = [n,p,K];
model.common_density = opts.common_density;
model.differential_density = opts.differential_density;
model.type = opts.type;

if numel(varargin) > 0
    DMAT = varargin{1};
end
switch model.type
  case 'common'
    tmp = [1 2;repmat([1 0],K-1,1)];
    model.differential_density = 0;
  case 'similar'
    K=K+1; % if all K-1 models are similar to the first one, the first one will not have any differential part
           % so we generate K+1 models and truncate the first model
    tmp = [1 2;repmat([1 1],K-1,1)];
  case 'differential' % must be generated after 'common' type
    tmp = [1 2;repmat([1 1],K-1,1)];
end

model.A = zeros(n,n,p,K); % allocate ground truth VAR parameter

% our strategy to generate K VAR models that
% 1. C type: all K models share same pattern but can have different value of VAR coefficient
% 2. D type: add differential part to C type
% 3. S type: make the parameter in C type similar
model.ind = cell(K,1);
model.seed = randi(100000000); % seed
rng(model.seed)
for kk=1:K
    fprintf('Group %d\n',kk)
    if exist('DMAT','var')
        [model.A(:,:,:,kk),model.ind{kk},model.VAR_spectrum{kk},model.ind_VAR{kk}] = gen_VAR(n,p,common_density,differential_density,1,DMAT(:,:,:,kk)); % look for generated C model
    else
        [model.A(:,:,:,kk),model.ind{kk},model.VAR_spectrum{kk},model.ind_VAR{kk}] = gen_VAR(n,p,common_density,differential_density,tmp(kk,2),model.A(:,:,:,tmp(kk,1)));

    end
    model.GC(:,:,kk) = sqrt(sum(model.A(:,:,:,kk).^2,3));
end
if strcmp(opts.type,'similar') % we generated first model and the rest is similar to the first one
                               % the first one will not have any
                               % differential part compare to others.
                               % To generalize, we remove first model.
    model.A = model.A(:,:,:,2:end);
    model.GC = model.GC(:,:,2:end);
    model.ind = model.ind(2:end);
    model.ind_VAR = model.ind_VAR(2:end);
    model.VAR_spectrum = model.VAR_spectrum(2:end);
    K=K-1;
end

% generating common index, differential index, removing off diagonal index
diag_ind = 1:n+1:n^2;
model.ind_common = setdiff(1:n^2,diag_ind);
model.ind_differential = cell(K,1);
for kk=1:K
  model.ind_common = intersect(model.ind_common,model.ind{kk});
end
for kk=1:K
  model.ind{kk} = setdiff(model.ind{kk},diag_ind); % remove diagonal parts
  model.ind_differential{kk} = setdiff(model.ind{kk},model.ind_common); % remove common part, so that
                                                  % the remaining is differential part
end
model = orderfields(model,[1 4 2 3 9 6 11 12 8 10 5 7]);
end
