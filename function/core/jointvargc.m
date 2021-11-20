function M = jointvargc(y,parameter, ALG_PARAMETER)
% This function computes the formulation CGN, DGN, FGN on the give data y.
% Input
% y: size(y)=[n, T, K] -> [dimension, timepoints, number of groups]
% parameter.p_var: VAR order
% parameter.GridSize: Resolution of the regularization grid search
% parameter.noisecov: the constraint on noise covariance when computing 
% Log-likelihood to be dense matrix, diagonal matrix or multiple of 
% identity matrix ['full', 'diag', 'identity']
% parameter.penalty_weight: ['LS', 'uniform'] adaptive weight or equal weight
% parameter.formulation: ['cgn', 'dgn', 'fgn']
% parameter.qnorm: ['cvx', 'ncvx']
% Output
% M.model: the resulting models in each pair of regularization
% M.flag: divergence flag, 0 is converge, -1 is diverge.
% M.time: inference time (seconds)
% M.GridSize: parameter.GridSize
% M.noisecov: parameter.noisecov
% M.penalty_weight: parameter.penalty_weight
% M.formulation: parameter.formulation
% M.qnorm: parameter.qnorm
%
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce
[n,T,K] = size(y);
p_var = parameter.varorder;
GridSize = parameter.GridSize;
data_concat = parameter.data_concat;
noisecov = parameter.noisecov; %'full', 'diag', 'identity'
penalty_weight = parameter.penalty_weight; % 'LS', 'uniform'
formulation = parameter.formulation; % cgn, dgn, fgn
qnorm = parameter.qnorm; % cvx, ncvx
if isfield(parameter, 'is_YH')
    is_YH = parameter.is_YH;
else
    is_YH = 0;
end
%======= vectorization
if (data_concat)
    Ktmp = K/2;K=2; % divides to 2 groups and set K=2
    eff_T = Ktmp*(T-p_var);
    y1 = y(:,:,1:Ktmp);
    y2 = y(:,:,(Ktmp+1):end);
    H = zeros(n*p_var,eff_T,2);
    Y = zeros(n,eff_T,2);
    disp('Concatenating H, Y matrix')
    for kk=1:Ktmp
        [H(:,(T-p_var)*(kk-1)+1:(T-p_var)*(kk),1),Y(:,(T-p_var)*(kk-1)+1:(T-p_var)*(kk),1)] = H_gen(y1(:,:,kk),p_var);
        [H(:,(T-p_var)*(kk-1)+1:(T-p_var)*(kk),2),Y(:,(T-p_var)*(kk-1)+1:(T-p_var)*(kk),2)] = H_gen(y2(:,:,kk),p_var);
    end
    disp('vectorizing model')
    [b,G] = vectorize_VAR(Y,H,[n,p_var,2,eff_T]);
elseif (is_YH)
    Y = parameter.Y;
    H = parameter.H;
    eff_T = parameter.eff_T;
    disp('vectorizing model')
    [b,G] = vectorize_VAR(Y,H,[n,p_var,K,eff_T]);
else
    eff_T = T-p_var;
    H = zeros(n*p_var,eff_T,K);
    Y = zeros(n,eff_T,K);
    disp('Generating H, Y matrix')
    for kk=1:K
        [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p_var);
    end
    disp('vectorizing model')
    [b,G] = vectorize_VAR(Y,H,[n,p_var,K,eff_T]);
end
xLS = G\b;
[P,~] = offdiagJSS(n,p_var,K);
if ~strcmp(formulation, 'fgn')
    L1 = sparse(P);
    L2 = sparse(P);
else
    Dtmp = diffmat(n,p_var,K);
    D = sparse(Dtmp*P);
    L1 = sparse(P);
    L2 = sparse(D);
end
% Algorithm parameters initialization
if ~exist('ALG_PARAMETER','var')
    ALG_PARAMETER = gen_alg_params(qnorm, formulation);
end
switch formulation
    case 'cgn'
        ALG_PARAMETER.dim = [n,p_var,K,p_var,p_var*K];
    case 'dgn'
        ALG_PARAMETER.dim = [n,p_var,K,p_var,p_var*K];
    case 'fgn'
        ALG_PARAMETER.dim = [n,p_var,K,p_var,p_var];
end
% Lambda bound computation
[Lambdacrit_1,Lambdacrit_2] = gen_critical_lambdas(G,b, xLS,ALG_PARAMETER.dim,qnorm,penalty_weight,formulation);
Lambda = logspace(-6,0,GridSize);
ALG_PARAMETER.x0 = xLS;


% algorithm selection
switch qnorm
    case 'cvx'
        fitmodel = @(a1,a2) spectral_ADMM(G, b,a1,a2,L1,L2,ALG_PARAMETER);
    case 'ncvx'
        fitmodel = @(a1,a2) adaptive_ADMM(G, b,a1,a2,L1,L2,ALG_PARAMETER);
end



switch formulation
    case 'cgn'
        A_reg = zeros(n,n,p_var,K,GridSize);
        A = zeros(n,n,p_var,K,GridSize);
        ls_flag = zeros(1,GridSize);
        ind_common = cell(1,GridSize);
        ind_differential = cell(1,GridSize);
        flag = zeros(1,GridSize);
        ind = cell(1,GridSize);
        total_time = zeros(1,GridSize);
        t1 = tic;
        parfor (ii=1:GridSize, 8) % this can be changed to parfor loops
            fprintf('Grid : (%d)/(%d) \n',ii,GridSize)
            a2 = Lambdacrit_2*Lambda(ii);
            [x_reg,~,~,history] = fitmodel(0,a2);
            total_time(ii) = history.total_time;
            A_reg_tmp = devect(full(x_reg),n,p_var,K); % convert to (n,n,p,K) format
            A_reg(:,:,:,:,ii) = A_reg_tmp; % this is for arranging result into parfor format
            [x_cls,ls_flag(1,ii)] = jointvargc_constrainedLS_DGN(G,b,find(x_reg));
            A_cls =devect(full(x_cls),n,p_var,K);
            A(:,:,:,:,ii) = A_cls;
            tmp_reg = reshape(x_reg,[p_var*K,n^2]);
            tmp_ls = reshape(x_cls,[p_var*K,n^2]);
            tmp_ls =sqrt(sum(tmp_ls.^2,1));
            tmp_reg = sqrt(sum(tmp_reg.^2,1));
            df = length(find(tmp_reg))+(p_var*K-1)*sum(tmp_reg./tmp_ls,'omitnan');
            score(1,ii) = model_selection_S(Y,H,A_cls,df,noisecov);
            tmp_ind = cell(1,K);
            diag_ind=1:n+1:n^2;
            for kk=1:K
                tmp_ind{kk} = setdiff(find(squeeze(A_reg_tmp(:,:,1,kk))),diag_ind);
            end
            ind(1,ii) = {tmp_ind};
            [ind_common(1,ii),ind_differential(1,ii)] = split_common_diff(tmp_ind,[n,p_var,K]); % find common and differential off-diagonal nonzero index
            flag(1,ii) = history.flag;
            if flag(1,ii) ==-1
                fprintf('max iteration exceed at grid (%d)\n',ii)
            end
        end
        GIC_LIST = fieldnames(score);
        for nn=1:length(GIC_LIST)
            if strcmp(GIC_LIST{nn},'LLH_full')||strcmp(GIC_LIST{nn},'LLH_diag')||strcmp(GIC_LIST{nn},'LLH_identity')||strcmp(GIC_LIST{nn},'df')||strcmp(GIC_LIST{nn},'SSE')||strcmp(GIC_LIST{nn},'flag')
                continue;
            end
            [~,M.index.(GIC_LIST{nn})] = min([score.(GIC_LIST{nn})]);
        end
        for ii=1:GridSize
            M.model(ii).stat.model_selection_score = score(ii);
            M.model(ii).A_reg = A_reg(:,:,1:p_var,:,ii);
            M.model(ii).A = A(:,:,1:p_var,:,ii);
            M.model(ii).GC = squeeze(sqrt(sum(M.model(ii).A.^2,3)));
            for kk=1:K
                M.model(ii).ind_VAR{kk} = find(M.model(ii).A(:,:,:,kk));
            end
            M.model(ii).ind = ind(ii);
            M.model(ii).ind_common = ind_common(ii);
            M.model(ii).ind_differential = ind_differential(ii);
            M.model(ii).flag = flag(ii);
            M.model(ii).ls_flag = ls_flag(ii);
            M.model(ii).total_time = total_time(ii);
            
        end
        M.flag = reshape([M.model.flag],[1,GridSize]);
        M.time = toc(t1);
        
        
    case 'dgn'
        t1 = tic;
        for ii=1:GridSize
            a1 = Lambdacrit_1*Lambda(ii);
            A_reg = zeros(n,n,p_var,K,1,GridSize);
            A = zeros(n,n,p_var,K,1,GridSize);
            ls_flag = zeros(1,GridSize);
            ind_common = cell(1,GridSize);
            ind_differential = cell(1,GridSize);
            flag = zeros(1,GridSize);
            ind = cell(1,GridSize);
            total_time = zeros(1,GridSize);
            parfor (jj=1:GridSize, 12)  % this can be changed to parfor loops
                fprintf('Grid : (%d,%d)/(%d, %d) \n',ii,jj,GridSize,GridSize)
                a2 = Lambdacrit_2*Lambda(jj);
                [x_reg,~,~,history] = fitmodel(a1,a2);
                total_time(1, jj) = history.total_time;
                A_reg_tmp = devect(full(x_reg),n,p_var,K); % convert to (n,n,p,K) format
                A_reg(:,:,:,:,1,jj) = A_reg_tmp; % this is for arranging result into parfor format
                [x_cls,ls_flag(1,jj)] = jointvargc_constrainedLS_DGN(G,b,find(x_reg));
                A_cls =devect(full(x_cls),n,p_var,K);
                A(:,:,:,:,1,jj) = A_cls;
                score(1,jj) = model_selection(Y,H,A_cls,noisecov);
                tmp_ind = cell(1,K);
                diag_ind=1:n+1:n^2;
                for kk=1:K
                    tmp_ind{kk} = setdiff(find(squeeze(A_reg_tmp(:,:,1,kk))),diag_ind);
                end
                ind(1,jj) = {tmp_ind};
                [ind_common(1,jj),ind_differential(1,jj)] = split_common_diff(tmp_ind,[n,p_var,K]); % find common and differential off-diagonal nonzero index
                flag(1,jj) = history.flag;
                if flag(1,jj) ==-1
                    fprintf('max iteration exceed at grid (%d,%d)\n',ii,jj)
                end
            end
            tmp_struct.stat.model_selection_score(ii,:) = score;
            tmp_struct.A_reg(:,:,1:p_var,1:K,ii,:) = A_reg;
            tmp_struct.A(:,:,1:p_var,1:K,ii,:) = A;
            tmp_struct.ind(ii,:) = ind;
            tmp_struct.ind_common(ii,:) = ind_common;
            tmp_struct.ind_differential(ii,:) = ind_differential;
            tmp_struct.flag(ii,:) = flag;
            tmp_struct.ls_flag(ii,:) = ls_flag;
            tmp_struct.total_time(ii, :) = total_time;
            
        end
        GIC_LIST = fieldnames(tmp_struct.stat.model_selection_score);
        for nn=1:length(GIC_LIST)
            if strcmp(GIC_LIST{nn},'LLH_full')||strcmp(GIC_LIST{nn},'LLH_diag')||strcmp(GIC_LIST{nn},'LLH_identity')||strcmp(GIC_LIST{nn},'df')||strcmp(GIC_LIST{nn},'SSE')||strcmp(GIC_LIST{nn},'flag')
                continue;
            end
            [~,M.index.(GIC_LIST{nn})] = min([tmp_struct.stat.model_selection_score.(GIC_LIST{nn})]);
        end
        for ii=1:GridSize
            for jj=1:GridSize  % this can be changed to parfor loops
                M.model(ii,jj).stat.model_selection_score = tmp_struct.stat.model_selection_score(ii,jj);
                M.model(ii,jj).A_reg = tmp_struct.A_reg(:,:,1:p_var,1:K,ii,jj);
                M.model(ii,jj).A = tmp_struct.A(:,:,1:p_var,1:K,ii,jj);
                M.model(ii,jj).GC = squeeze(sqrt(sum(M.model(ii,jj).A.^2,3)));
                for kk=1:K
                    M.model(ii,jj).ind_VAR{kk} = find(M.model(ii,jj).A(:,:,:,kk));
                end
                M.model(ii,jj).ind = tmp_struct.ind(ii,jj);
                M.model(ii,jj).ind_common = tmp_struct.ind_common(ii,jj);
                M.model(ii,jj).ind_differential = tmp_struct.ind_differential(ii,jj);
                M.model(ii,jj).flag = tmp_struct.flag(ii,jj);
                M.model(ii,jj).total_time = tmp_struct.total_time(ii, jj);
            end
        end
        M.flag = reshape([M.model.flag],[GridSize,GridSize]);
        M.time = toc(t1);
        
    case 'fgn'
        t1 = tic;
        Ind = (1:1:(size(D,2)))';
        Dplus=D;Dminus=D;
        Dplus(D==-1) = 0;
        Dminus(D==1) = 0;
        Dminus = abs(Dminus);
        Indplus_in = Dplus*Ind;
        Indminus_in = abs(Dminus*Ind);
        for ii=1:GridSize
            a1 = Lambdacrit_1*Lambda(ii);
            A_reg = zeros(n,n,p_var,K,1,GridSize);
            A = zeros(n,n,p_var,K,1,GridSize);
            ls_flag = zeros(1,GridSize);
            ind_common = cell(1,GridSize);
            ind_differential = cell(1,GridSize);
            flag = zeros(1,GridSize);
            ind = cell(1,GridSize);
            total_time = zeros(1,GridSize);
            for jj=1:GridSize
                fprintf('Grid : (%d,%d)/(%d, %d) \n',ii,jj,GridSize,GridSize)
                Indplus = Indplus_in;
                Indminus = Indminus_in;
                a2 = Lambdacrit_2*Lambda(jj);
                [x_reg,Px,Dx,history] = fitmodel(a1,a2);
                total_time(1, jj) = history.total_time;
                A_reg_tmp = devect(full(x_reg),n,p_var,K); % convert to (n,n,p,K) format
                A_reg(:,:,:,:,1,jj) = A_reg_tmp; % this is for arranging result into parfor format
                x_cls = jointvargc_constrainedLS_FGN(G,b,D,Dx,P,Px,'off');
                A_cls =devect(full(x_cls),n,p_var,K);
                fused_index=intersect( ...
                    union(unique(Indplus(Dplus*x_cls~=0)), ...
                    unique(Indminus(Dminus*x_cls~=0))), ...
                    union(Indplus(D*x_cls==0), ...
                    Indminus(D*x_cls==0))); % intuitively, this operation is to find indices of nonzero variables but with zero differences
                tmp =length(find(diff(x_cls(fused_index))==0));
                df = length(find(x_cls))-tmp;
                A(:,:,:,:,1,jj) = A_cls;
                score(1,jj) = model_selection_S(Y,H,A_cls,df,noisecov);
                tmp_ind = cell(1,K);
                diag_ind=1:n+1:n^2;
                for kk=1:K
                    tmp_ind{kk} = setdiff(find(squeeze(A_reg_tmp(:,:,1,kk))),diag_ind);
                end
                ind(1,jj) = {tmp_ind};
                [ind_common(1,jj),ind_differential(1,jj)] = split_common_diff(tmp_ind,[n,p_var,K]); % find common and differential off-diagonal nonzero index
                flag(1,jj) = history.flag;
                if flag(1,jj) ==-1
                    fprintf('max iteration exceed at grid (%d,%d)\n',ii,jj)
                end
                
            end
            tmp_struct.stat.model_selection_score(ii,:) = score;
            tmp_struct.A_reg(:,:,1:p_var,:,ii,:) = A_reg;
            tmp_struct.A(:,:,1:p_var,:,ii,:) = A;
            tmp_struct.ind(ii,:) = ind;
            tmp_struct.ind_common(ii,:) = ind_common;
            tmp_struct.ind_differential(ii,:) = ind_differential;
            tmp_struct.flag(ii,:) = flag;
            tmp_struct.ls_flag(ii,:) = ls_flag;
            tmp_struct.total_time(ii, :) = total_time;
        end
        GIC_LIST = fieldnames(tmp_struct.stat.model_selection_score);
        for nn=1:length(GIC_LIST)
            if strcmp(GIC_LIST{nn},'LLH_full')||strcmp(GIC_LIST{nn},'LLH_diag')||strcmp(GIC_LIST{nn},'LLH_identity')||strcmp(GIC_LIST{nn},'df')||strcmp(GIC_LIST{nn},'SSE')||strcmp(GIC_LIST{nn},'flag')
                continue;
            end
            [~,M.index.(GIC_LIST{nn})] = min([tmp_struct.stat.model_selection_score.(GIC_LIST{nn})]);
        end
        for ii=1:GridSize
            for jj=1:GridSize
                M.model(ii,jj).stat.model_selection_score = tmp_struct.stat.model_selection_score(ii,jj);
                M.model(ii,jj).A_reg = tmp_struct.A_reg(:,:,1:p_var,:,ii,jj);
                M.model(ii,jj).A = tmp_struct.A(:,:,1:p_var,:,ii,jj);
                M.model(ii,jj).GC = squeeze(sqrt(sum(M.model(ii,jj).A.^2,3)));
                for kk=1:K
                    M.model(ii,jj).ind_VAR{kk} = find(M.model(ii,jj).A(:,:,:,kk));
                end
                M.model(ii,jj).ind = tmp_struct.ind(ii,jj);
                M.model(ii,jj).ind_common = tmp_struct.ind_common(ii,jj);
                M.model(ii,jj).ind_differential = tmp_struct.ind_differential(ii,jj);
                M.model(ii,jj).flag = tmp_struct.flag(ii,jj);
                M.model(ii,jj).total_time = tmp_struct.total_time(ii, jj);
            end
        end
        M.flag = reshape([M.model.flag],[GridSize,GridSize]);
        M.time = toc(t1);
end



M.GridSize = GridSize;
M.noisecov = noisecov;
M.penalty_weight = penalty_weight;
M.formulation = formulation;
M.qnorm = qnorm;
end