function [x,L1x,L2x, history] = adaptive_ADMM(G, b,a1,a2,L1,L2,ALG_PARAMETER)
% This function solve the problem 
% min_x (1/2)||Gx-b||_2^2 + a1||L_1x||_{2,q} + a2||L_2x||_{2,q}
% by ADMM with adaptive penalty parameter \rho
% The update rule is stated in manuscript: https://arxiv.org/abs/2105.07196
% Input  G: vectorized H matrix for all k
%        b: vectorized Y matrix for all k
%       a1: control differential density
%       a2: control common density
%       L1: desired affine transformation
%       L2: desired affine transformation
%       ALG_PARAMETER: generated from gen_alg_params.m
% Output 
% x: converged solution
% L1x: L1*x
% L2x: L2*x
% history: history of algorithm
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce

pp=2;  % L2 norm for the regularization
qq=0.5;  % L1/2 for regularization
FREQ_PRINT = ALG_PARAMETER.FREQ_PRINT;  % print to console every FREQ_PRINT iteration
MAXITERS = ALG_PARAMETER.MAXITERS;
ABSTOL = ALG_PARAMETER.ABSTOL; %1e-7
RELTOL = ALG_PARAMETER.RELTOL; %1e-5
PRINT_RESULT = ALG_PARAMETER.PRINT_RESULT;
IS_ADAPTIVE = ALG_PARAMETER.IS_ADAPTIVE;
n = ALG_PARAMETER.dim(1);
p = ALG_PARAMETER.dim(2);
K = ALG_PARAMETER.dim(3);
m1 = ALG_PARAMETER.dim(4);
m2 = ALG_PARAMETER.dim(5);
L1 = sparse(L1);
L2 = sparse(L2);
is_chol = ALG_PARAMETER.is_chol;
rho = ALG_PARAMETER.rho_init;
Ts = ALG_PARAMETER.Ts;
multiplier = ALG_PARAMETER.multiplier;
x0 = ALG_PARAMETER.x0;
t_start = tic;



% store variables

nn = size(L1,2); % dimension of primal variable x
nd = size(L2,1); % Dimension of Px + Dimension of Dx
np = size(L1,1);
Gtb = G'*b;

A = [L1 ; L2];
At= A';
L2t = sparse(L2');
L1t = sparse(L1');
L2tL2 = sparse(L2'*L2);
GtG = sparse(G'*G);
L1tL1 = sparse(L1'*L1);
L1tL1pL2tL2 = L2tL2+L1tL1;
% ADMM solver

if PRINT_RESULT
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n','iter', ...
        'r norm','eps pri','s norm','eps dual', 'objective');
end

% initial start
x=x0;

y1 = zeros(np,1); y2 = zeros(nd,1);
L1x = L1*x;
L2x = L2*x;
z1 = L1x;
z2 = L2x;
z1old = z1;
z2old = z2;

objx0 = 0.5*norm(G*x-b)^2 + normpq_adaptive(z1,pp,qq,m1,a1)+normpq_adaptive(z2,pp,qq,m2,a2);
if is_chol
    L = chol(GtG+rho*(L1tL1pL2tL2),'lower');
    L = sparse(L); U = L';
end
t_iteration = tic;
for k=1:MAXITERS
    
    % x1-update
    %     c = Gtb+ P'*(rhotilde*(x2+x3)-(z1+z2));
    c = Gtb+ L1t*(rho*z1+y1) + L2t*(rho*z2+y2);
    if is_chol
        
        x = U \ (L \ c);
    else
        x = (GtG+rho*(L1tL1pL2tL2))\c;
    end
    % x2-update
    
    L1x = L1*x;
    z1 = prox_pq_eff_adaptive(L1x-y1/rho,a1/rho,[pp,qq,m1]);
    
    % x3-update
    L2x = L2*x;
    z2 = prox_pq_eff_adaptive(L2x-y2/rho,a2/rho,[pp,qq,m2]);
    y1 = y1 + rho*(-L1x+z1);
    y2 = y2 + rho*(-L2x+z2);
    if  (k>3) && (mod(k,Ts)==0)&&(IS_ADAPTIVE)
        if (history.r_norm(k-1) < history.eps_pri(k-1))
            IS_ADAPTIVE = 0;
            rho_change = 1;
            rho = rho*multiplier;
            %             disp('----------------TURN ADAPTIVE OFF----------------')
        else
            rho_change = 1;
            rho = rho*multiplier;
        end
        if (rho_change)&&(is_chol)
            L = chol(GtG+rho*(L1tL1pL2tL2),'lower');
            L = sparse(L); U = L';
        end
        
    end
    
    
    
    
    
    % stopping criterion
    
    history.t(k) = toc(t_iteration);
    obj = 0.5*norm(G*x-b)^2 + normpq_adaptive(L1x,pp,qq,m1,a1)+normpq_adaptive(L2x,pp,qq,m2,a2);
    history.nz_count(k) = length(find(z1))+length(find(z2));
    history.sp_difference(k) = length(setdiff(find(z1),find(z1old)))+length(setdiff(find(z1old),find(z1)))+length(setdiff(find(z2),find(z2old)))+length(setdiff(find(z2old),find(z2)));
    history.objval(k) = obj;
    history.r_norm(k) = norm([L1x-z1;L2x-z2]);
    history.s_norm(k)  = norm(rho*At*([z1 - z1old;z2-z2old]));
    if k==1
        history.reldiff_norm(k)  = 1;
    else
        history.reldiff_norm(k)  = norm(x-xold)/norm(xold);
    end
    
    history.eps_pri(k) = sqrt(np+nd)*ABSTOL + RELTOL*max(norm(A*x), norm([z1;z2]));
    history.eps_dual(k)= sqrt(nn)*ABSTOL + RELTOL*norm(At*[y1;y2]);
    history.Lagrange(k) = obj + rho/2*norm([L1x-z1;L2x-z2]+[y1/rho;y2/rho])^2;
    history.rho(k) = rho;
    
    if (PRINT_RESULT && mod(k,FREQ_PRINT) == 0)
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n',k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), ...
            history.objval(k));
    end
    
    if ((history.r_norm(k) < history.eps_pri(k) && ...
            history.s_norm(k) < history.eps_dual(k) ))
        history.fit = 0.5*norm(G*x-b)^2;
        history.flag = 0;
        break;
    end
    z1old = z1;
    z2old = z2;
    xold = x;
end
if (k==MAXITERS)
    history.flag = -1;
end

% if PRINT_RESULT
t_end = toc(t_start);
history.total_time = t_end;
history.tpi = t_end/k;
% end
% return sparse result
if length(z1)==length(z2) % for DGN, CGN
    z1((z2==0)) = 0;
end

if all(z1==0,'all')
    history.zero_flag = 1;
else
    history.zero_flag = 0;
end

tmp = z1;
tmp_2 = double(L1~=0)*x;
tmp(tmp~=0) = tmp_2(tmp~=0);
X = reshape(x,p,K,n^2); % x1 is not sparse
X2 = reshape(tmp,p,K,n^2-n); % tmp is offdiagonal entries only and sparse
IND_offdiag = setdiff((1:n^2)',(1:n+1:n^2)','rows');
X(:,:,IND_offdiag) = X2;
x = X(:);
L2x = z2;
L1x = z1;
history.objval = [objx0 history.objval];
end

