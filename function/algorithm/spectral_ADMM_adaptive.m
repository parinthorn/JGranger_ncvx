function [x,L1x,L2x, history] = spectral_ADMM_adaptive(G, b,a1,a2,pp,qq,PARAMETER,varargin)
% Spectral ADMM for Sparse Granger Causality with non-convex penalties (formulation S)
% This program solve the problem
%  \min_{x} (1/2)||Gx-b||_{2}^{2} + a1 * \mynorm{L1x}{pp}{qq}{m1} + a2 * \mynorm{L2x}{pp}{qq}{m2}
% The input parameters are
%           PARAMETER.dim : [n,p,K,m1,m2]
%                    .rho       : initial ADMM parameter
%                    .epscor    : correlation of spectral step-size and the
%                    difference of 
%                    .Ts        :
%                    .varargin  : initial solution
% The output are
%           x         : solution of the problem in sparse format
%           Px        : A projected solution by matrix P
%           history   : Optimization history which has the data of
%                  .r_norm     : primal residual norm
%                  .s_norm     : dual residual norm
%                  .eps_pri     : primal residual norm convergence threshold
%                  .eps_dual     : dual residual norm convergence threshold
%                  .Lagrange     : value of augmented Lagrangian
%                  .rho     : ADMM penalty in each iteration
%                  .fit     : the sum square loss of final solution
%                  .tpi     : primal residual norm


FREQ_PRINT = 100;
MAXITERS = 50000;
ABSTOL = 1e-7; %-7
RELTOL = 1e-5; %-5
% ABSTOL = 1e-6; %-7
% RELTOL = 1e-4; %-5
eps_rate = 1e-2;
itr_same_index = 0;
SPARSITY_CONVERGE = 0;
PRINT_RESULT = PARAMETER.PRINT_RESULT;
IS_ADAPTIVE = PARAMETER.IS_ADAPTIVE;
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);
m1 = PARAMETER.dim(4);
m2 = PARAMETER.dim(5);
L1 = sparse(PARAMETER.L1);
L2 = sparse(PARAMETER.L2);
is_chol = PARAMETER.is_chol;
multiplier = PARAMETER.multiplier;
toggle = PARAMETER.toggle;
is_spectral = PARAMETER.is_spectral;

rho = PARAMETER.rho_init;
epscor = PARAMETER.epscor;
Ts = PARAMETER.Ts;
MAX_CHECK = max([1000,5*Ts]);
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
optargin = size(varargin,2);
if optargin == 0
    x = zeros(nn,1);
else
    x = varargin{1};
end

y1 = zeros(np,1); y2 = zeros(nd,1);
y1hat = y1;
y2hat = y2;
L1x = L1*x;
L2x = L2*x;
z1 = L1x;
z2 = L2x;
xold = x;
z1old = z1;
z2old = z2;
y1old = y1;
y2old = y2;
objx0 = 0.5*norm(G*x-b)^2 + normpq_adaptive(z1,pp,qq,m1,a1)+normpq_adaptive(z2,pp,qq,m2,a2);
if is_chol
L = chol(GtG+rho*(L1tL1pL2tL2),'lower'); 
L = sparse(L); U = L';
end
y1hatk0 = y1hat;
y2hatk0 = y2hat;
y1k0 = y1;
y2k0 = y2;
z1k0 = z1;
z2k0 = z2;
xk0 = x;

k0 = 2;
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
    if  (k>3) && (mod(k,Ts)==0)&& (IS_ADAPTIVE||(history.r_norm(k-1) > history.eps_pri(k-1))) %&& (history.r_norm(k-1) > history.eps_pri(k-1))
        primal_rate = abs((history.r_norm(k-1)-history.r_norm(k0))/history.r_norm(k0));
        dual_rate = abs((history.s_norm(k-1)-history.s_norm(k0))/history.s_norm(k0));
        y1hat =  y1old + rho*(-L1x+z1old);
        y2hat = y2old + rho*(-L2x+z2old);
        delta_H = A*(x-xk0);
        delta_yhat = [y1hat-y1hatk0;y2hat-y2hatk0];
        delta_y = [y1-y1k0;y2-y2k0];
        delta_G = -([z1-z1k0;z2-z2k0]);
        tmp1 = (delta_H'*delta_yhat);
        akMG = tmp1/(norm(delta_H,2)^2);
        akSD = norm(delta_yhat,2)^2/tmp1;
        if 2*akMG>akSD
            ak=akMG;
        else
            ak = akSD-akMG/2;
        end
        bkSD = norm(delta_y,2)^2/(delta_G'*delta_y);
        bkMG = (delta_G'*delta_y)/norm(delta_G,2)^2;
        if 2*bkMG>bkSD
            bk=bkMG;
        else
            bk = bkSD-bkMG/2;
        end
        akcor = delta_H'*delta_yhat/(norm(delta_H,2)*norm(delta_yhat,2));
        bkcor = delta_G'*delta_y/(norm(delta_G,2)*norm(delta_y,2));
        history.rho_corr(k,:) = [akcor,bkcor];
        if (history.r_norm(k-1) < history.eps_pri(k-1))&&((~is_spectral))
            IS_ADAPTIVE = 0;
            rho_change = 1;
            rho = rho*multiplier;
%             disp('----------------TURN ADAPTIVE OFF----------------')
        elseif (history.r_norm(k-1) > history.eps_pri(k-1))&& ...
                (history.s_norm(k-1) > history.eps_dual(k-1))&&((~is_spectral))% && (abs((history.rho(k-1)-history.rho(k0-1))/history.rho(k0-1))<1e-3)
%             disp('--------------PRIMAL DUAL MAY DIVERGE------------')
            rho = rho*multiplier;
            rho_change = 1;
        elseif (primal_rate<=eps_rate)&&(~is_spectral)
%             disp('---------------PENALTY NOT ADAPTIVE--------------')
            rho = rho*multiplier;
            rho_change = 1;
        elseif ((history.r_norm(k-1) > history.eps_pri(k-1))&&(history.s_norm(k-1) < history.eps_dual(k-1)))&&(~is_spectral||qq==1)
%             disp('----------------PRIMAL INFEASIBLE----------------')
            rho = rho*multiplier;
            rho_change = 1;
        elseif (akcor > epscor) &&  (bkcor>epscor)&&(is_spectral)
            rho = sqrt(ak*bk);
            rho_change = 1;
        elseif (akcor > epscor) &&  (bkcor<=epscor)&&(is_spectral)
            rho = ak;
            rho_change = 1;
        elseif (akcor <= epscor) &&  (bkcor>epscor)&&(is_spectral)
            rho = bk;
            rho_change = 1;
        else
            rho_change = 0;

        end
        
        if (rho_change)&&(is_chol)
            L = chol(GtG+rho*(L1tL1pL2tL2),'lower');
            L = sparse(L); U = L';
            
        end
        
        y1hatk0 = y1hat;
        y2hatk0 = y2hat;
        xk0 = xold;
        z1k0 = z1;
        z2k0 = z2;
        y1k0 = y1;
        y2k0 = y2;
        k0=k;
    end
    
    
    
    
    
    % stopping criterion

    obj = 0.5*norm(G*x-b)^2 + normpq_adaptive(z1,pp,qq,m1,a1)+normpq_adaptive(z2,pp,qq,m2,a2);
    history.nz_count(k) = length(find(z1))+length(find(z2));
    history.sp_difference(k) = length(setdiff(find(z1),find(z1old)))+length(setdiff(find(z1old),find(z1)))+length(setdiff(find(z2),find(z2old)))+length(setdiff(find(z2old),find(z2)));
    history.objval(k) = obj;
    history.r_norm(k) = norm([L1x-z1;L2x-z2]);
    history.s_norm(k)  = norm(rho*At*([z1 - z1old;z2-z2old]));
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
        break;
    end
    
    if (k>2)&&(history.sp_difference(k)==0)
        itr_same_index = itr_same_index+1;
    else
        itr_same_index = 0;
    end
    if (itr_same_index>=MAX_CHECK)
        SPARSITY_CONVERGE = 1;
    end
    if (SPARSITY_CONVERGE)&&(qq==0.5)
        history.fit = 0.5*norm(G*x-b)^2;
        break;
    end
    
    
    %     counting_nz_index = k:-1:max([k-10*Ts,1]);
    %     if (IS_ADAPTIVE==0)&&(PARAMETER.IS_ADAPTIVE==1)&&((sum(history.nz_count(counting_nz_index))-history.nz_count(k)*length(counting_nz_index))<1e-6)
    %         history.fit = 0.5*norm(G*x-b)^2;
    %         break
    %     end
    xold = x;
    z1old = z1;
    z2old = z2;
    y1old = y1;
    y2old = y2;
%     z1hatold = y1hat;
%     z2hatold = y2hat;
end
if (k==MAXITERS)
    history.flag = -1;
else
    history.flag = SPARSITY_CONVERGE;
end

if PRINT_RESULT
    t_end = toc(t_start);
    history.tpi = t_end/k;
end
% return sparse result
switch toggle
    case 'formulationD'
        z1((z2==0)) = 0;
    case 'formulationS'
        % Do nothing because the similarity will eventually shown by
        % constraint LS solver.
%         D = L2;
%         Ind = (1:1:(size(D,2)))';
%         Dplus=D;Dminus=D;
%         Dplus(D==-1) = 0;
%         Dminus(D==1) = 0;
%         Dminus = abs(Dminus);
%         Indplus = Dplus*Ind;
%         Indminus = abs(Dminus*Ind);
%         
%         change = 1;
%         count_change = 0;
%         while change
%             count_change = count_change+1;
%             tmp = x(Indminus(z2==0));
%             x(Indminus(z2==0))=x(Indplus(z2==0));
%             if all(tmp==x(Indminus(z2==0)),'all')
%                 change = 0;
%             end
%         end
%         for fused_index=(z2==0)
%             x(Indminus(fused_index))=x(Indplus(fused_index));
%         end
            
        % the script below extract the fused index out of the vector x

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

