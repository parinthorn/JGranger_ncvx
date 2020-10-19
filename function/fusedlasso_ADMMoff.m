function[x,Dx,Px, history] = fusedlasso_ADMMoff(G, b, P,D,a1, a2, PARAMETER, rho,Num,varargin)

n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);

t_start = tic;

PRINT_RESULT = 1;
FREQ_PRINT = 100;
MAXITERS = 50000; 
ABSTOL = 1e-7;
RELTOL = 1e-5;

% store variables
nn = n^2*p*K; nd = (n^2-n)*p*(K-1); np = (n^2-n)*p*K;
Gtb = G'*b;
A = [P ; D];
At= A';
Dt = D';
Pt = P';
DtD = sparse(D'*D); 
GtG = sparse(G'*G);
PtP = sparse(P'*P);
L = chol(GtG+rho*PtP+rho*DtD,'lower'); 
L = sparse(L); U = L';

% ADMM solver

if PRINT_RESULT
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n','iter', ...
        'r norm','eps pri','s norm','eps dual', 'objective');
end

% initial start 
optargin = size(varargin,2);
if optargin == 0,
    x1 = zeros(nn,1);
else
    x1 = varargin{1};
end

if a2 == 0,
    [x,Px,history] = grouplasso(G, b, P,a1, [n p 1],rho,Num);
end

z1 = zeros(np,1); %z2 = zeros(nd,1);
z2 = zeros((n*n-n)*nchoosek(K,2)*p,1);
Dx1 = D*x1;
Px1 = P*x1;

x3 = Dx1;
x2 = Px1;

obj = 0.5*norm(G*x1-b)^2/Num + a1*norm21(x2,p)+a2*norm21(x3,p);

for k=1:MAXITERS,
    
    % x1-update
    
    c = Gtb+ Pt*(rho*x2-z1) + Dt*(rho*x3-z2); 
    x1 = U \ (L \ c);
    
    % x2-update
    x2old = x2;
    Px1 = P*x1;
    x2 = prox_sumof2norm(Px1+z1/rho,p,a1/rho);
    
    
    % x3-update
    x3old = x3;
    Dx1 = D*x1;
    x3 = prox_sumof2norm(Dx1+z2/rho,p,a2/rho);
    
    % z-update
    z1 = z1 + rho*(Px1-x2);
    z2 = z2 + rho*(Dx1-x3);
    
    
    % stopping criterion
    
    obj = 0.5*norm(G*x1-b)^2/Num+a1*norm21(x2,p)+a2*norm21(x3,p);
  
    history.objval(k) = obj;  
    history.r_norm(k) = norm([Px1-x2;Dx1-x3]);
    history.s_norm(k)  = norm(rho*At*([x2 - x2old;x3-x3old]));    
      
    history.eps_pri(k) = sqrt(np+nd)*ABSTOL + RELTOL*max(norm(A*x1), norm([x2;x3]));
    history.eps_dual(k)= sqrt(nn)*ABSTOL + RELTOL*norm(At*[z1;z2]);
%     history.x{k} = x2;
    
    if (PRINT_RESULT && mod(k,FREQ_PRINT) == 0)
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n',k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), ...
            history.objval(k));
    end
    
    if (history.r_norm(k) < history.eps_pri(k) && ...
            history.s_norm(k) < history.eps_dual(k) )
        break;
    end
%         if (history.r_norm(k) < ABSTOL)
%             break;
%         end

end

% if PRINT_RESULT
%     toc(t_start);
% end
 toc(t_start);
% return sparse result
X = reshape(x1,p,K,n^2); % x1 is not sparse
X2 = reshape(x2,p,K,n^2-n); % x2 is offdiagonal entries only and sparse
IND_offdiag = setdiff((1:n^2)',(1:n+1:n^2)','rows');
X(:,:,IND_offdiag) = X2;
x = X(:);

Dx = x3;
Px = x2;
end



