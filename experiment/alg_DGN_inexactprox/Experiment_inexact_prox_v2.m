clear
clc
clf
p=50;
K = 10;

q = 0.5;
x0 = 0.1*ones(p*K,1);
x0 = reshape(x0,[p,K]);
z = 1+1*randn(p,K);

T = 2000;
obj_equal = zeros(T,1);


%% Parinthorn fixed point Derivation [CONVERGE]
a1 = 18;
a2 = 0;
% obj_prox = @(x) a2*norm(x,2)^(q)+a1*sum(norms(reshape(x,[p,K]),2,1).^(q))+ ...
%     0.5*norm(x(:)-z(:),2)^2;
obj_prox = @(x) a2*norm(x,2)^(q)+a1*sum(sqrt(sum(reshape(x,[p,K]).^2,1)).^(q))+ ...
    0.5*norm(x(:)-z(:),2)^2;
equality_obj = @(x,kk) norm((a2*q*sqrt(norm(x,'fro')^2)^(q-2)+ ...
    a1*q*sqrt(norm(x(:,kk),'fro')^2)^(q-2)+1)*(x(:,kk))-z(:,kk),2)^2;
PARAMETERS  = [2,q,p];
prox_x = prox_pq_eff(z(:),a1,PARAMETERS);
x = x0;
for tt=1:T
    for kk=1:K
        rcp = a2*q*(norm(x(:),2))^(q-2)+a1*q*(norm(x(:,kk),2))^(q-2);
        x(:,kk) = (z(:,kk)-rcp*x(:,kk));
    end
    
    for kk=1:K
        obj_equal(tt) = obj_equal(tt) + equality_obj(x,kk);
    end
    obj(tt) = obj_prox(x(:));
    if (tt>2) && (abs((obj_equal(tt) - obj_equal(tt-1))) <1e-6)
        disp('equal')
        break
    end 
end
obj_equal = obj_equal(1:tt);
%% test



%%