function [x,obj_prox,obj_sse] = fix_point_prox_post(z,nz_ind,a1,a2,q)
[p,K] = size(z);
% size(nz_ind) = K
MAX_ITER = 5000;
x0 = 0.1*ones(p,K);
% x0(:,nz_ind) = 0.1;
x = x0;
sse_fun = @(x,kk) norm((a2*q*sqrt(norm(x,'fro')^2)^(q-2)+ ...
    a1*q*sqrt(norm(x(:,kk),'fro')^2)^(q-2)+1)*(x(:,kk))-z(:,kk),'fro')^2; % sum square error of optimality condition
prox_obj_fun = @(x) a2*norm(x(:),'fro')^(q)+a1*sum(norms(x,2,1).^(q))+ ...
    0.5*norm(x(:)-z(:),2)^2;
% prox_obj_fun = @(x) a2*norm(x,2)^(q)+a1*sum(norms(reshape(x,[p,K]),2,1).^(q))+ ...
%     0.5*norm(x(:)-z(:),2)^2;
obj_sse = zeros(MAX_ITER,1);
obj_prox = zeros(MAX_ITER,1);
for tt=1:MAX_ITER
%     x_old = x;
    for kk=1:K % update only given nonzero group
        rcp = a2*q*(norm(x(:),2))^(q-2)+a1*q*(norm(x(:,kk),2))^(q-2);
        x(:,kk) = (z(:,kk)-rcp*x(:,kk));
    end
    
    for kk=1:K
        obj_sse(tt) = obj_sse(tt) + sse_fun(x,kk);
    end
    obj_prox(tt) = prox_obj_fun(x);
    if (tt>2) && (abs((obj_sse(tt) - obj_sse(tt-1))) <1e-6)
%         disp('equal')
        obj_sse = obj_sse(obj_sse~=0);
        obj_prox = obj_prox(obj_prox~=0);
        break
    end 
end
x(:,setdiff(1:K,nz_ind)) = 0;
obj_prox = prox_obj_fun(x);
end