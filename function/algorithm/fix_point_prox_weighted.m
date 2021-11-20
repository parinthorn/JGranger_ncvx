function [x,prox_objval,flag,obj_prox_array] = fix_point_prox_weighted(z,nz_ind,w1,w2,q)
% This function numerically solve the problem
% min_{x} g(x) + 1/2v ||x-z||_2^2
% where g(x) is reported in the thesis[1], with a given sparsity of x
% Input
% z: input vector
% nz_ind: given sparsity of x
% w1: parameter of g(x)
% w2: parameter of g(x)
% Output
% x: converged solution
% prox_objval: converged objective of proximal operator
% flag: convergence flag 0 is converge, -1 is not converge
% obj_prox_array: history of proximal objective
% 
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce
tol = 1e-3;
[p,K] = size(z);
% size(nz_ind) = K
MAX_ITER = 100000;
SAFEGUARD_STEP = 20;
OSCILLATE_CHECK_START = 50;
x0 = zeros(p,K);
x0(:,nz_ind) = 0.1;
x = x0;
prox_obj_fun = @(x) w2*norm(x(:),'fro')^(q)+sum(w1'.*sqrt(sum(x.^2,1)).^(q))+ ...
    0.5*norm(x(:)-z(:),'fro')^2;
obj_prox_array = zeros(MAX_ITER,1);

tmp_vec = zeros(p,K,MAX_ITER);
for tt=1:MAX_ITER
    tmp_vec(:,:,tt) = x;
    for kk=nz_ind % update only given nonzero group
        rcp = w2*q*(norm(x(:),2))^(q-2)+w1(kk)*q*(norm(x(:,kk),2))^(q-2);
        x(:,kk) = (z(:,kk)-rcp*x(:,kk));
    end
    obj_prox_array(tt) = prox_obj_fun(x);
    if (tt>2) && ((abs((obj_prox_array(tt) - obj_prox_array(tt-1))/obj_prox_array(tt-1)) <tol))
        obj_prox_array = obj_prox_array(1:tt);
        flag = 0;
        prox_objval = obj_prox_array(end);
        return
    elseif (tt>OSCILLATE_CHECK_START)
        std_vec = std(tmp_vec(:,:,SAFEGUARD_STEP:tt),0,3);
        mean_vec = mean(tmp_vec(:,:,SAFEGUARD_STEP:tt),3);
        std_vec(mean_vec~=0) = std_vec(mean_vec~=0)./abs(mean_vec((mean_vec~=0)));
        if any(std_vec>1e-2,'all')
            flag = -1;
            x = nan;
            obj_prox_array = inf;
            prox_objval = obj_prox_array(end);
            return
        end
    end
end
flag = -1;
end