function [x,prox_objval,flag,obj_prox_array] = fix_point_prox(z,nz_ind,a1,a2,q)
tol = 1e-3;
[p,K] = size(z);
% size(nz_ind) = K
MAX_ITER = 100000;
SAFEGUARD_STEP = 20;
OSCILLATE_CHECK_START = 50;
x0 = zeros(p,K);
x0(:,nz_ind) = 0.1;
x = x0;
prox_obj_fun = @(x) a2*norm(x(:),'fro')^(q)+a1*sum(sqrt(sum(x.^2,1)).^(q))+ ... 
    0.5*norm(x(:)-z(:),2)^2;
obj_prox_array = zeros(MAX_ITER,1);

tmp_vec = zeros(p,K,MAX_ITER);
for tt=1:MAX_ITER
    tmp_vec(:,:,tt) = x;
    for kk=nz_ind % update only given nonzero group
        rcp = a2*q*(norm(x(:),2))^(q-2)+a1*q*(norm(x(:,kk),2))^(q-2);
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
%         d1 = obj_prox(tt-18:tt-3);
%         d2 = flipud(obj_prox(tt-16:tt-1));
%         std_prox = std(obj_prox(tt-10:tt));
%         cond = std_prox/mean(obj_prox(tt-10:tt));
        if any(std_vec>1e-2,'all')
%             plot(obj_prox_array(1:tt))
            
%             for kk=1:K
%                 subplot(1,K,kk)
%                 plot(squeeze(tmp_vec(:,kk,SAFEGUARD_STEP:tt))')
%             end
%             pause(0.1)
            flag = -1;
            x = nan;
            obj_prox_array = inf;
            prox_objval = obj_prox_array(end);
            return
        end
    end
end
% for kk=1:K
%     subplot(1,K,kk)
%     plot(squeeze(tmp_vec(:,kk,SAFEGUARD_STEP:tt))')
% end
% pause(0.1)
flag = -1;
end