function [x,obj_prox] = select_prox(z,a1,a2,qq)
prox_obj_fun = @(x) a2*norm(x(:),'fro')^(qq)+a1*sum(sqrt(sum(x.^2,1)).^(qq))+ ...
    0.5*norm(x(:)-z(:),'fro')^2;
[p,K] = size(z);
reference_loss = prox_obj_fun(zeros(p,K));
h = zeros(p,K);
for kk=1:K % choose kk from K
    ind = combnk(1:K,kk); % all combination
    obj_prox_val = zeros(size(ind,1),1);
    x = zeros(p,K,size(ind,1));
    flag = zeros(size(ind,1),1);
    for jj = 1:size(ind,1)
        nz_ind = ind(jj,:);
        
        [x(:,:,jj),prox_obj,flag(jj)] = fix_point_prox(z,nz_ind,a1,a2,qq);
        
        obj_prox_val(jj) =prox_obj;
    end
    [~,choose_ind] = min(obj_prox_val);
    if obj_prox_val(choose_ind)<reference_loss
        h = x(:,:,choose_ind);
        reference_loss = prox_obj_fun(h);
    end
end
x = h;
obj_prox = prox_obj;

end