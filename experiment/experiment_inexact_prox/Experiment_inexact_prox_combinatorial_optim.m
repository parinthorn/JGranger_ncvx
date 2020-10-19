clear
clc
clf
for tt=1:1
p=10;
K = 5;
q = 0.5;
a1 = 5;a2 = 2;
z = 1+1*randn(p,K);
%%
PARAMETERS  = [2,q,p];
tic;
prox_x = prox_pq_eff(z(:),a1,PARAMETERS);
toc;
prox_x = reshape(prox_x,[p,K]);


prox_obj_fun = @(x) a2*norm(x(:),'fro')^(q)+a1*sum(norms(x,2,1).^(q))+ ...
    0.5*norm(x(:)-z(:),'fro')^2;
% prox_obj_fun = @(x) a2*norm(x,2)^(q)+a1*sum(norms(reshape(x,[p,K]),2,1).^(q))+ ...
%     0.5*norm(x(:)-z(:),2)^2;
reference_loss = prox_obj_fun(zeros(p*K,1)); % begin loss with zero solution
zero_loss = prox_obj_fun(zeros(p*K,1));
h = zeros(p,K);
t1 = tic;
for kk=1:K % choose kk from K
    ind = combnk(1:K,kk); % all combination
%     x_tmp = h;
    obj_prox_val = zeros(size(ind,1),1);
    x = zeros(p,K,size(ind,1));
    for jj = 1:size(ind,1)
        nz_ind = ind(jj,:);
        [x(:,:,jj),prox_obj,obj_sse] = fix_point_prox(z,nz_ind,a1,a2,q);
        obj_prox_val(jj) =prox_obj(end);
    end
    [tmp,choose_ind] = min(obj_prox_val);
%     disp(tmp)
    if obj_prox_val(choose_ind)<=reference_loss
%         h = zeros(p,K);
        h = x(:,:,choose_ind);
%         h(:,ind(choose_ind,:)) = x(:,ind(choose_ind,:));
        reference_loss = prox_obj_fun(h);
%         disp(reference_loss)

    end
end
t_combi = toc(t1);
% disp(norm(h-prox_x,'fro'))
x_combi = h;
%% compare combinatorial with linear time greedy algorithm
t1 = tic;
[x_greedy,prox_obj_greedy,obj_sse_greedy] = fix_point_prox_post(z,1:1:K,a1,a2,q);
% Find Sparse solution
reference_loss = prox_obj_fun(x_greedy);
% zero_loss = obj_prox(zeros(p*K,1));
z_ind = [];
h = x_greedy;
for kk=1:K
    tmp = h(:,kk);
    h(:,kk) = 0;
    if prox_obj_fun(h)<reference_loss
        z_ind = [z_ind kk];
        h(:,kk) = tmp;
    else
        h(:,kk) = tmp;
    end
end
nz_ind = setdiff(1:K,z_ind);
h=x_greedy;
h(:,z_ind) = 0;

if zero_loss<prox_obj_fun(h)
    h = zeros(p,K);
    z_ind = 1:1:K;
    x_greedy = h;
else
    [x_greedy,prox_obj_greedy,obj_sse_greedy] = fix_point_prox(z,setdiff(1:K,z_ind),a1,a2,q);
end

t_greedy = toc(t1);
% x_greedy = h;
sol_diff(tt) = norm(x_greedy-x_combi,'fro')/norm(x_combi,'fro');
fprintf('realization: %d, relative diff %.3f \n',tt,sol_diff(tt))
if sol_diff(tt)>0.1
    error('investigate this')
end
[X,FVAL] = fminsearch(prox_obj_fun,[ones(p,K)]);
end

