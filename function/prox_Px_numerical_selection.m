function r = prox_Px_numerical_selection(z,PARAMETER,selection_type)

PInd = PARAMETER.proj_ind;
a1 = PARAMETER.a1;
a2 = PARAMETER.a2;
pp = PARAMETER.p;
qq = PARAMETER.q;
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);
if (a1==0) || (a2==0)
    r = prox_Px(z,PARAMETER);
    return
end
r = zeros(n^2*p*K,1);
tmp_out = zeros(p,K,n^2-n);
oo = 1:1:n^2*p*K;
tmp = z(PInd);
tmp = reshape(tmp,[p,K,n^2-n]);
x0 = reshape(prox_pq_eff(z(PInd),a2,[2 qq p*K]),[p*K,n^2-n]); % add comment
nz_col = find(x0(1,:));
% tic
for ii= nz_col
    tmp_out(:,:,ii) = select_prox_selection(tmp(:,:,ii),a1,a2,qq,selection_type);
end
% toc;
r(PInd) = tmp_out(:);

r((~ismember(oo,PInd))) = z((~ismember(oo,PInd)));
end