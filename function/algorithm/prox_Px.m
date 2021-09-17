function r = prox_Px(z,a1,a2,PARAMETER)
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);
pp = PARAMETER.p;
qq = PARAMETER.q;
PInd = PARAMETER.proj_ind;
if (a1==0)&&(a2~=0)
    gLen = p*K;
elseif (a1~=0)&&(a2==0)
    gLen = p;
elseif (a1==0)&&(a2==0)
    r=z;
    return;
else
    error('no closed-form solution for general value of a1, a2')
end


r = zeros(n^2*p*K,1);
oo = 1:1:n^2*p*K;
r(PInd) = prox_pq_eff(z(PInd),a,[pp qq gLen]);

r((~ismember(oo,PInd))) = z((~ismember(oo,PInd)));
end