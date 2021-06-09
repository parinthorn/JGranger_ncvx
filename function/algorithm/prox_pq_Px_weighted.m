function r = prox_pq_Px_weighted(z,PInd,v,reg,pp,qq,PARAMETER,gLen)
n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);
% gLen = pK
% gLen = p*K;
% tmpz = z;
r = zeros(n^2*p*K,1);
oo = 1:1:n^2*p*K;
% r(PInd) = prox_pq_eff(z(PInd),v*reg,[pp qq gLen]);
r(PInd) = prox_pq_eff_adaptive(z(PInd),v*reg,[pp qq gLen]);

r((~ismember(oo,PInd))) = z((~ismember(oo,PInd)));
end