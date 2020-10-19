function r = prox_pq_Px_efficient(z,PInd,v,reg,pp,qq,PARAMETER,gLen)
n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);
% gLen = pK
% gLen = p*K;
% tmpz = z;
r = zeros(n^2*p*K,1);
oo = 1:1:n^2*p*K;
r(PInd) = prox_pq_eff(z(PInd),v*reg,[pp qq gLen]);

r((~ismember(oo,PInd))) = z((~ismember(oo,PInd)));
end
% 
% function prox_x = prox_pq_Internal(z,v,PARAMETERS)
% % This function computes proximal operator 
% % of composite norm pq reg*(||x||_p,q)^q
% % for some p, q
% p = PARAMETERS(1);
% q = PARAMETERS(2);
% gLen = PARAMETERS(3);
% n = length(z);
% gNo = n/gLen;
% tmp_prox = zeros(gLen,gNo);
% normz_vect = ((sum((reshape(z,[gLen,gNo])).^2,1))).^(0.5)';
% if (p==2 && q==0.5)
% gTag = (normz_vect)>1.5*v^(2/3);
% ngTag = (normz_vect)<=1.5*v^(2/3);
% coeff_num = @(zz) (16.*zz.^(1.5).*(cos(pi/3-(1/3)*(acos(v/4*(3./zz).^(3/2))))).^3);
% coeff_den = @(zz) (3*sqrt(3)*v+16.*zz.^(1.5).*cos(pi/3-(1/3)*(acos(v/4*(3./zz).^(3/2)))).^3);
% tmp_prox(gTag) = (coeff_num(normz_vect(gTag))./coeff_den(normz_vect(gTag)));
% tmp_prox(ngTag) = 0;
% tmp_prox = kron(tmp_prox',ones(gLen,1));
% prox_x = tmp_prox.*z;
% elseif (p==2 && q==1)
%     gTag = (normz_vect)>v;
%     ngTag = (normz_vect)<=v;
%     tmp_prox(gTag) = (1-v./(normz_vect(gTag)));
%     tmp_prox(ngTag) = 0;
%     tmp_prox = kron(tmp_prox',ones(gLen,1));
%     prox_x = tmp_prox.*z;
% end
% end
