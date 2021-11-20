function prox_x = prox_pq_eff_adaptive(z,v,PARAMETERS)
% This function computes proximal operator
% of composite norm pq v*(||x||_p,q)^q
% for p=2, q=1 or p=2, q=0.5
%
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce
p = PARAMETERS(1);
q = PARAMETERS(2);
gLen = PARAMETERS(3);
n = length(z);
gNo = n/gLen;
normz_vect = ((sum((reshape(z,[gLen,gNo])).^2,1))).^(0.5)';
if length(v)==1
    if (p==2 && q==0.5)
        nz_ind = (normz_vect)>1.5*v^(2/3);
        z_ind = (normz_vect)<=1.5*v^(2/3);
        coeff_num = @(zz) (16.*zz.^(1.5).*(cos(pi/3-(1/3)*(acos(v/4*(3./zz).^(3/2))))).^3);
        coeff_den = @(zz) (3*sqrt(3)*v+16.*zz.^(1.5).*cos(pi/3-(1/3)*(acos(v/4*(3./zz).^(3/2)))).^3);
        tmp_prox(nz_ind) = (coeff_num(normz_vect(nz_ind))./coeff_den(normz_vect(nz_ind)));
        tmp_prox(z_ind) = 0;
        tmp_prox = kron(tmp_prox',ones(gLen,1));
        prox_x = tmp_prox.*z;
    elseif (p==2 && q==1)
        nz_ind = (normz_vect)>v;
        z_ind = (normz_vect)<=v;
        tmp_prox(nz_ind) = (1-v./(normz_vect(nz_ind)));
        tmp_prox(z_ind) = 0;
        tmp_prox = kron(tmp_prox',ones(gLen,1));
        prox_x = tmp_prox.*z;
    end
else
    if length(v)~=gNo
        error('adaptive penalization dimension error')
    end
    
    
    if (p==2 && q==0.5)
        nz_ind = (normz_vect)>1.5.*v.^(2/3);
        z_ind = (normz_vect)<=1.5.*v.^(2/3);
        v = v(nz_ind);
        coeff_num = @(zz) (16.*zz.^(1.5).*(cos(pi/3-(1/3).*(acos(v./4.*(3./zz).^(3/2))))).^3);
        coeff_den = @(zz) (3.*sqrt(3).*v+16.*zz.^(1.5).*cos(pi/3-(1/3).*(acos(v./4.*(3./zz).^(3/2)))).^3);
        tmp_prox(nz_ind) = (coeff_num(normz_vect(nz_ind))./coeff_den(normz_vect(nz_ind)));
        tmp_prox(z_ind) = 0;
        tmp_prox = kron(tmp_prox',ones(gLen,1));
        prox_x = tmp_prox.*z;
    elseif (p==2 && q==1)
        nz_ind = (normz_vect)>v;
        z_ind = (normz_vect)<=v;
        v = v(nz_ind);
        tmp_prox(nz_ind) = (1-v./(normz_vect(nz_ind)));
        tmp_prox(z_ind) = 0;
        tmp_prox = kron(tmp_prox',ones(gLen,1));
        prox_x = tmp_prox.*z;
    end
end
end

