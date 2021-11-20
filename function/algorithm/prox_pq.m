function prox_x = prox_pq(z,v,reg,PARAMETERS)
% This function computes proximal operator 
% of composite norm pq reg*(||x||_p,q)^q
% for some p, q
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce
p = PARAMETERS(1);
q = PARAMETERS(2);
gLen = PARAMETERS(3);
n = length(z);
gNo = n/gLen;
tmp_prox = zeros(gLen,gNo);
for LL=1:gNo
    tmp_prox(:,LL) = prox_norm_pq(z(gLen*(LL-1)+1:gLen*LL),p,q,reg,v);
end
prox_x = real(tmp_prox(:));
    function x = prox_norm_pq(z,p,q,reg,v)
        normz = norm(z,2);
        dimz = size(z,1);
        vL = v*reg;
        Q = @(p,q,x) reg*norm(x,p)^(q)+(1/(2*v))*norm(x-z,2)^2;
        if (p==2 && q==1)
            if normz>vL
                x = (1-vL/normz)*z;
            else
                x = zeros(dimz,1);
            end
        elseif (p==2 && q==0)
            if normz>sqrt(2*vL)
                x=z;
            else
                x=zeros(dimz,1);
            end
        elseif (p==2 && q==0.5)
            if normz>(1.5*(vL)^(2/3))
                phi = @(zz) acos(vL/4*(3/norm(zz,2))^(3/2));
                c = (cos(pi/3-phi(z)/3))^3;
                x = (16*normz^(3/2)*c)/(3*sqrt(3)*vL+16*normz^(3/2)*c)*z;
            else
                x=zeros(dimz,1);
            end
        elseif (p==1 && q==0.5)
            E = @(z) acos(vL*dimz/4*(3/norm(z,1))^(3/2));
            ztilde = z- (sqrt(3)*vL/(4*sqrt(norm(z,1))*cos(pi/3-E(z)/3))).*sign(z);
            if Q(1,0.5,ztilde)<Q(1,0.5,0)
                x = ztilde;
            else
                x = zeros(dimz,1);
            end
        end
    end
end

