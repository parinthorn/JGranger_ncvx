
function [x,cvx_optval] = grouplasso_sharedsp_cvxoff(G, b, P,a1, a2, PARAMETER)

n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);

nn = n^2*p*K;
np = (n^2-n)*p*K;

cvx_begin
variables u(nn,1)
variables v(p*K,n^2-n)
variables w(p,np/p)
minimize ( 0.5*sum((G*u(:)-b).^2)+a1*sum(norms(w,2))+a2*sum(norms(v,2)))
vtmp = P*u(:);
wtmp = P*u(:);
v == reshape(vtmp,p*K,n^2-n);
w == reshape(wtmp,p,np/p);
cvx_end
x = u(:);