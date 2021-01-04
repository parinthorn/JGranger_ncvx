
function [x,cvx_optval] = fusedlasso_cvxoff(G, b, P,D,a1, a2, PARAMETER)

n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);

nn = n^2*p*K;
nd = (n^2-n)*p*(K-1);
np = (n^2-n)*p*K;

cvx_begin
variables u(p,nn/p)
variables v(p,nd/p)
variables w(p,np/p)
minimize ( 0.5*sum((G*u(:)-b).^2)+a1*sum(norms(w,2))+a2*sum(norms(v,2)))
vtmp = D*u(:);
wtmp = P*u(:);
v == reshape(vtmp,p,nd/p);
w == reshape(wtmp,p,np/p);
cvx_end
x = u(:);