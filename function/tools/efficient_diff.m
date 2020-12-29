function [idx1,idx2] = efficient_diff(PARAMETER)
n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);
[P,~] = offdiagJSS(n,p,K);
D = sparse(diffmat(n,p,K))*P;
x=(1:1:n^2*p*K)';
D1 = D;
D2 = D;
D1(D==-1) = 0;
D2(D==1) = 0;
D2 = (-1)*D2;
idx1 = D1*x;
idx2 = D2*x;
end