
function[y,G,x] = veccoefmat(Y,H,A)
% 
% VECCOEFMAT vectorizes the AR coefficients used in the LS problem:
% minimize || Y - AH ||_2^2 = || y - Gx ||_2^2
% where A is n x np, H is np x N, Y is n x N.
% A = [A1 A2 ... Ap], H = [H1; H2; ...; Hp], Y = [y1 y2 ... yN]
% After the vectorization,  y is nN x 1, x is pn^2 x 1, G is nN x pn^2
% (see the description of y, x, G from the paper)
% USAGE [y,G,x] = veccoefmat(Y,H,A);

[n,m] = size(A);
if mod(m,n)~= 0,
    error('Check the dimension of A');
else
    p = floor(m/n);
end
[r,s] = size(H);
if r~= m
    error('The dimensions of A and H must agree');
else
    N = s;
end
[q,t] = size(Y);
if q~=n
    error('The dimensions of Y and A must agree');
elseif t~=N
    error('The dimensions of Y and H must agree');
end

y = reshape(Y',n*N,1);
A0 = reshape(A,n,n,p);
ind = linindex(n,n,p,'row');
x = A0(ind);

H0 = zeros(n,N,p);
for k=1:p,
    H0(:,:,k) = H((k-1)*n+1:k*n,:);
end
ind = linindex(n,N,p,'col');
vecH = H0(ind);

tmp = reshape(vecH,n*p,N);
blkH = tmp';

G = kron(speye(n),blkH); % not efficient when n is large
G = sparse(G);

tmp1 = Y-A*H; tmp2 = y-G*x;
if (norm(vec(tmp1))-norm(tmp2))>= 1e-3,
    error('There is some error in the calculation.');
end