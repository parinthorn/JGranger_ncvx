
function [x] = prox_sumof2norm(u,p,a)
% 
% PROX_SUMOF2NORM computes the proximal operator of a*f where
% 
% f(x) = sum_{k=1}^K ||xk ||_2
% where xk is p x 1 and 'a' is a scalar
% 
% The proximal operator of f is the block soft thresholding
% 
% prox_{af}(u) _{kth block} = max(1- a/||uk||_2, 0 )* uk
% 
% USAGE: [x] = prox_sumof2norm(u,p,a)
% u = (u1,u2,...,uN) where uk has size p x 1
%


n = length(u);
if mod(n,p)~= 0
    error('mod(n,p) must be zero. Enter a new p');
else 
    M = floor(n/p);
end
z = reshape(u,p,M); % z = [u1 u2 ... uM]
w = norms(z,2); % use norms in CVX

z = z.*repmat(pos(1 - a./w),p,1);
x = z(:);
