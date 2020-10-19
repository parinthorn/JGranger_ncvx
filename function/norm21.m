function[y,w] = norm21(x,m)
% NORM21 return the group norm of an n-dimensional vector x  
% where x = (x1,x2,...,xM) and x1,x2,...,xM are vectors of size m
%     In other words, we chop x in to k subvectors and each subvector has size m
% 
% Therefore, mod(n,m) must be zero
% 
% y = NORM21(x,m) = sum_{j=1}^M || xj ||_2
% w = [ ||x1|| ||x2|| ... ||xM|| ]


n = length(x);
if mod(n,m)~= 0
    error('mod(n,m) must be zero. Enter a new m');
else 
    M = floor(n/m);
end

if m==1
    y = norm(x,1); % typical l1-norm
    w = abs(x);
else
    z = reshape(x,m,M); % z = [x1 x2 ... xM]
    w = norms(z,2)'; % use norms in CVX
    y = sum(w); 
end