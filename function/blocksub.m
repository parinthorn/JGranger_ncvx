
function[y,IND] = blocksub(x,p,m)
% USAGE: y = blocksub(x,p,m)
% BLOCKSUB performs a linear spacing of step size 'm' in a blockwise manner 
% where each block has 'p' entries
% FOR EXAMPLE, n=20 , p = 2, m = 3
% y = (x1,x2,x6,x7,x11,x12,x16,x17)
% IND : linear indices of the corresponding entries after block spacing

n = length(x);

K = floor(n/(m+p));
if K*(m+p)+p > n % if the last block size exceeds the length of x
    K=K-1;
end

t = (0:K)'*(m+p);
s = 1:p;
IND = repmat(t,1,p)+repmat(s,K+1,1);
IND = reshape(IND',(K+1)*p,1);
y = x(IND);