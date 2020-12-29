function [x,H,Y] = gen_timeseries(A,Num,varargin)
% generate x from x(:,t)=A(1)*x(:,t-1)+...+A(p)*x(:,t-p)+noise
% x = the generated timeseries
%   = [x(:,1) x(:,2) ... x(:,Num)]
% H = [y(p)     y(p+1)  ...     y(N-1)
%      y(p-1)   y(p)    ...     y(N-2)
%       ...     ...     ...     ...
%      y(1)     y(2)    ...     y(N-p)]
% Y = x(:,p+1:N)
%       [x,H,Y] = gent_timeseries(phi,Num)
% A = the coefficient in 3D-array (A(:,:,1),A(:,:,2),...,A(:,:,p))
% Num = amount of data in the generated timeseries

n = size(A,1);
p = size(A,3);

optvargin = size(varargin,2);

if optvargin == 0,
    noise_var=1;
else
    noise_var = varargin{1};
end

u = sqrt(noise_var)*randn(n,Num);

x = zeros(n,Num);
matA = reshape(A,n,n*p);
for k=p+1:Num,
    x(:,k) = matA*reshape(x(:,k-1:-1:k-p),n*p,1);
    x(:,k) = x(:,k) + u(:,k);
end

x = detrend(x','constant')';

H = zeros(n*p,Num-p);
for k=1:p,
    H((k-1)*n+1:k*n,:) = reshape(x(:,p-k+1:Num-k),n,Num-p);
end
Y = x(:,p+1:end);
