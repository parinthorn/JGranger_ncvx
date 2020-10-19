function [x] = sim_VAR(A,T,noise_var,seed,is_double)
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
K = size(A,4);

% optvargin = size(varargin,2);

if nargin<3
    noise_var = 1;
    seed = ('default');
    is_double = 1;
elseif nargin<4
    seed = ('default');
    is_double = 1;
elseif nargin<5
    is_double = 1;
elseif nargin>6
    error('exceed number of input')
end
rng(seed);
if is_double
    u = sqrt(noise_var)*randn(n,T);
    x = zeros(n,T,K);
    matA = reshape(A,[n,n*p,K]);
else
    u = (single(sqrt(noise_var)*randn(n,T)));
    x = (single(zeros(n,T,K)));
    matA = (single(reshape(A,[n,n*p,K])));
end

for kk=1:K
    AK = matA(:,:,kk);
for tt=p+1:T
    x(:,tt,kk) = AK*reshape(x(:,tt-1:-1:tt-p),n*p,1);
    x(:,tt,kk) = x(:,tt) + u(:,tt);
end
end
x = x-mean(x,2);
