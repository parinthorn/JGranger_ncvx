function [ind] = linindex(m,n,p,type)
% LININDEX gives the linear indices of the vectorized 3D array
% Z = (Z1,Z2,...,Zp) where Zk is m x n
% if type == 'row' ind is the linear indices of
% [(Z1)_{11} (Z2)_{11} ... (Zp)_{11}
% (Z1)_{12} (Z2)_{12} ... (Zp)_{12} ...
% (Z1)_{1n} (Z2)_{1n} ... (Zp)_{1n}
% (Z1)_{21} (Z2)_{21} ... (Zp)_{21} ...
% (Z1)_{2n} (Z2)_{2n} ... (Zp)_{2n} ...
% (Z1)_{m1} (Z2)_{m1} ... (Zp)_{m1} ...
% (Z1)_{mn} (Z2)_{mn} ... (Zp)_{mn}
% ]
% if type == 'col' ind is the linear indices of
% [(Z1)_{11} (Z2)_{11} ... (Zp)_{11}
% (Z1)_{21} (Z2)_{21} ... (Zp)_{21} ...
% (Z1)_{m1} (Z2)_{m1} ... (Zp)_{m1}
% (Z1)_{12} (Z2)_{12} ... (Zp)_{12} ...
% (Z1)_{m2} (Z2)_{m2} ... (Zp)_{m2} ...
% (Z1)_{1m} (Z2)_{1m} ... (Zp)_{1m} ...
% (Z1)_{mn} (Z2)_{mn} ... (Zp)_{mn}
% ]
K = repmat([1:p]',m*n,1); % same as K = kron(ones(m*n,1),[1:p]');
switch type
case 'row'
I = kron([1:m]',ones(n*p,1));
Jtmp = kron([1:n]',ones(p,1));
J = repmat(Jtmp,m,1); % same as J = kron(ones(m,1),Jtmp);
ind= sub2ind([m n p],I,J,K);
case 'col'
Itmp = kron([1:m]',ones(p,1));
I = repmat(Itmp,n,1); % same as I = kron(ones(n,1),Itmp);
J = kron([1:n]',ones(m*p,1));
ind = sub2ind([m n p],I,J,K);
end