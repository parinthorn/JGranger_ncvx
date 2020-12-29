function[D] = genDIFFmat(q,r)

% create a block forward difference matrix of size q(r-1)xqr
% D = [-I  I   0 0 ...   0
%       0 -I   I 0 ...   0
%       0  0  -I I ...   0 
%       .  .   . .       0
%       0  0   0  0  -I  I]
% 
% where I = identity matrix of size q

% D1 = [-sparse(eye(q*(r-1))) zeros(q*(r-1),q)];

D1 = [-speye(q*(r-1)) sparse(q*(r-1),q)];
D2 = [sparse(q*(r-1),q)  speye(q*(r-1))];

D = sparse(D1+D2);

end