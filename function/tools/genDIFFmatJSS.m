function[D] = genDIFFmatJSS(n,p,k)

% create a block forward difference matrix of size q(r-1)xqr
% D = [-I  I   0 0 ...   0
%       0 -I   I 0 ...   0
%       0  0  -I I ...   0 
%       .  .   . .       0
%       0  0   0  0  -I  I]
% 
% where I = identity matrix of size p

% D1 = [-sparse(eye(q*(r-1))) zeros(q*(r-1),q)];


% D1 = [-speye(p*(k-1)) sparse(p*(k-1),p)];
% D2 = [sparse(p*(k-1),p)  speye(p*(k-1))];
% 
% D = sparse(D1+D2);
tmp = speye(n^2-n);
Da = [speye(k-1) zeros(k-1,1)];
Db = [zeros(k-1,1) speye(k-1)];
tmpD = sparse(Da-Db);
DD = kron(tmpD,speye(p));
D = kron(tmp,DD);

end