function [P_p,P,P1_p,D,w1,w2] = weighted_projection_S(xLS,q,PARAMETER)
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);
gamma = PARAMETER.gamma;

IND_DIAG = 1:n+1:n^2; % indices of diagonal elements
P1 = speye(n^2); 
P1(IND_DIAG,:) = [];
P = sparse(kron(P1,speye(p*K)));
Dtmp = diffmat(n,p,K);
D = sparse(Dtmp*P);

w1 = (1./(sqrt(sum(reshape(xLS,[p,n^2*K]).^2,1)))).^(q*gamma);
P1_p = sparse(kron(P1,speye(K)))*diag(w1); % project w1 to off block-diagonal size K

w1 = sparse(kron(P1,speye(K)))*w1(:);

P_p = kron(P1_p,speye(p)); % Expand w1 to original size

w2 = (1./(sqrt(sum(reshape(D*xLS,[p,(n^2-n)*nchoosek(K,2)]).^2,1)))).^(q*gamma);
w2 = w2(:);
end