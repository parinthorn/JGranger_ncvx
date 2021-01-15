function [P_p,P_pK,P1_p,P1_pK,P] = weighted_projection(xLS,q,PARAMETER)
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);
gamma = PARAMETER.gamma;

IND_DIAG = 1:n+1:n^2; % indices of diagonal elements
P1 = speye(n^2); 
P1(IND_DIAG,:) = [];
P = sparse(kron(P1,speye(p*K)));

w1 = (1./(sqrt(sum(reshape(xLS,[p,n^2*K]).^2,1)))).^(q*gamma);

P1_p = sparse(kron(P1,speye(K)))*diag(w1); % project w1 to off block-diagonal size K
P_p = kron(P1_p,speye(p)); % Expand w1 to original size

w2 = (1./(sqrt(sum(reshape(xLS,[p*K,n^2]).^2,1)))).^(q*gamma);
P1_pK = P1*diag(w2); % project w2 to off diagonal block
P_pK = kron(P1_pK,speye(p*K));

% P_m1 = P*sparse(kron(P1_m1,speye(m1)));
% tmp = unique(max(P_m1,[],2),'stable');
% P1_m1(P1_m1~=0) = tmp;

% P_m2 = P*sparse(kron(P1_m2,speye(m2)));
% tmp = unique(max(P_m2,[],2),'stable');
% P1_m2(P1_m2~=0) = tmp;
end