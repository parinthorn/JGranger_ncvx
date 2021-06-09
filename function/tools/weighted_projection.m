function [P_p,P_pK,P1_p,P1_pK,P] = weighted_projection(xLS,q,PARAMETER,gamma)
n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);

IND_DIAG = 1:n+1:n^2; % indices of diagonal elements
P1 = speye(n^2); 
P1(IND_DIAG,:) = [];
P = sparse(kron(P1,speye(p*K)));

w1 = (1./(sqrt(sum(reshape(xLS,[p,n^2*K]).^2,1)))).^(q*gamma);
tmp_matrix = sparse((length(w1)),(length(w1)));
tmp_matrix(1:length(w1)+1:length(w1)^2)=w1;
P1_p = sparse(kron(P1,speye(K)))*tmp_matrix; % project w1 to off block-diagonal size K
P_p = kron(P1_p,speye(p)); % Expand w1 to original size

w2 = (1./(sqrt(sum(reshape(xLS,[p*K,n^2]).^2,1)))).^(q*gamma);
tmp_matrix = sparse((length(w2)),(length(w2)));
tmp_matrix(1:length(w2)+1:length(w2)^2)=w2;
P1_pK = P1*tmp_matrix; % project w2 to off diagonal block
P_pK = kron(P1_pK,speye(p*K));

% P_m1 = P*sparse(kron(P1_m1,speye(m1)));
% tmp = unique(max(P_m1,[],2),'stable');
% P1_m1(P1_m1~=0) = tmp;

% P_m2 = P*sparse(kron(P1_m2,speye(m2)));
% tmp = unique(max(P_m2,[],2),'stable');
% P1_m2(P1_m2~=0) = tmp;
end