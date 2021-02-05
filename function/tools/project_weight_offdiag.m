function [P,P1_p,P1_pK] = project_weight_offdiag(PARAMETER)
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);

IND_DIAG = 1:n+1:n^2; % indices of diagonal elements
P1_pK = speye(n^2); 
P1_pK(IND_DIAG,:) = [];
P = sparse(kron(P1_pK,speye(p*K)));

P1_p = sparse(kron(P1_pK,speye(K))); % project w1 to off block-diagonal size K
end