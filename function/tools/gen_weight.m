function [weight_1,weight_2] = gen_weight(xLS,PARAMETER,q,formulation)
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);
gamma = PARAMETER.gamma;
IND_DIAG = 1:n+1:n^2; % indices of diagonal elements
P1 = speye(n^2);
P1(IND_DIAG,:) = [];
if ~strcmp(formulation,'fgn')
    
    w1 = (1./(sqrt(sum(reshape(xLS,[p,n^2*K]).^2,1)))).^(q*gamma);
    tmp_matrix = sparse((length(w1)),(length(w1)));
    tmp_matrix(1:length(w1)+1:length(w1)^2)=w1;
    P1_p = sparse(kron(P1,speye(K)))*tmp_matrix; % project w1 to off block-diagonal size K
    w2 = (1./(sqrt(sum(reshape(xLS,[p*K,n^2]).^2,1)))).^(q*gamma);
    tmp_matrix = sparse((length(w2)),(length(w2)));
    tmp_matrix(1:length(w2)+1:length(w2)^2)=w2;
    P1_pK = P1*tmp_matrix; % project w2 to off diagonal block
    weight_1 = max(P1_p,[],2);
    weight_2 = max(P1_pK,[],2);
else
    P = sparse(kron(P1,speye(p*K)));
    Dtmp = diffmat(n,p,K);
    D = sparse(Dtmp*P);
    w1 = (1./(sqrt(sum(reshape(xLS,[p,n^2*K]).^2,1)))).^(q*gamma);
    tmp_matrix = sparse((length(w1)),(length(w1)));
    tmp_matrix(1:length(w1)+1:length(w1)^2)=w1;
    P1_p = sparse(kron(P1,speye(K)))*tmp_matrix; % project w1 to off block-diagonal size K
    weight_1 =  max(P1_p,[],2);
    w2 = (1./(sqrt(sum(reshape(D*xLS,[p,(n^2-n)*nchoosek(K,2)]).^2,1)))).^(q*gamma);
    weight_2 = w2(:);
    
end

end