function df = eval_S_df(A)
n = size(A,1);
p= size(A,3);
K = size(A,4);
PARAMETER = [n,p,K];
IND_DIAG = 1:n+1:n^2; % indices of diagonal elements
P1 = speye(n^2); 
P1(IND_DIAG,:) = [];
P = sparse(kron(P1,speye(p*K)));
Dtmp = diffmat(n,p,K);
D = sparse(Dtmp*P);

Ind = (1:1:(size(D,2)))';
Dplus=D;Dminus=D;
Dplus(D==-1) = 0;
Dminus(D==1) = 0;
Dminus = abs(Dminus);
Indplus = Dplus*Ind;
Indminus = abs(Dminus*Ind);
idx = efficient_vect(PARAMETER);
x = A(idx);
fused_index=intersect( ...
    union(unique(Indplus(Dplus*x~=0)), ...
    unique(Indminus(Dminus*x~=0))), ...
    union(Indplus(D*x==0), ...
    Indminus(D*x==0))); % intuitively, this operation is to find indices of nonzero variables but with zero differences
%         tmp = (reshape(x_cls(fused_index),[p,length(x_cls(fused_index))/p]));
tmp =length(find(diff(x(fused_index))==0));
df = length(find(x))-tmp;

end