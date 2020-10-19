function D = diffmat(n,p,k)
% This function create difference transformation matrix
% which transforms vectorized VAR parameter into differenced of all groups
% in all off-diagonal
tmp = eye(k);
tmpD = zeros(nchoosek(k,2),k);
r=0;
for i=1:k
    for j=i+1:k
        r = r+1;
        tmpD(r,:) = tmp(i,:)-tmp(j,:);
    end
end
tmpD = sparse(kron(tmpD,eye(p)));
D = sparse(kron(speye(n^2-n),tmpD));
% disp(size(D))
end