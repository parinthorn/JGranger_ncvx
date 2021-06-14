function[matP,matD] = constructPD(PARAMETER)
% USAGE: PARAMETER.n = 3; PARAMETER.p = 2; PARAMETER.K = 10;
%   [matP,matD] = constructPD(PARAMETER)
% matP is a structure consisting of P and P'P
% matD is a structure consisting of D, D'D, and D1 ( check kronecker
% product of Dk with I_{pk} in the paper.
% 
% Notes: 
% 
% 0) If 'n' or 'K' are sufficiently, like n >= 100, K >= 100 (at the same
% time), it has memory limit problem.
% 1) Kronecker product of two sparse matrices are efficient in most cases.
% 2) Computing quadratic product (like P'P) directly is more efficient 
% than assigning values to certain indices (if implemented in MATLAB).
% 3) Creating a sparse matrix from 'sparse' command should be used but not always fastest.
% 4) Saving sparse arrays takes much less space than saving indices of sparse matrix 
% (and hope to exploit indexing instead; this hypothesis was false). MATLAB
% manages sparse arrays better than us handling the indices.
% 5) Sometimes, if matrices A, B are already sparse, no need to command 'sparse( kron(A,B))' 
%     or sparse ( operations of A,B) again. It takes some time to compute.
% MATLAB returns sparse structure of operations of two sparse matrices. 

n = PARAMETER.n ;
p = PARAMETER.p ;
K = PARAMETER.K ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct P and P'P
IND_DIAG = 1:n+1:n^2; % indices of diagonal elements (n x n)
P1 = speye(n^2);  P1(IND_DIAG,:) = [];
P = sparse(kron(P1,speye(p*K))); 
matP.P = P;
matP.PtP = P'*P;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct D and D'D

% The goal is to create (I,J) indices of element 1 and -1 and use 'sparse' command to create Dk

IJplus = zeros(nchoosek(K,2),2); % (I,J) indices of 1
IJminus = zeros(nchoosek(K,2),2); % (I,J) indices of -1
rowind_init = 0;
for colind = 1:(K-1)
    tmp = repmat(colind,(K-colind),1); % 
    rowind = rowind_init + (1: (K-colind))';    
    IJplus(rowind,:) = [rowind tmp] ;
    IJminus(rowind,:) = [rowind (colind+1:K)'];
    rowind_init = rowind(end);
end
numrowDk = nchoosek(K,2);  plus1 = ones(numrowDk,1); neg1 = -ones(numrowDk,1);
Dk = sparse([IJplus(:,1) ; IJminus(:,1)],[IJplus(:,2) ; IJminus(:,2)], [plus1 ; neg1],numrowDk,K); % it's better to construct sparse matrices first
D1 = kron(Dk,speye(p)); % pK x pK 
D = kron(P1,D1);

% In fact, D is created just for computing D'D, so we can remove it. 
matD.D = D;
matD.DtD = D'*D;
matD.D1 = sparse(D1); % used when computing Dx; see the following code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% The following code is to compute Dx (needs only P and D1)
% numrowDk = nchoosek(K,2);
% x = randn(n^2*p*K,1);
% z = P*x; blockz = reshape(z,p*K,n^2-n); % split z into blocks each of size pK
% D1z = D1*blockz; % [D1*block 1 , D1*block 2 , ... ]
% y = reshape(D1z,numrowDk*p*(n^2-n),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes for using these matrices in optimization: 
% P and D1 are used for iterative methods (when computing Px, Dx)
% P'P and D'D are cached only once. In fact, we should cache P'D+D'D 