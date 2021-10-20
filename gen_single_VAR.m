function [A,ind,VAR_spectrum,ind_VAR] = gen_single_VAR(n,p,VAR_MATRIX_DENSITY,ADJUST_VAR_COEFF)
THRESH = ADJUST_VAR_COEFF(1);
MIN_VAR = ADJUST_VAR_COEFF(2);
S = sprand(n,n,VAR_MATRIX_DENSITY)+eye(n);
UNSTABLE = 1;
diag_ind = find(eye(n));
k = length(diag_ind);
diag_ind3D = kron(n^2*(0:p-1)',ones(k,1))+kron(ones(p,1),diag_ind);
off_ind3D = (1-repmat(eye(n),1,1,p));
A = zeros(n,n,p);
ii = 0;
while UNSTABLE
    ii = ii+1;
    for k=1:p
        A(:,:,k) = 0.1*sprandn(S(:,:,1));
    end
    poles = -0.7+2*0.7*rand(n,p); % make the poles inside the unit circle
    characeq = zeros(n,p+1);
    for jj=1:n
        characeq(jj,:) = poly(poles(jj,:)); % each row is [1 -a1 -a2 ... -ap]
    end
    aux = -characeq(:,2:end);
    A(diag_ind3D) = aux(:); % replace the diagonal entries with stable polynomial
    A(((abs(A)<THRESH)&(A~=0))&(off_ind3D))=MIN_VAR.*sign(A(((abs(A)<THRESH)&(A~=0))&(off_ind3D)));
    AA = zeros(n,n*p);
    for k=1:p
        AA(1:n,k*n-(n-1):k*n) = A(:,:,k);
    end
    AA = sparse([AA ; [eye(n*(p-1)) zeros(n*(p-1),n)]]);
    if (abs(eigs(AA,1)) < 1)
        VAR_spectrum = eig(full(AA));
        UNSTABLE = 0;
    elseif ii>100
        error("CANNOT FIND A STABLE MODEL")
    end
end
ind = find(A(:,:,1));
ind_VAR = find(A);
end
