function [ind_nz,A,VAR_spectrum,seed] = gen_VAR(n,p,density,diffdensity,GENERATION_TYPE,S)
  THRESH = 0.4;
  GAIN = 0.4;
seed = randi(1000000);
rng(seed)
if GENERATION_TYPE==2
  S = sprand(n,n,density)+eye(n);
end
if (nargin<5)
  S = sprand(n,n,density)+eye(n);
end
%% Randomize AR coefficients
if (GENERATION_TYPE==0)||(GENERATION_TYPE==2) % Generate common structure with given matrix 'S'
    MAX_EIG = 1;
    diag_ind = find(eye(n));
    k = length(diag_ind);
    diag_ind3D = kron(n^2*(0:p-1)',ones(k,1))+kron(ones(p,1),diag_ind);
    off_ind3D = (1-repmat(eye(n),1,1,p));
    A = zeros(n,n,p);
    ii = 0;
    while MAX_EIG
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
        A(((abs(A)<THRESH)&(A~=0))&(off_ind3D))=GAIN.*sign(A(((abs(A)<THRESH)&(A~=0))&(off_ind3D)));
        AA = zeros(n,n*p);
        for k=1:p
            AA(1:n,k*n-(n-1):k*n) = A(:,:,k);
        end
        AA = sparse([AA ; [eye(n*(p-1)) zeros(n*(p-1),n)]]);
        if (max(abs(eigs(AA))) < 1)
            VAR_spectrum = eig(full(AA));
            MAX_EIG = 0;
        end
    end
elseif (GENERATION_TYPE==1) % Generate Similar model
    MAX_EIG = 1;
    off_ind3D = (1-repmat(eye(n),1,1,p));
    diag_ind = (1:n+1:n^2)';
    A = zeros(n,n,p);
    R = sprand(n,n,diffdensity);
    R(diag_ind) = 0;
    R(S(:,:,1)~=0) = 0;
    ii = 0;
    while MAX_EIG
        ii = ii+1;
        for k=1:p
            A(:,:,k) = S(:,:,k)+ 0.2*sprandn(R); % add differential part
        end
        A(((abs(A)<THRESH)&(A~=0))&(off_ind3D))=GAIN.*sign(A(((abs(A)<THRESH)&(A~=0))&(off_ind3D)));
        AA = zeros(n,n*p);
        for k=1:p
            AA(1:n,k*n-(n-1):k*n) = A(:,:,k);
        end
        AA = sparse([AA ; [speye(n*(p-1)) zeros(n*(p-1),n)]]);
        if max(abs(eigs(AA))) < 1
            VAR_spectrum = eig(full(AA));
            MAX_EIG = 0;
        end
    end
end
ind_nz = find(A(:,:,1));
end
