function [y,G,x] = vectorize_VAR(Y,H,PARAMETER)

n = PARAMETER(1); p = PARAMETER(2); K = PARAMETER(3);
N = PARAMETER(4);

%% Create G

G = [];
blkH = zeros(N,n*p,K);BLKH = zeros(N,n*p*K,K);
E1 = zeros(K,1); E1(1) = 1;
for k=1:K,
    H0 = zeros(n,N,p);
    for j=1:p,
        H0(:,:,j) = H((j-1)*n+1:j*n,:,k);
    end
    ind = linindex(n,N,p,'col');
    vecH = H0(ind);

    TMP = zeros(n*p*K*N,1);
    [tmp2,IND] = blocksub(TMP,p,p*(K-1));
    TMP(IND) = vecH;
    TMP = reshape(TMP,n*p*K,N);
    TMP = TMP';

    BLKH(:,:,k) = sparse([zeros(N,(k-1)*p) TMP(:,1:end-(k-1)*p)]);
    BLKG = kron(speye(n),BLKH(:,:,k)); % not efficient when n is large
    % size of BLKG is nN x pn^2K
    G = [G;BLKG];
end

G = sparse(G);

%% Create y
ind = linindex(n,n,p,'row');
y =zeros(n*N,K);
for k=1:K,
    Ek = zeros(K,1); Ek(k) = 1;
    y(:,k) = reshape(Y(:,:,k)',n*N,1);
end
y = reshape(y,[n*N*K,1]);
