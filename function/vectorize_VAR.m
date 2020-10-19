function [y,G] = vectorize_VAR(Y,H,PARAMETER)

n = PARAMETER(1); p = PARAMETER(2); K = PARAMETER(3);
T = PARAMETER(4);
Num = T-p;
%% Create G

G = [];
blkH = zeros(Num,n*p,K);BLKH = zeros(Num,n*p*K,K);
E1 = zeros(K,1); E1(1) = 1;
for k=1:K,
    H0 = zeros(n,Num,p);
    for j=1:p,
        H0(:,:,j) = H((j-1)*n+1:j*n,:,k);
    end
    ind = linindex(n,Num,p,'col');
    vecH = H0(ind);

    TMP = zeros(n*p*K*Num,1);
    [tmp2,IND] = blocksub(TMP,p,p*(K-1));
    TMP(IND) = vecH;
    TMP = reshape(TMP,n*p*K,Num);
    TMP = TMP';

    BLKH(:,:,k) = sparse([zeros(Num,(k-1)*p) TMP(:,1:end-(k-1)*p)]);
    BLKG = kron(speye(n),BLKH(:,:,k)); % not efficient when n is large
    % size of BLKG is nN x pn^2K
    G = [G;BLKG];
end

G = sparse(G);

%% Create y
ind = linindex(n,n,p,'row');
y =zeros(n*Num,K);
for k=1:K,
    Ek = zeros(K,1); Ek(k) = 1;
    y(:,k) = reshape(Y(:,:,k)',n*Num,1);
end
y = reshape(y,[n*Num*K,1]);
