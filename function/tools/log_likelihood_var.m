function [LLH_full,LLH_hetero,LLH_homo,SSE] = log_likelihood_var(data,A,n,p,K)
LLH_full = 0;
LLH_hetero = 0;
LLH_homo = 0;
SSE = 0;
tmpA = reshape(A,[n,n*p,K]);
for kk=1:K
    [H,Y] = H_gen(data(:,:,kk),p);
    Num = size(Y,2);
    Ek = Y - tmpA(:,:,kk)*H; % Error term
    SSE = SSE+sum(Ek.^2,'all');
    Sigma = Ek*Ek'/(Num);
    Sigma = (Sigma+Sigma')/2;
    try
        L = chol(Sigma,'lower');
        logdetSigma = 2*sum(log((diag(L))));
    catch
        logdetSigma = NaN;
    end
    
    %     disp(size(Sigma))
    LLH_full = LLH_full+(Num)*logdetSigma;
    LLH_hetero = LLH_hetero+(Num)*sum(log(abs(diag(Sigma)))); % assume diagonal
    LLH_homo = LLH_homo + Num*n*log(sum(diag(Sigma))/n); % assume multiple of Identity
end
end