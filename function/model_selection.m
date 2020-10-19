function score = model_selection(Y,A)

[n,~,p,K] = size(A);
df = length(find(A));
fitting = log_likelihood_var(Y,A,n,p,K);
Num = size(Y,2);
score.bic = fitting + log(Num)*df;
score.aic = fitting + 2*df;
score.aicc = score.aic + (2*df^2+2*df)/(Num-df-1);
end
function LLH = log_likelihood_var(data,A,n,p,K)
T = size(data,2);
LLH = 0;
tmpA = reshape(A,[n,n*p,K]);
for kk=1:K
    [H,Y] = H_gen(data(:,:,kk),p);
    Ek = Y - tmpA(:,:,kk)*H; % Error term
    Sigma = Ek*Ek'/(T-p);
    Sigma = (Sigma+Sigma)/2;
    L = chol(Sigma,'lower');
    logdetSigma = 2*sum(log(diag(L)));
%     disp(size(Sigma))
    LLH = LLH+(T-p)*logdetSigma;%+(-1/2)*trace(Ek'*(Sigma\eye(n))*Ek);
end

end