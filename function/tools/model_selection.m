function score = model_selection(Y,A)
[n,~,p,K] = size(A);
df = length(find(A));
fitting = log_likelihood_var(Y,A,n,p,K);
Num = size(Y,2);

gamma = log(n^2*p*K)/log(n*(Num));
kappa = 1.5*(1-1/(2*gamma));
binom_term = arrayfun(@(x) log_stirling_approx(n^2*p*K)-log_stirling_approx(n^2*p*K-x)-log_stirling_approx(x) , df);
score.eBIC = fitting+log(Num).*df + 2*kappa.*binom_term;
score.GIC_2 =  fitting+df*(n^2*p*K)^(1/3);
score.GIC_3 =  fitting+df*(2*log(n^2*p*K));
score.GIC_4 =  fitting+df*(2*(log(n^2*p*K)+log(log(n^2*p*K))));
score.GIC_5 = fitting+df*log(log(Num)).*log(n^2*p*K);
score.GIC_6 = fitting+df*log(Num).*log(n^2*p*K);

score.bic = fitting + log(Num)*df;
score.aic = fitting + 2*df;
score.aicc = score.aic + (2*df^2+2*df)/(Num-df-1);
score.L = fitting;
score.df = df;
end
function LLH = log_likelihood_var(data,A,n,p,K)
Num = size(data,2);
LLH = 0;
tmpA = reshape(A,[n,n*p,K]);
for kk=1:K
    [H,Y] = H_gen(data(:,:,kk),p);
    Ek = Y - tmpA(:,:,kk)*H; % Error term
    Sigma = Ek*Ek'/(Num);
    Sigma = (Sigma+Sigma')/2;
    L = chol(Sigma,'lower');
    logdetSigma = 2*sum(log(diag(L)));
    %     disp(size(Sigma))
    LLH = LLH+(Num)*logdetSigma;%+(-1/2)*trace(Ek'*(Sigma\eye(n))*Ek);
end

end