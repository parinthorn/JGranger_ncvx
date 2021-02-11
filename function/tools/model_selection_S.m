function score = model_selection_S(Y,A,df)
[n,~,p,K] = size(A);
[score.LLH_full,score.LLH_hetero,score.LLH_homo,score.SSE] = log_likelihood_var(Y,A,n,p,K);% fitting = -2*Log-Likelihood
df_lasso = length(find(A));
Num = size(Y,2);
gamma = log(n^2*p*K)/log(n*(Num));
kappa = 1.5*(1-1/(2*gamma));
if df==n^2*p*K
    binom_term = 0;
else
    binom_term = arrayfun(@(x) log_stirling_approx(n^2*p*K)-log_stirling_approx(n^2*p*K-x)-log_stirling_approx(x) , df);
end
fitting = score.LLH_full;
score.eBIC = fitting+log(Num).*df + 2*kappa.*binom_term;
score.GIC_2 =  fitting+df*(n^2*p*K)^(1/3);
score.GIC_3 =  fitting+df*(2*log(n^2*p*K));
score.GIC_4 =  fitting+df*(2*(log(n^2*p*K)+log(log(n^2*p*K))));
score.GIC_5 = fitting+df*log(log(Num)).*log(n^2*p*K);
score.GIC_6 = fitting+df*log(Num).*log(n^2*p*K);

score.bic_lasso = fitting + log(Num)*df_lasso;
score.bic = fitting + log(Num)*df;
score.aic = fitting + 2*df;
score.aicc = score.aic + (2*df^2+2*df)/(Num-df-1);
score.df = df;
score.df_lasso = df_lasso;
end
function [LLH_full,LLH_hetero,LLH_homo,SSE] = log_likelihood_var(data,A,n,p,K)
Num = size(data,2);
LLH_full = 0;
LLH_hetero = 0;
LLH_homo = 0;
SSE = 0;
tmpA = reshape(A,[n,n*p,K]);
for kk=1:K
    [H,Y] = H_gen(data(:,:,kk),p);
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
    LLH_hetero = LLH_hetero+(Num)*sum(log(abs(diag(Sigma))));
    LLH_homo = LLH_homo + Num*sum(diag(Sigma))/n;
end
end