function M = fix_loglikelihood(M_in,data,toggle)
[n,~,p,K] = size(M_in.model(1).A);
for kk=1:K
    [H(:,:,kk),Y(:,:,kk)] = H_gen(data(:,:,kk),size(M_in.model(1).A,3));
end
GridSize = M_in.GridSize ;
M = M_in;
if prod(size(M.model))==M.GridSize
    for ii=1:GridSize
        A = M.model(ii).A;
        df = M.model(ii).stat.model_selection_score.df;
        [M.model(ii).stat.model_selection_score.LLH_full,M.model(ii).stat.model_selection_score.LLH_hetero,M.model(ii).stat.model_selection_score.LLH_homo,M.model(ii).stat.model_selection_score.SSE] = log_likelihood_var(Y,A,n,p,K);
        Num = size(Y,2);
        
        gamma = log(n^2*p*K)/log(n*(Num)*K);
        kappa = min([1,1.5*(1-1/(2*gamma))]);
        if df==n^2*p*K
            binom_term = 0;
        else
            binom_term = arrayfun(@(x) log_stirling_approx(n^2*p*K)-log_stirling_approx(n^2*p*K-x)-log_stirling_approx(x) , df);
        end
        fitting = M.model(ii).stat.model_selection_score.LLH_full; % default
        M.model(ii).stat.model_selection_score.flag = 0;
        M.model(ii).stat.model_selection_score.eBIC =   fitting+log(Num)*df + 2*kappa*binom_term;
        M.model(ii).stat.model_selection_score.GIC_2 =  fitting+df*(n^2*p*K)^(1/3);
        M.model(ii).stat.model_selection_score.GIC_3 =  fitting+df*(2*log(n^2*p*K));
        M.model(ii).stat.model_selection_score.GIC_4 =  fitting+df*(2*(log(n^2*p*K)+log(log(n^2*p*K))));
        M.model(ii).stat.model_selection_score.GIC_5 =  fitting+df*log(log(Num))*log(n^2*p*K);
        M.model(ii).stat.model_selection_score.GIC_6 =  fitting+df*log(Num)*log(n^2*p*K);
        
        M.model(ii).stat.model_selection_score.bic = fitting + log(Num)*df;
        M.model(ii).stat.model_selection_score.aic = fitting + 2*df;
        M.model(ii).stat.model_selection_score.aicc = M.model(ii).stat.model_selection_score.aic + (2*df^2+2*df)/(Num-df-1);
        
        M.model(ii).stat.model_selection_score.df = df;
    end
    return;
end


for ii=1:GridSize
    for jj=1:GridSize
        A = M.model(ii,jj).A;
        df = length(find(A));
        [M.model(ii,jj).stat.model_selection_score.LLH_full,M.model(ii,jj).stat.model_selection_score.LLH_hetero,M.model(ii,jj).stat.model_selection_score.LLH_homo,M.model(ii,jj).stat.model_selection_score.SSE] = log_likelihood_var(Y,H,A,n,p,K);
        Num = size(Y,2);
        
        gamma = log(n^2*p*K)/log(n*(Num)*K);
        kappa = min([1,1.5*(1-1/(2*gamma))]);
        if df==n^2*p*K
            binom_term = 0;
        else
            binom_term = arrayfun(@(x) log_stirling_approx(n^2*p*K)-log_stirling_approx(n^2*p*K-x)-log_stirling_approx(x) , df);
        end
        fitting = M.model(ii,jj).stat.model_selection_score.(toggle); % default
        M.model(ii,jj).stat.model_selection_score.flag = 0;
        M.model(ii,jj).stat.model_selection_score.eBIC =   fitting+log(Num)*df + 2*kappa*binom_term;
        M.model(ii,jj).stat.model_selection_score.GIC_2 =  fitting+df*(n^2*p*K)^(1/3);
        M.model(ii,jj).stat.model_selection_score.GIC_3 =  fitting+df*(2*log(n^2*p*K));
        M.model(ii,jj).stat.model_selection_score.GIC_4 =  fitting+df*(2*(log(n^2*p*K)+log(log(n^2*p*K))));
        M.model(ii,jj).stat.model_selection_score.GIC_5 =  fitting+df*log(log(Num))*log(n^2*p*K);
        M.model(ii,jj).stat.model_selection_score.GIC_6 =  fitting+df*log(Num)*log(n^2*p*K);
        
        M.model(ii,jj).stat.model_selection_score.bic = fitting + log(Num)*df;
        M.model(ii,jj).stat.model_selection_score.aic = fitting + 2*df;
        M.model(ii,jj).stat.model_selection_score.aicc = M.model(ii,jj).stat.model_selection_score.aic + (2*df^2+2*df)/(Num-df-1);
        
        M.model(ii,jj).stat.model_selection_score.df = df;
    end
end
end