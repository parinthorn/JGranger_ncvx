function score = group_lasso_model_selection(Y,A,A_reg,toggle)
[n,~,p,K] = size(A);


switch toggle
    case 'yuanlin' %2006
        idx = efficient_vect([n,p,K]);
        x_reg = A_reg(idx);
        x_cls = A(idx);
        tmp_reg = reshape(x_reg,[p,n^2*K]);
        tmp_ls = reshape(x_cls,[p,n^2*K]);
        tmp_ls =sqrt(sum(tmp_ls.^2,1));
        tmp_reg = sqrt(sum(tmp_reg.^2,1));
        df = length(find(tmp_reg))+(p-1)*sum(tmp_reg./tmp_ls,'omitnan');
%     case 'brehenyhuang' % 2009
%     case 'vaiter' %2012
    case 'lasso'
        df = length(find(A));
end




fitting = log_likelihood_var(Y,A,n,p,K);
Num = size(Y,2);
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