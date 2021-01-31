function [M_out] = augment_score(M,T,toggle)
name_list = {'bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
tmp = [M.model.stat]; tmp = [tmp.model_selection_score]; tmp=(reshape(([tmp.df]),[M.GridSize,M.GridSize]));
n = size(M.model(1).A,1);
p = size(M.model(1).A,3);
K = size(M.model(1).A,4);
% idx = efficient_vect([n,p,K]);
M_out = M;
for ii=1:M.GridSize
    for jj=1:M.GridSize
%         x_reg = M.model(ii,jj).A_reg(idx);
%         x_cls = M.model(ii,jj).A(idx);
%         tmp_reg_pK = reshape(x_reg,[p*K,n^2]);
%         pK_index = (find(all(tmp_reg_pK~=0,1)));
%         tmp_reg_p = tmp_reg_pK;
%         tmp_reg_p(:,pK_index) = 0;
%         tmp_reg_p = reshape(tmp_reg_p,[p,K*n^2]);
%         tmp_ls_pK = reshape(x_cls,[p*K,n^2]);
%         tmp_ls_pK =sqrt(sum(tmp_ls_pK.^2,1));
%         tmp_reg_pK = sqrt(sum(tmp_reg_pK.^2,1));
%         df_pK = length(find(tmp_reg_pK(pK_index)))+(p*K-1)*sum(tmp_reg_pK(pK_index)./tmp_ls_pK(pK_index),'omitnan');
%         tmp_ls_p = reshape(x_cls,[p,K*n^2]);
%         tmp_ls_p =sqrt(sum(tmp_ls_p.^2,1));
%         tmp_reg_p = sqrt(sum(tmp_reg_p.^2,1));
%         df_p = length(find(tmp_reg_p))+(p-1)*sum(tmp_reg_p./tmp_ls_p,'omitnan');
%         df = df_pK+df_p;
        df=tmp(ii,jj);
        switch toggle
            case 'sse'
                Num = n*K*(T-p);
                SSE = M.model(ii,jj).stat.model_selection_score.SSE;
                fitting = Num*log(SSE/Num);
            case 'llh'
                Num = T-p;
                LLH = M.model(ii,jj).stat.model_selection_score.L;
                fitting = LLH;
        end
        
        
        
        
%         fitting = LLH;
        gamma = log(n^2*p*K)/log(n*(T-p)*K);
        kappa = min([1,1.5*(1-1/(2*gamma))]);
        binom_term = arrayfun(@(x) log_stirling_approx(n^2*p*K)-log_stirling_approx(n^2*p*K-x)-log_stirling_approx(x) , df);
        score.eBIC =   fitting+log(Num)*df + 2*kappa*binom_term;
        score.GIC_2 =  fitting+df*(n^2*p*K)^(1/3);
        score.GIC_3 =  fitting+df*(2*log(n^2*p*K));
        score.GIC_4 =  fitting+df*(2*(log(n^2*p*K)+log(log(n^2*p*K))));
        score.GIC_5 =  fitting+df*log(log(Num))*log(n^2*p*K);
        score.GIC_6 =  fitting+df*log(Num)*log(n^2*p*K);
        score.bic = fitting + log(Num)*df;
        score.aic = fitting + 2*df;
        score.aicc = score.aic + (2*df^2+2*df)/(Num-df-1);
        M_out.model(ii,jj).stat.model_selection_score.df = df;
        for nn=1:length(name_list)
        M_out.model(ii,jj).stat.model_selection_score.(name_list{nn}) =  score.(name_list{nn});
        end
    end
end

for nn=1:length(name_list)
    tmp = [M_out.model.stat]; tmp = [tmp.model_selection_score];tmp=(reshape(([tmp.(name_list{nn})]),[30,30]));
    [~,M_out.index.(name_list{nn})] = min(tmp(:));
end


end