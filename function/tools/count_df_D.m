function df = count_df_D(x_reg,x_cls,PARAMETER)
n=  PARAMETER(1);
p = PARAMETER(2);
K=PARAMETER(3);
tmp_reg_pK = reshape(x_reg,[p*K,n^2]);
pK_index = (find(all(tmp_reg_pK~=0,1)));
tmp_reg_p = tmp_reg_pK;
tmp_reg_p(:,pK_index) = 0;
tmp_reg_p = reshape(tmp_reg_p,[p,K*n^2]);
tmp_reg_p = sqrt(sum(tmp_reg_p.^2,1));
tmp_ls_pK = reshape(x_cls,[p*K,n^2]);
tmp_ls_pK =sqrt(sum(tmp_ls_pK.^2,1));
tmp_reg_pK = sqrt(sum(tmp_reg_pK.^2,1));
df_pK = length(find(tmp_reg_pK(pK_index)))+(p*K-1)*sum(tmp_reg_pK(pK_index)./tmp_ls_pK(pK_index),'omitnan');
tmp_ls_p = reshape(x_cls,[p,K*n^2]);
tmp_ls_p =sqrt(sum(tmp_ls_p.^2,1));
df_p = length(find(tmp_reg_p))+(p-1)*sum(tmp_reg_p./tmp_ls_p,'omitnan');
df = df_pK+df_p;
        
end