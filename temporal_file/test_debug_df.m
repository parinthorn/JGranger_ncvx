%% fix degrees of freedom & model selection score

%input: matrix D, vectorized form x_cls
clear
clc

inpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';

n=20;p=1;K=5;
idx = efficient_vect([n,p,K]);
[P,~] = offdiagJSS(n,p,K);
Dtmp = diffmat(n,p,K);
D = sparse(Dtmp*P);
Ind = (1:1:(size(D,2)))';
Dplus=D;Dminus=D;
Dplus(D==-1) = 0;
Dminus(D==1) = 0;
Dminus = abs(Dminus);
Indplus = Dplus*Ind;
Indminus = abs(Dminus*Ind);

dd=2;
realz = 100;
mname ={'1','5'};
for ii=1:dd
    for jj=1:realz
        load([inpath,'result_formulationS_',mname{ii},'percent','_lag',int2str(p),'_K',int2str(K),'_',int2str(jj)],'M')
        % change M.index,
        % M.model(a1,a2).stat.model_selection_score.bic,aic,aicc,df
        for a1=1:M.GridSize
            for a2=1:M.GridSize
                x_cls = M.model(a1,a2).A(idx);
                fused_index=intersect( ...
                    union(unique(Indplus(Dplus*x_cls~=0)), ...
                    unique(Indminus(Dminus*x_cls~=0))), ...
                    union(Indplus(D*x_cls==0), ...
                    Indminus(D*x_cls==0))); % intuitively, this operation is to find indices of nonzero variables but with zero differences
                %         tmp = (reshape(x_cls(fused_index),[p,length(x_cls(fused_index))/p]));
                tmp =length(find(diff(x_cls(fused_index))==0));
                df = length(find(x_cls))-tmp;
                
            end
        end
        
        
    end
end
