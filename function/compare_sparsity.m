function [stat] = compare_sparsity(ind_true,ind,n,K,toggle)

switch toggle
    case 'commonROC'
        ind_true = setdiff(ind_true,1:n+1:n^2);
        GridSize = size(ind,1);
        TP = zeros(GridSize);
        TN = zeros(GridSize);
        FP = zeros(GridSize);
        FN = zeros(GridSize);
        for ii=1:GridSize
            for jj=1:GridSize
                %                 for kk=1:K
                ind{ii,jj} = setdiff(ind{ii,jj},1:n+1:n^2);
                TP(ii,jj) = TP(ii,jj)+length(intersect(ind_true,ind{ii,jj}));
                FN(ii,jj) = FN(ii,jj)+length(setdiff(ind_true,ind{ii,jj}));
                FP(ii,jj) = FP(ii,jj)+length(setdiff(ind{ii,jj},ind_true));
                TN(ii,jj) = n^2-n-TP(ii,jj)-FN(ii,jj)-FP(ii,jj);
                %                 end
            end
        end
    case 'differentialROC'
        GridSize = size(ind,1);
        TP = zeros(GridSize);
        TN = zeros(GridSize);
        FP = zeros(GridSize);
        FN = zeros(GridSize);
        for ii=1:GridSize
            for jj=1:GridSize
                for kk=1:K
                    ind_true = setdiff(ind_true{kk},1:n+1:n^2);
                    ind{ii,jj} = setdiff(ind{ii,jj}{kk},1:n+1:n^2);
                    TP(ii,jj) = TP(ii,jj)+length(intersect(ind_true{kk},ind{ii,jj}{kk}));
                    FN(ii,jj) = FN(ii,jj)+length(setdiff(ind_true{kk},ind{ii,jj}{kk}));
                    FP(ii,jj) = FP(ii,jj)+length(setdiff(ind{ii,jj}{kk},ind_true{kk}));
                end
                TN(ii,jj) = K*(n^2-n)-TP(ii,jj)-FN(ii,jj)-FP(ii,jj);
            end
        end
    case 'single_common'
        ind_true = setdiff(ind_true,1:n+1:n^2);
        ind = setdiff(ind,1:n+1:n^2);
        TP = length(intersect(ind_true,ind));
        FN = length(setdiff(ind_true,ind));
        FP = length(setdiff(ind,ind_true));
        TN = n^2-n-TP-FN-FP;
    case 'single_differential'
        TP = 0;
        FN = 0;
        FP = 0;
        for kk=1:K
            ind_true{kk} = setdiff(ind_true{kk},1:n+1:n^2);
            ind{kk} = setdiff(ind{kk},1:n+1:n^2);
            TP = TP+length(intersect(ind_true{kk},ind{kk}));
            FN = FN+length(setdiff(ind_true{kk},ind{kk}));
            FP = FP+length(setdiff(ind{kk},ind_true{kk}));
        end
        
        TN = K*(n^2-n)-TP-FN-FP;

end
stat(:,:,1) = TP;
stat(:,:,2) = TN;
stat(:,:,3) = FP;
stat(:,:,4) = FN;
stat = squeeze(stat);
end