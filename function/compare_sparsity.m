function [stat] = compare_sparsity(ind_true,ind,n,K,toggle)
switch toggle
    case 'commonROC'
        GridSize = size(ind,1);
        TP = zeros(GridSize);
        TN = zeros(GridSize);
        FP = zeros(GridSize);
        FN = zeros(GridSize);
        for ii=1:GridSize
            for jj=1:GridSize
                for kk=1:K
                    TP(ii,jj) = TP(ii,jj)+length(intersect(ind_true,ind{ii,jj}));
                    FN(ii,jj) = FN(ii,jj)+length(setdiff(ind_true,ind{ii,jj}));
                    FP(ii,jj) = FP(ii,jj)+length(setdiff(ind{ii,jj},ind_true));
                    TN(ii,jj) = n^2-n-TP(ii,jj)-FN(ii,jj)-FP(ii,jj);
                end
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
                    TP(ii,jj) = TP(ii,jj)+length(intersect(ind_true{kk},ind{ii,jj}{kk}));
                    FN(ii,jj) = FN(ii,jj)+length(setdiff(ind_true{kk},ind{ii,jj}{kk}));
                    FP(ii,jj) = FP(ii,jj)+length(setdiff(ind{ii,jj}{kk},ind_true{kk}));
                    TN(ii,jj) = n^2-n-TP(ii,jj)-FN(ii,jj)-FP(ii,jj);
                end
            end
        end
    case 'single_common'
        TP = 0;
        FN = 0;
        FP = 0;
        TN = 0;
        for kk=1:K
            TP = TP+length(intersect(ind_true,ind));
            FN = FN+length(setdiff(ind_true,ind));
            FP = FP+length(setdiff(ind,ind_true));
            TN = TN+n^2-n-TP-FN-FP;
        end
    case 'single_differential'
        TP = 0;
        FN = 0;
        FP = 0;
        TN = 0;
        for kk=1:K
            TP = TP+length(intersect(ind_true{kk},ind{kk}));
            FN = FN+length(setdiff(ind_true{kk},ind{kk}));
            FP = FP+length(setdiff(ind{kk},ind_true{kk}));
            TN = TN+n^2-n-TP-FN-FP;
        end
end
stat(:,:,1) = TP;
stat(:,:,2) = TN;
stat(:,:,3) = FP;
stat(:,:,4) = FN;
stat = squeeze(stat);
end