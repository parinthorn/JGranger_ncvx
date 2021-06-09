function s = compare_sparsity(ind_true,ind)
TP = length(intersect(ind_true,ind));
FP = length(setdiff(ind,ind_true));
FN = length(setdiff(ind_true,ind));
s.F1 = TP/(TP+0.5*(FP+FN));
end