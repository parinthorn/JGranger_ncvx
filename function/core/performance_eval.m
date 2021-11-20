function model_acc = performance_eval(M,GTmodel)
% This function compute confusion matrix for the estimated model, 
% performance index in all regularization pairs.
% Input
% M: output of jointvargc.m
% GTmodel: Groundtruth model
% Output
% model_acc: performance index [common, differential, total, VAR, estimation error]
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce

n = GTmodel.dim(1);
p = GTmodel.dim(2);
K = GTmodel.dim(3);
GridSize = M.GridSize;
% code below are matlab tricks that extract field from multidimensional
% struct variable. e.g. M(1:2,1:2).A. [M.A] will return [M(1,1).A (2,1).A M(1,2).A M(2,2).A]
% The bracket affect the type of M.A, [M.A] is the same type as A, {M.A}
% collects all [M.A] in cell.
if numel((M.model))==GridSize
    ind_common = reshape([M.model.ind_common],[GridSize 1]);
    ind_differential = reshape([M.model.ind_differential],[GridSize 1]);
    ind_total = reshape([M.model.ind],[GridSize 1]);
    ind_VAR = reshape({M.model.ind_VAR},[GridSize 1]);
    
    for ii=1:GridSize
        
        
        [stat.accuracy.common.confusion_matrix] = compare_sparsity(GTmodel.ind_common,ind_common{ii},n,p,K,'single_common');
        
        [stat.accuracy.differential.confusion_matrix] = compare_sparsity(GTmodel.ind_differential,ind_differential{ii},n,p,K,'single_differential');
        
        
        
        
        [stat.accuracy.total.confusion_matrix] = compare_sparsity(GTmodel.ind,ind_total{ii},n,p,K,'single_differential');
        
        [stat.accuracy.VAR_coeff.confusion_matrix] = compare_sparsity(GTmodel.ind_VAR,ind_VAR{ii},n,p,K,'single_VAR');
        
        
        model_acc(ii).common = performance_score(squeeze(stat.accuracy.common.confusion_matrix));
        
        model_acc(ii).differential = performance_score(squeeze(stat.accuracy.differential.confusion_matrix));
        
        model_acc(ii).total = performance_score(squeeze(stat.accuracy.total.confusion_matrix));
        
        model_acc(ii).VAR_coeff = performance_score(squeeze(stat.accuracy.VAR_coeff.confusion_matrix));
        
        model_acc(ii).bias = sqrt(sum((GTmodel.A-M.model(ii).A).^2,'all')/sum(GTmodel.A.^2,'all'));
    end
    return
end



ind_common = reshape([M.model.ind_common],[GridSize GridSize]);
ind_differential = reshape([M.model.ind_differential],[GridSize GridSize]);
ind_total = reshape([M.model.ind],[GridSize GridSize]);
ind_VAR = reshape({M.model.ind_VAR},[GridSize GridSize]);



for ii=1:GridSize
    for jj=1:GridSize
        if sub2ind([30,30],ii,jj)==M.index.bic
            disp('stop')
        end
        
        [stat.accuracy.common.confusion_matrix] = compare_sparsity(GTmodel.ind_common,ind_common{ii,jj},n,p,K,'single_common');
        
        [stat.accuracy.differential.confusion_matrix] = compare_sparsity(GTmodel.ind_differential,ind_differential{ii,jj},n,p,K,'single_differential');
        
        
        
        
        [stat.accuracy.total.confusion_matrix] = compare_sparsity(GTmodel.ind,ind_total{ii,jj},n,p,K,'single_differential');
        
        [stat.accuracy.VAR_coeff.confusion_matrix] = compare_sparsity(GTmodel.ind_VAR,ind_VAR{ii,jj},n,p,K,'single_VAR');
        
        
        model_acc(ii,jj).common = performance_score(squeeze(stat.accuracy.common.confusion_matrix));
        
        model_acc(ii,jj).differential = performance_score(squeeze(stat.accuracy.differential.confusion_matrix));
        
        model_acc(ii,jj).total = performance_score(squeeze(stat.accuracy.total.confusion_matrix));
        
        model_acc(ii,jj).VAR_coeff = performance_score(squeeze(stat.accuracy.VAR_coeff.confusion_matrix));
        
        model_acc(ii,jj).bias = sqrt(sum((GTmodel.A-M.model(ii,jj).A).^2,'all')/sum(M.model(ii,jj).A.^2,'all'));
    end
end
% stat.bias = squeeze(sqrt(sum(bsxfun(@minus, M.A,GTmodel.A).^2,[1,2,3,4]))./sqrt(sum(model.A.^2,'all')));
end