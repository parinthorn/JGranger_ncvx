clear
clc
inpath = './data_compare/';
resultpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\';
mname = {'1','5'};
realization = 20;
load([inpath,'model_K5_p1'])
R.F1 = zeros(length(mname),realization);
R.MCC = zeros(length(mname),realization);
R.TPR = zeros(length(mname),realization);
R.FPR = zeros(length(mname),realization);
R.ACC = zeros(length(mname),realization);
for ii=1:length(mname)
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{2,3,ii,jj};
        fname = [resultpath,'result_formulationD_',mname{ii},'percent_lag1_K5_',int2str(jj)];
        load(fname)
        model_acc = performance_eval(M,GTmodel);
        R.F1(ii,jj) =model_acc(M.index.bic).total.F1;
        R.MCC(ii,jj) =model_acc(M.index.bic).total.MCC;
        R.TPR(ii,jj) =model_acc(M.index.bic).total.TPR;
        R.FPR(ii,jj) =model_acc(M.index.bic).total.FPR;
        R.ACC(ii,jj) =model_acc(M.index.bic).total.ACC;
        
%          model_acc(ii,jj).result = performance_eval(M,GTmodel);
%          model_acc(ii,jj).M = M;
        
    end
end
save([resultpath,'formulation_D_result_total'],'R')