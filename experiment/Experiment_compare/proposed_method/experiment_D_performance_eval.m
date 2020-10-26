clear
clc
inpath = './data_compare/';
mname = {'1','5'};
realization = 20;
load([inpath,'model_K5_p1'])
F1_avg = zeros(2,1);
for ii=2:length(mname)
    for jj=1:realization
        GTmodel = E{2,3,ii,jj};
        fname = ['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\result_formulationD_',mname{ii},'percent_lag1_K5_',int2str(jj)];
        load(fname)
        model_acc = performance_eval(M,GTmodel);
        disp(model_acc(M.index.bic).total.F1)
        F1_avg(ii) = F1_avg(ii)+model_acc(M.index.bic).total.F1/realization;
    end
end