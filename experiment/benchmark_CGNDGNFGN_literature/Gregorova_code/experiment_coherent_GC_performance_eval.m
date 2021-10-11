clear
clc
inpath = './experiment/model_parameters/';
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_C_magda/';
outpath = './results2plot/';
mname = {'10','20'};
realization = 100;
load([inpath,'model_K5_p1'])
R.F1 = zeros(length(mname),realization);
R.MCC = zeros(length(mname),realization);
R.TPR = zeros(length(mname),realization);
R.FPR = zeros(length(mname),realization);
R.ACC = zeros(length(mname),realization);
for ii=1:length(mname)
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{2,ii+2,2,jj};
        fname = [resultpath,'estim_magda_',mname{ii},'percent_',int2str(jj)];
        load(fname)
        model_acc = performance_eval(M,GTmodel);
        R.F1(ii,jj) =model_acc(M.index.eBIC).common.F1;
        R.MCC(ii,jj) =model_acc(M.index.eBIC).common.MCC;
        R.TPR(ii,jj) =model_acc(M.index.eBIC).common.TPR;
        R.FPR(ii,jj) =model_acc(M.index.eBIC).common.FPR;
        R.ACC(ii,jj) =model_acc(M.index.eBIC).common.ACC;
        R.bias(ii,jj) =model_acc(M.index.eBIC).bias;
               fprintf(' F1 avg:%.3f \n MCC avg:%.3f \n ACC avg:%.3f \n FPR avg:%.3f \n TPR avg:%.3f \n', ...
            mean(R.F1(ii,1:jj)),mean(R.MCC(ii,1:jj)),mean(R.ACC(ii,1:jj)),mean(R.FPR(ii,1:jj)),mean(R.TPR(ii,1:jj)))
 
    end
end

% save([resultpath,'magda_result'],'R')
save([outpath,'magda_result'],'R')