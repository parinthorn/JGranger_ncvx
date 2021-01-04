clear
clc
inpath = './data_compare/';
resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_S_result/';
mname = {'1','5'};
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
        GTmodel = E{3,3,ii,jj};
        fname = [resultpath,'result_formulationS_',mname{ii},'percent_lag1_K5_',int2str(jj)];
        load(fname)
        model_acc = performance_eval(M,GTmodel);
        toggle_list = {'total','common','differential'};
        ALL_RESULT(ii,jj).model_acc = model_acc;
        for tt = 1:length(toggle_list)
            toggle = toggle_list{tt};
            R.(toggle).F1(ii,jj) =model_acc(M.index.bic).(toggle).F1;
            R.(toggle).MCC(ii,jj) =model_acc(M.index.bic).(toggle).MCC;
            R.(toggle).TPR(ii,jj) =model_acc(M.index.bic).(toggle).TPR;
            R.(toggle).FPR(ii,jj) =model_acc(M.index.bic).(toggle).FPR;
            R.(toggle).ACC(ii,jj) =model_acc(M.index.bic).(toggle).ACC;
        end
        
    end
end
% save([resultpath,'formulation_S_result'],'R')