clear
clc
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\ALL_RESULT_formulation_D_result_regrid_K50all.mat')
dd=1;
realz = 20;

score_list = {'TPR','FPR','ACC','F1','MCC'};
for jj=1:length(score_list)
    avg.(score_list{jj}) = zeros(30,30);
end
for ii=1:realz
    tmp = [ALL_RESULT(dd,ii).model_acc.total];
    for jj=1:length(score_list)
    R(ii).(score_list{jj}) = reshape([tmp.(score_list{jj})],[30,30]);
    
    avg.(score_list{jj})= avg.(score_list{jj})+R(ii).(score_list{jj})/realz;
    end
%     [~,best_index(ii,1)] = max((1-[tmp.FPR]).*[tmp.TPR]);
    [~,best_index(ii,1)] = max([tmp.F1]);
end



