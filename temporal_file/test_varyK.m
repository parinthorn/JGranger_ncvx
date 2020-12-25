clear
clc
inpath = './data_compare/';

% resultpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
% resultpath = 'D:/parinthorn_thesis/formulation_D_result/';
resultpath = 'C:/Users/CU_EE_LAB408/Desktop/tmp/'; 
% load('G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/best_index_K50.mat')
mname = {'1','5'};
dd = length(mname);
dd=1;
% dd=2;
realization = 5;
load([inpath,'model_K50_p1'])
K_list = [5,15,25,35,50];
for ii=1:dd
    k_count = 0;
    for K=K_list
        k_count = k_count+1;
        for jj=1:realization
            fprintf('(%d,%d)\n',ii,jj)
            original_model = E{2,3,ii,jj};
            GTmodel = extract_group(original_model,1:K);
%             fname = [resultpath,'varyK_result_formulationD_',mname{ii},'percent_lag1_K',int2str(K),'_',int2str(jj)];
            fname = [resultpath,'cvx_varyK_result_formulationD_',mname{ii},'percent_lag1_K',int2str(K),'_',int2str(jj)];

            load(fname)
            model_acc = performance_eval(M,GTmodel);
            toggle_list = {'total','common','differential'};
%             M.index.bic=best_index(jj);
            ALL_RESULT(ii,k_count,jj).model_acc = model_acc;
            ALL_RESULT(ii,k_count,jj).index = M.index;
            for tt = 1:length(toggle_list)
                toggle = toggle_list{tt};
                R.(toggle).F1(ii,k_count,jj) =model_acc(M.index.bic).(toggle).F1;
                R.(toggle).MCC(ii,k_count,jj) =model_acc(M.index.bic).(toggle).MCC;
                R.(toggle).TPR(ii,k_count,jj) =model_acc(M.index.bic).(toggle).TPR;
                R.(toggle).FPR(ii,k_count,jj) =model_acc(M.index.bic).(toggle).FPR;
                R.(toggle).ACC(ii,k_count,jj) =model_acc(M.index.bic).(toggle).ACC;
                
            end
        end
        
    end
end
% 
% save([resultpath,'cvx_varyK_formulation_D_result_acc'],'R')
% save([resultpath,'cvx_varyK_formulation_D_result_ALL_RESULT'],'ALL_RESULT')
