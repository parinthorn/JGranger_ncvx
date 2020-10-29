%% this script evaluate performance of ....
clear
clc
close all
modelpath = './data_compare/';
inpath = './experiment/Experiment_compare/skrip_code_2/data_R_formulationD/';
outpath = './experiment/Experiment_compare/skrip_code_2/data_R_formulationD/';
type = 2; %D type
cd = 3; %common density set to percent(cd); percent=[1%, 5%, 10%, 20%]
T = 100;
p = 1;
K = 5;
n = 20; % time-series channels
[P,~] = offdiagJSS(n,p,K);
load([modelpath,'model_K',int2str(K),'_p1']) % struct E
[~,~,dd,m] = size(E);
m=20;
GridSize = 30;
mname = {'1','5'};

for ii=1:dd
    score(ii).total.TPR = 0;
score(ii).total.FPR = 0;
score(ii).total.ACC = 0;
score(ii).total.F1 = 0;
score(ii).bias = 0;

    score(ii).common.TPR = 0;
score(ii).common.FPR = 0;
score(ii).common.ACC = 0;
score(ii).common.F1 = 0;

    score(ii).differential.TPR = 0;
score(ii).differential.FPR = 0;
score(ii).differential.ACC = 0;
score(ii).differential.F1 = 0;
    for jj=1:m
        % performance eval
        model = E{type,cd,ii,jj};
        D = readtable([inpath,'K5_Final_Est_',mname{ii},'percent_',int2str(jj),'.csv']);
        D = table2array(D);
        D = D(:,2:end); % D = size (n,n*K), p=1
        A = reshape(D,[n,n,K]);
        ind_nz = cell(K,1);
        for kk=1:K
            ind_nz{kk} = setdiff(find(A(:,:,kk)),1:n+1:n^2);
        end
        [ind_common,ind_differential] = split_common_diff(ind_nz,[n,p,K]); % find common and differential off-diagonal nonzero index
        common_confusion = compare_sparsity(model.ind_common,ind_common{1},n,K,'single_common');
        common_score = performance_score(common_confusion);
        differential_confusion = compare_sparsity(model.ind_differential,ind_differential{1},n,K,'single_differential');
        differential_score = performance_score(differential_confusion);
        total_confusion = compare_sparsity(model.ind,ind_nz,n,K,'single_differential');
        total_score = performance_score(total_confusion);
        total_score.bias = sqrt(sum((A-model.A).^2,'all')/sum(model.A.^2,'all'));
        score(ii).total.TPR = score(ii).total.TPR+total_score.TPR/m;
        score(ii).total.FPR = score(ii).total.FPR+total_score.FPR/m;
        score(ii).total.ACC = score(ii).total.ACC+total_score.ACC/m;
        score(ii).total.F1 = score(ii).total.F1+total_score.F1/m;
        score(ii).bias = score(ii).bias+total_score.bias/m;
        
        score(ii).common.TPR = score(ii).common.TPR+common_score.TPR/m;
        score(ii).common.FPR = score(ii).common.FPR+common_score.FPR/m;
        score(ii).common.ACC = score(ii).common.ACC+common_score.ACC/m;
        score(ii).common.F1 = score(ii).common.F1+common_score.F1/m;
        
        score(ii).differential.TPR = score(ii).differential.TPR+differential_score.TPR/m;
        score(ii).differential.FPR = score(ii).differential.FPR+differential_score.FPR/m;
        score(ii).differential.ACC = score(ii).differential.ACC+differential_score.ACC/m;
        score(ii).differential.F1 = score(ii).differential.F1+differential_score.F1/m;
    end
end

% save([outpath,'skrip_formulationD_accuracy'],'score')






