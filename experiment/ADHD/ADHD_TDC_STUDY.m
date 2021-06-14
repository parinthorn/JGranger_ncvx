%% Real data experiment
% The goal of this experiment is to study the effective connectivity
% differences between ADHD and TDC by using proposed formulations.
% Setting
% - use cvx-DGN, cvx-FGN to estimate concatenation of K=18 subjects in each of ADHD and TDC data sets
% - use cvx-CGN to estimate two models from K=18 ADHD and TDC respectively
% the common part of GC matrix in each model will be considered as group
% level GC.

%% Data concatenation [Run this first]
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0); % load data in format (n,T,K)
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
y_total = cat(3,y_TDC,y_ADHD_C);
y_total = y_total-mean(y_total,2); % detrend in time axis
%% cvx-DGN estimation
parameter.varorder = 1;
parameter.formulation = 'dgn';
parameter.penalty_weight = 'LS';
parameter.GridSize = 30; % resolution of regularization grid
parameter.data_concat = 1; % set to 1 for time-series concatenation without dependency between patients
parameter.noisecov = 'diag'; % diagonal covariance assumptions
parameter.qnorm = 'cvx';
ALG_PARAMETER = gen_alg_params(parameter.qnorm, parameter.formulation);
M = jointvargc(y_total,parameter,ALG_PARAMETER); % data with dimension (n,T*K,2), K is # subjects in each TDC, ADHD
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_D2K','M') % verified for reproduce
%% cvx-FGN estimation
parameter.varorder = 1;
parameter.formulation = 'fgn';
parameter.penalty_weight = 'LS';
parameter.GridSize = 30; % resolution of regularization grid
parameter.data_concat = 1; % set to 1 for time-series concatenation without dependency between patients
parameter.noisecov = 'diag'; % diagonal covariance assumptions
parameter.qnorm = 'cvx';
ALG_PARAMETER = gen_alg_params(parameter.qnorm, parameter.formulation);
M = jointvargc(y_total,parameter,ALG_PARAMETER); % data with dimension (n,T*K,2), K is # subjects in each TDC, ADHD
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_F2K','M') % verified for reproduce
%% cvx-CGN estimation
clear M
parameter.varorder = 1;
parameter.formulation = 'cgn';
parameter.penalty_weight = 'LS';
parameter.GridSize = 30; % resolution of regularization grid
parameter.data_concat = 0;
parameter.noisecov = 'diag'; % diagonal covariance assumptions
parameter.qnorm = 'cvx';
ALG_PARAMETER = gen_alg_params(parameter.qnorm, parameter.formulation);
M.TDC = jointvargc(y_TDC,parameter,ALG_PARAMETER); % data with dimension (n,T*K,2), K is # subjects in each TDC, ADHD
M.ADHD_C = jointvargc(y_ADHD_C,parameter,ALG_PARAMETER); % data with dimension (n,T*K,2), K is # subjects in each TDC, ADHD
save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_C18K','M') % verified for reproduce
%% Convert to weighted Adjacency matrix
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
outpath = './experiment/result_to_plot/';
% D2K
load([inpath,'estim_D2K'])
AdjTDC = M.model(M.index.eBIC).GC(:,:,1)';
AdjADHD = M.model(M.index.eBIC).GC(:,:,2)';
writematrix(AdjTDC,[outpath,'AdjTDC_D2K.txt'],'Delimiter','tab')
writematrix(AdjADHD,[outpath,'AdjADHD_D2K.txt'],'Delimiter','tab')
% F2K
clear M
load([inpath,'estim_F2K'])
AdjTDC = M.model(M.index.eBIC).GC(:,:,1)';
AdjADHD = M.model(M.index.eBIC).GC(:,:,2)';
writematrix(AdjTDC,[outpath,'AdjTDC_F2K.txt'],'Delimiter','tab')
writematrix(AdjADHD,[outpath,'AdjADHD_F2K.txt'],'Delimiter','tab')
% C18K
clear M
load([inpath,'estim_C18K'])
AdjTDC = mean(M.TDC.model(M.TDC.index.eBIC).GC,3)';
AdjADHD = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3)';
writematrix(AdjTDC,[outpath,'AdjTDC_C18K.txt'],'Delimiter','tab')
writematrix(AdjADHD,[outpath,'AdjADHD_C18K.txt'],'Delimiter','tab')
%% Edges Centrality: D2K
clear
clc
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('.\experiment\experiment_real_data\AAL_116.mat')
load([inpath,'estim_D2K'])
GC.total = M.model(M.index.eBIC).GC;
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
[~,R.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)')); % transpose for adjacency
[~,R.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));
E_TDC = R.TDC';
E_ADHD = R.ADHD';

E_in = R.TDC'-R.ADHD'; % convert to GC
E=E_in;
THRESH =5*std(E_in(E_in~=0));
E(E_in~=0) = E(E_in~=0)-mean(E(E_in~=0));
E_test = E;

E(abs(E_test)>THRESH)=1;
E(abs(E_test)<=THRESH)=0;


selected_index = find(E);
[~,I] = sort(abs(E_in(selected_index)),'descend');


for ii=1:length(selected_index)
    [ee,cc] = ind2sub([116,116],selected_index(I(ii)));
    Cause_name{ii} = ([AAL_116.name_full{cc}]);
    Effect_name{ii} = ([AAL_116.name_full{ee}]);
    
    
    Centrality_Diff(ii) = -E_in(selected_index(I(ii)));
    ADHD_Centrality(ii) = E_ADHD(selected_index(I(ii)));
    TDC_Centrality(ii) = E_TDC(selected_index(I(ii)));
    
    if Centrality_Diff(ii) <0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
end
T = table(Cause_name',Effect_name',ADHD_Centrality',TDC_Centrality',Centrality_Diff',Type_Diff','VariableNames',{'Cause','Effect','ADHD_centrality','TDC_centrality','Centrality Diff(ADHD-TDC)','Type'});

%% Edges Centrality: F2K
clear
clc
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load([inpath,'estim_F2K'])
GC.total = M.model(M.index.eBIC).GC;
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
[~,R.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)')); % transpose for adjacency
[~,R.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));
E_TDC = R.TDC';
E_ADHD = R.ADHD';

E_in = R.TDC'-R.ADHD'; % convert to GC
E=E_in;
THRESH =5*std(E_in(E_in~=0));
E(E_in~=0) = E(E_in~=0)-mean(E(E_in~=0));
E_test = E;

E(abs(E_test)>THRESH)=1;
E(abs(E_test)<=THRESH)=0;


selected_index = find(E);
[~,I] = sort(abs(E_in(selected_index)),'descend');


for ii=1:length(selected_index)
    [ee,cc] = ind2sub([116,116],selected_index(I(ii)));
    Cause_name{ii} = ([AAL_116.name_full{cc}]);
    Effect_name{ii} = ([AAL_116.name_full{ee}]);
    
    
    Centrality_Diff(ii) = -E_in(selected_index(I(ii)));
    ADHD_Centrality(ii) = E_ADHD(selected_index(I(ii)));
    TDC_Centrality(ii) = E_TDC(selected_index(I(ii)));
    
    if Centrality_Diff(ii) <0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
end
T = table(Cause_name',Effect_name',ADHD_Centrality',TDC_Centrality',Centrality_Diff',Type_Diff','VariableNames',{'Cause','Effect','ADHD_centrality','TDC_centrality','Centrality Diff(ADHD-TDC)','Type'});

%% Edges Centrality: C18K
clear
clc
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load([inpath,'estim_C18K'])
GC.total(:,:,1) = mean(M.TDC.model(M.TDC.index.eBIC).GC,3);

GC.total(:,:,2) = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3);
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
[~,R.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)')); % transpose for adjacency
[~,R.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));

E_TDC = R.TDC';
E_ADHD = R.ADHD';
E_in = R.TDC'-R.ADHD'; % convert to GC
E=E_in;
THRESH =5*std(E(E~=0));
E(E_in~=0) = E(E_in~=0)-mean(E(E_in~=0));
E_test = E;

E(abs(E_test)>THRESH)=1;
E(abs(E_test)<=THRESH)=0;


selected_index = find(E);
[~,I] = sort(abs(E_in(selected_index)),'descend');


for ii=1:length(selected_index)
    [ee,cc] = ind2sub([116,116],selected_index(I(ii)));
    Cause_name{ii} = ([AAL_116.name_full{cc}]);
    Effect_name{ii} = ([AAL_116.name_full{ee}]);
    
    
    Centrality_Diff(ii) = -E_in(selected_index(I(ii)));
    ADHD_Centrality(ii) = E_ADHD(selected_index(I(ii)));
    TDC_Centrality(ii) = E_TDC(selected_index(I(ii)));
    
    if Centrality_Diff(ii) <0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
end
T = table(Cause_name',Effect_name',ADHD_Centrality',TDC_Centrality',Centrality_Diff',Type_Diff','VariableNames',{'Cause','Effect','ADHD_centrality','TDC_centrality','Centrality Diff(ADHD-TDC)','Type'});




