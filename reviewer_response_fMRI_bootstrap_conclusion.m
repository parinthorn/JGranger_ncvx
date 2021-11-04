%% processing bootstrap fMRI data

%% Edges Centrality: D2K
clear
clc
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('.\experiment\ADHD\AAL_116.mat')

bootstrap_path = dir('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_D2K_seed_*_bootstrap_*_DIAGCOV.mat');
GC_bagging = zeros(116,116);
bootstrap_sample = length(bootstrap_path);

for xx=1:bootstrap_sample
    fprintf("bootstrap sample No. %d\n", xx)
    
load([inpath,bootstrap_path(xx).name])



GC.total = M.model(M.index.eBIC).GC;
GC_bagging = GC_bagging + GC.total/bootstrap_sample;

bootstrap_inference(xx).GC = GC.total;

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
bootstrap_inference(xx).table = table(Cause_name',Effect_name',ADHD_Centrality',TDC_Centrality',Centrality_Diff',Type_Diff','VariableNames',{'Cause','Effect','ADHD_centrality','TDC_centrality','Centrality Diff(ADHD-TDC)','Type'});

end

%%
clc
for ii =1:bootstrap_sample
    fprintf("bootstrap sample No. %d\n", ii)
    disp(bootstrap_inference(ii).table)
end

%% ensemble average results

[~,R.TDC] = betweenness_centrality(sparse(GC_bagging(:,:,1)')); % transpose for adjacency
[~,R.ADHD] = betweenness_centrality(sparse(GC_bagging(:,:,2)'));
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
T_bagging = table(Cause_name',Effect_name',ADHD_Centrality',TDC_Centrality',Centrality_Diff',Type_Diff','VariableNames',{'Cause','Effect','ADHD_centrality','TDC_centrality','Centrality Diff(ADHD-TDC)','Type'});

