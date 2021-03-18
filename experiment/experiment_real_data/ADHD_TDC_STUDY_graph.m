%% D2K
clear
clc
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
G.TDC = digraph(GC.total(:,:,1)',AAL_116.name);
G.ADHD = digraph(GC.total(:,:,2)',AAL_116.name);
% 
C1 = (centrality(G.TDC,'betweenness','Cost',G.TDC.Edges.Weight));
C2 = (centrality(G.ADHD,'betweenness','Cost',G.ADHD.Edges.Weight));

C_in = (C1-C2);
C=C_in;
THRESH =2*std(C_in);

C(abs(C_in-mean(C_in))>THRESH)=1;
C(abs(C_in-mean(C_in))<=THRESH)=0;
% 
% figure(1)
% subplot(1,2,1)
% p1 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
% p1.NodeCData = C;
% p1.MarkerSize = 7;
% colormap(1-gray)
% 
% subplot(1,2,2)
% p2 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
% p2.NodeCData = C;
% p2.MarkerSize = 7;
% colormap(1-gray)

Network_list = {'DMN','DMNa','DMNv','SM','Visual','SN','Cerebellar','Vermis'};
index.DMN = [59,60,61,62,85,86];
index.DMNa = [29,30,31,32,87,88];
index.DMNv = [35,36,37,38,39,40,55,56,65,66,67,68];
index.SM = [1,2,7,8,19,20,57,58,63,64,69,70];
index.Visual = [43,44,45,46,47,48,49,50,51,52,53,54];
index.SN = [7,8,9,10, 29,30, 23,24, 19,20];
index.Cerebellar = [91:108];
index.Vermis = [109:116];
for nn=1:116
    clear Network
    cnt = 0;
    Network = [];
    for ss = 1:8
        if any(nn == index.(Network_list{ss}),'all')
            cnt=cnt+1;
            if cnt>1
           Network = [Network,',',Network_list{ss}];
            else
                Network = [Network_list{ss}];
            end
        end
    end
    if cnt==0
        AAL_116.Network{nn} = 'unassigned';
    else
        AAL_116.Network{nn} = Network;
    end
end

selected_index = find(C);
[~,I] = sort(abs(C_in(selected_index)),'descend');
% fprintf('Difference nodes are,\n')
for ii=1:length(selected_index)
    ROI_name{ii} = ([AAL_116.name_full{selected_index(I(ii))}]);
    Centrality_Diff(ii) = C_in(selected_index(I(ii)));
    
    Network_type{ii} = AAL_116.Network(selected_index(I(ii)));
    
    if C_in(selected_index(I(ii))) >0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
    
end
T = table(ROI_name',Centrality_Diff',Type_Diff','VariableNames',{'ROI','Centrality Diff(TDC-ADHD)','Type'});
disp(T)
%% S2K
clear
clc
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
G.TDC = digraph(GC.total(:,:,1)',AAL_116.name);
G.ADHD = digraph(GC.total(:,:,2)',AAL_116.name);
C1 = (centrality(G.TDC,'betweenness','Cost',G.TDC.Edges.Weight));
C2 = (centrality(G.ADHD,'betweenness','Cost',G.ADHD.Edges.Weight));

C_in = (C1-C2);
C=C_in;
THRESH =2*std(C_in);

C(abs(C_in-mean(C_in))>THRESH)=1;
C(abs(C_in-mean(C_in))<=THRESH)=0;

% figure(1)
% subplot(1,2,1)
% p1 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
% p1.NodeCData = C;
% p1.MarkerSize = 7;
% colormap(1-gray)
% 
% subplot(1,2,2)
% p2 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
% p2.NodeCData = C;
% p2.MarkerSize = 7;
% colormap(1-gray)
Network_list = {'DMN','DMNa','DMNv','SM','Visual','SN','Cerebellar','Vermis'};
index.DMN = [59,60,61,62,85,86];
index.DMNa = [29,30,31,32,87,88];
index.DMNv = [35,36,37,38,39,40,55,56,65,66,67,68];
index.SM = [1,2,7,8,19,20,57,58,63,64,69,70];
index.Visual = [43,44,45,46,47,48,49,50,51,52,53,54];
index.SN = [7,8,9,10, 29,30, 23,24, 19,20];
index.Cerebellar = [91:108];
index.Vermis = [109:116];
for nn=1:116
    clear Network
    cnt = 0;
    Network = [];
    for ss = 1:8
        if any(nn == index.(Network_list{ss}),'all')
            cnt=cnt+1;
            if cnt>1
           Network = [Network,',',Network_list{ss}];
            else
                Network = [Network_list{ss}];
            end
        end
    end
    if cnt==0
        AAL_116.Network{nn} = 'unassigned';
    else
        AAL_116.Network{nn} = Network;
    end
end




selected_index = find(C);
[~,I] = sort(abs(C_in(selected_index)),'descend');
% fprintf('Difference nodes are,\n')
for ii=1:length(selected_index)
    ROI_name{ii} = ([AAL_116.name_full{selected_index(I(ii))}]);
    Centrality_Diff(ii) = C_in(selected_index(I(ii)));
    Network_type{ii} = AAL_116.Network(selected_index(I(ii)));
    if C_in(selected_index(I(ii))) >0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
end
T = table(ROI_name',Centrality_Diff',Type_Diff','VariableNames',{'ROI','Centrality Diff(TDC-ADHD)','Type'});
disp(T)

%% C18K
clear
clc
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered_LLHcorrection','M')
GC.total(:,:,1) = mean(M.TDC.model(M.TDC.index.eBIC).GC,3);
GC.total(:,:,2) = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3);
GC.total(GC.total>0) = 1./GC.total(GC.total>0);

G.TDC = digraph(GC.total(:,:,1)',AAL_116.name);
G.ADHD = digraph(GC.total(:,:,2)',AAL_116.name);
C1 = (centrality(G.TDC,'betweenness','Cost',G.TDC.Edges.Weight));
C2 = (centrality(G.ADHD,'betweenness','Cost',G.ADHD.Edges.Weight));

C_in = (C1-C2);
C=C_in;
THRESH =2*std(C_in);

C(abs(C_in-mean(C_in))>THRESH)=1;
C(abs(C_in-mean(C_in))<=THRESH)=0;

% figure(1)
% subplot(1,2,1)
% p1 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
% p1.NodeCData = C;
% p1.MarkerSize = 7;
% colormap(1-gray)
% 
% subplot(1,2,2)
% p2 = plot(G.TDC,'XData',AAL_116.MNI_coord(:,1),'YData',AAL_116.MNI_coord(:,2),'ZData',AAL_116.MNI_coord(:,3));
% p2.NodeCData = C;
% p2.MarkerSize = 7;
% colormap(1-gray)
selected_index = find(C);
[~,I] = sort(abs(C_in(selected_index)),'descend');
% fprintf('Difference nodes are,\n')
for ii=1:length(selected_index)
    ROI_name{ii} = ([AAL_116.name_full{selected_index(I(ii))}]);
    Centrality_Diff(ii) = C_in(selected_index(I(ii)));
    if C_in(selected_index(I(ii))) >0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
end
T = table(ROI_name',Centrality_Diff',Type_Diff','VariableNames',{'ROI','Centrality Diff(TDC-ADHD)','Type'});
disp(T)
%% Edges Centrality
clear
clc
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
[~,R.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)')); % transpose for adjacency
[~,R.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));

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
    
    
    Centrality_Diff(ii) = E_in(selected_index(I(ii)));
    
    if Centrality_Diff(ii) >0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
end
T = table(Cause_name',Effect_name',Centrality_Diff',Type_Diff','VariableNames',{'Cause','Effect','Centrality Diff(TDC-ADHD)','Type'});

%% Edges Centrality
clear
clc
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
[~,R.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)')); % transpose for adjacency
[~,R.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));

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
    
    
    Centrality_Diff(ii) = E_in(selected_index(I(ii)));
    
    if Centrality_Diff(ii) >0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
end
T = table(Cause_name',Effect_name',Centrality_Diff',Type_Diff','VariableNames',{'Cause','Effect','Centrality Diff(TDC-ADHD)','Type'});

%% Edges Centrality
clear
clc
load('E:\JGranger_ncvx\experiment\experiment_real_data\AAL_116.mat')
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered_LLHcorrection')
GC.total(:,:,1) = mean(M.TDC.model(M.TDC.index.eBIC).GC,3);

GC.total(:,:,2) = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3);
GC.total(GC.total>0) = 1./GC.total(GC.total>0);
[~,R.TDC] = betweenness_centrality(sparse(GC.total(:,:,1)')); % transpose for adjacency
[~,R.ADHD] = betweenness_centrality(sparse(GC.total(:,:,2)'));

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
    
    
    Centrality_Diff(ii) = E_in(selected_index(I(ii)));
    
    if Centrality_Diff(ii) >0
        toggle = 'TDC>ADHD(Missing)';
    else
        toggle = 'TDC<ADHD(Extra)';
    end
    Type_Diff{ii} = toggle;
    
end
T = table(Cause_name',Effect_name',Centrality_Diff',Type_Diff','VariableNames',{'Cause','Effect','Centrality Diff(TDC-ADHD)','Type'});
