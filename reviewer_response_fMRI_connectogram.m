clear
clc
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('.\experiment\ADHD\AAL_116.mat')
load([inpath,'estim_D2K'])
%%
groups = 1:116;
connections = [];
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



kk=2;
for zz=1:length(I)
    jj = I(zz);
    GC_causee_index = find(actual_GC(:, jj, kk)~=0);
    GC_causee = actual_GC(GC_causee_index, jj,kk);
    GC_causal_index = jj*ones(length(GC_causee),1);
    
    connections = [connections; [GC_causal_index GC_causee_index GC_causee ones(length(GC_causee),1)]];
    
end
% connections(:,1) = randi(4,50,1);
% connections(:,2) = randi(4,50,1);
% connections(:,3) = 2*rand(50,1);
% connections(5,3) = 20;
% connections(:,4) = sign((rand(50,1) > 0.5) - 0.5);
connectionsT = array2table(connections);
chordPlot(groups,connectionsT)
%%  circular graph  [most difference]
clf
Right_Left_permutation = [116,115,114,113,112,111,110,109,107,105,103,101,99,97,95,93,91,89,87,85,83,81,79,77,75,73,71,69,67,65,63,61,59,57,55,53,51,49,47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108];
actual_GC = M.model(M.index.eBIC).GC;
% actual_GC = actual_GC(Right_Left_permutation,Right_Left_permutation,:);
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

% E(abs(E_test)>THRESH)=1;
% E(abs(E_test)<=THRESH)=0;


selected_index = find(E);
[~,I] = sort(abs(E_in(selected_index)),'descend');
GC_TDC = zeros(116,116);
GC_ADHD = zeros(116,116);
for ii=1:50  %length(selected_index)
    [ee,cc] = ind2sub([116,116],selected_index(I(ii)));
    Cause_name{ii} = ([AAL_116.name_full{cc}]);  % column index
    Effect_name{ii} = ([AAL_116.name_full{ee}]);  % row index
    
    GC_TDC(ee, cc) = actual_GC(ee, cc, 1);
    GC_ADHD(ee, cc) = actual_GC(ee, cc, 2);
    
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


myLabel = AAL_116.name(Right_Left_permutation);
for jj = 1:116
    newStr = strrep(myLabel{jj},'_',' ');
    myLabel{jj} = newStr;
    
end

tt = tiledlayout(1, 2);
nexttile;
circularGraph(GC_TDC(Right_Left_permutation, Right_Left_permutation)','Label',myLabel);

nexttile;
circularGraph(GC_ADHD(Right_Left_permutation, Right_Left_permutation)','Label',myLabel);


%%  circular graph  [High GC]
clear
clc
clf
close all
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('.\experiment\ADHD\AAL_116.mat')


filename_list = {'estim_C18K', 'estim_D2K', 'estim_F2K'};
for file_index = 1:length(filename_list)
    clf
close all
    file_name  = filename_list{file_index};
    load([inpath,file_name])
Right_Left_permutation = [116,115,114,113,112,111,110,109,107,105,103,101,99,97,95,93,91,89,87,85,83,81,79,77,75,73,71,69,67,65,63,61,59,57,55,53,51,49,47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108];

if strcmp(file_name, 'estim_C18K')
    GC.total(:,:,1) = mean(M.TDC.model(M.TDC.index.eBIC).GC,3);
    GC.total(:,:,2) = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3);
else
    GC.total = M.model(M.index.eBIC).GC;
end
% actual_GC = actual_GC(Right_Left_permutation,Right_Left_permutation,:);

% E(abs(E_test)>THRESH)=1;
% E(abs(E_test)<=THRESH)=0;

GC_TDC = GC.total(:,:,1).*(1-eye(116));
GC_ADHD = GC.total(:,:,2).*(1-eye(116));
selected_index_TDC = find(GC_TDC);
selected_index_ADHD = find(GC_ADHD);
[~,I_TDC] = sort(abs(GC_TDC(selected_index_TDC)),'descend');
[~,I_ADHD] = sort(abs(GC_ADHD(selected_index_ADHD)),'descend');

GC_TDC_toplot = zeros(116,116);
GC_ADHD_toplot = zeros(116,116);

first_n_rank = 50;

for ii=1:first_n_rank
% for ii=1:length(selected_index_TDC)
    [ee,cc] = ind2sub([116,116],selected_index_TDC(I_TDC(ii)));
    
    
    GC_TDC_toplot(ee, cc) = GC.total(ee, cc, 1);
   
    
end

for ii=1:first_n_rank
% for ii=1:length(selected_index_ADHD)
    
    [ee,cc] = ind2sub([116,116],selected_index_ADHD(I_ADHD(ii)));
    GC_ADHD_toplot(ee, cc) = GC.total(ee, cc, 2);
    
    
end
ROI_color = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1];

Right_Left_permutation = Right_Left_permutation([89:116, 1:88]);
ROI_color = ROI_color([89:116, 1:88]);
myColorMap = zeros(116,3);

myLabel = AAL_116.name(Right_Left_permutation);

ROI_name = {'Posterior Fossa', ...
'Temporal lobe', ...
'Central Structures', ...
'Parietal Lobe', ...
'Occipital Lobe', ...
'Insula and Cingulate Gyri', ...
'Frontal Lobe'};

for jj = 1:116
    newStr = strrep(myLabel{jj},'_',' ');
    myLabel{jj} = newStr;
    scale = 0.8;
    switch ROI_color(jj)
        case 1
            xx = [1 0 0].*scale;
        case 2
            xx = [0 1 0].*scale;
        case 3
            xx = [0 0 1].*scale;
        case 4
            xx = [1 1 0].*scale;
        case 5
            xx = [1 0 1].*scale;
        case 6
            xx = [0 1 1].*scale;
        case 7
            xx = [0.5 0.6 0.8].*scale;
    end
    myColorMap(jj, :) = xx;
    
end
figurepath = './results2plot/figures/';
figure(1)
% nexttile;
% myColorMap = lines(length(GC_TDC_toplot));
h1=circularGraph(GC_TDC_toplot(Right_Left_permutation, Right_Left_permutation)','Colormap',myColorMap,'Label',myLabel);
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
delete(h1.ShowButton)
delete(h1.HideButton)
% print([figurepath,'reviewer_response_fMRI_D2K_circular_TDC'],'-painters','-dpng','-r300')
% print([figurepath,'reviewer_response_fMRI_D2K_circular_TDC'],'-painters','-depsc','-r300')
exportgraphics(gcf,[figurepath,'reviewer_response_',file_name,'_circular_TDC.png'])
exportgraphics(gcf,[figurepath,'reviewer_response_',file_name,'_circular_TDC.eps'])

% nexttile;
figure(2)
h2=circularGraph(GC_ADHD_toplot(Right_Left_permutation, Right_Left_permutation)','Colormap',myColorMap,'Label',myLabel);
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
delete(h2.ShowButton)
delete(h2.HideButton)
set(gcf, 'Position', pp);
% print([figurepath,'reviewer_response_fMRI_D2K_circular_ADHD'],'-painters','-dpng','-r300')
% print([figurepath,'reviewer_response_fMRI_D2K_circular_ADHD'],'-painters','-depsc','-r300')
exportgraphics(gcf,[figurepath,'reviewer_response_',file_name,'_circular_ADHD.png'])
exportgraphics(gcf,[figurepath,'reviewer_response_',file_name,'_circular_ADHD.eps'])
end

figure(3)
hold on
for jj=1:7
    switch jj
        case 1
            xx = [1 0 0].*scale;
        case 2
            xx = [0 1 0].*scale;
        case 3
            xx = [0 0 1].*scale;
        case 4
            xx = [1 1 0].*scale;
        case 5
            xx = [1 0 1].*scale;
        case 6
            xx = [0 1 1].*scale;
        case 7
            xx = [0.5 0.6 0.8].*scale;
    end
plot(NaN,NaN,'.','Color',xx, 'MarkerSize',69);
end
hold off
lh = legend(ROI_name);
lh.Position(1) = 0.5 - lh.Position(3)/2; 
lh.Position(2) = 0.518 - lh.Position(4)/2;
set(gca,'FontSize', 32, 'xticklabel',[], 'yticklabel',[])
% print([figurepath,'reviewer_response_fMRI_legend'],'-painters','-dpng','-r300')
% print([figurepath,'reviewer_response_fMRI_legend'],'-painters','-depsc','-r300')
exportgraphics(gcf,[figurepath,'reviewer_response_fMRI_legend.png'])
exportgraphics(gcf,[figurepath,'reviewer_response_fMRI_legend.eps'])
%% Bootstrap D2K

clear
clc
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\';
load('.\experiment\ADHD\AAL_116.mat')

bootstrap_path = dir('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_D2K_seed_*_bootstrap_*_DIAGCOV.mat');
GC.total = zeros(116,116, 2);
GC.repetition = zeros(116, 116, 2);

bootstrap_sample = length(bootstrap_path);

for xx=1:bootstrap_sample
    fprintf("bootstrap sample No. %d\n", xx)
    
load([inpath,bootstrap_path(xx).name])

GC.total = GC.total +  (M.model(M.index.eBIC).GC)/bootstrap_sample;
GC.repetition = GC.repetition + (M.model(M.index.eBIC).GC~=0)/bootstrap_sample;

end

%%
clf
close all
GC_TDC = GC.total(:,:,1).*(1-eye(116)).*(GC.repetition(:, :, 1)>0.95);
GC_ADHD = GC.total(:,:,2).*(1-eye(116)).*(GC.repetition(:, :, 2)>0.95);
% GC_TDC = GC.total(:,:,1).*(1-eye(116));
% GC_ADHD = GC.total(:,:,2).*(1-eye(116));

selected_index_TDC = find(GC_TDC);
selected_index_ADHD = find(GC_ADHD);
[~,I_TDC] = sort(abs(GC_TDC(selected_index_TDC)),'descend');
[~,I_ADHD] = sort(abs(GC_ADHD(selected_index_ADHD)),'descend');

GC_TDC_toplot = zeros(116,116);
GC_ADHD_toplot = zeros(116,116);

first_n_rank = 50;

% for ii=1:first_n_rank
% % for ii=1:length(selected_index_TDC)
%     [ee,cc] = ind2sub([116,116],selected_index_TDC(I_TDC(ii)));
%     
%     
%     GC_TDC_toplot(ee, cc) = GC.total(ee, cc, 1);
%    
%     
% end
% 
% for ii=1:first_n_rank
% % for ii=1:length(selected_index_ADHD)
%     
%     [ee,cc] = ind2sub([116,116],selected_index_ADHD(I_ADHD(ii)));
%     GC_ADHD_toplot(ee, cc) = GC.total(ee, cc, 2);
%     
%     
% end

GC_TDC_toplot = GC_TDC;
GC_ADHD_toplot = GC_ADHD;


ROI_color = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1];
Right_Left_permutation = [116,115,114,113,112,111,110,109,107,105,103,101,99,97,95,93,91,89,87,85,83,81,79,77,75,73,71,69,67,65,63,61,59,57,55,53,51,49,47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108];

Right_Left_permutation = Right_Left_permutation([89:116, 1:88]);
ROI_color = ROI_color([89:116, 1:88]);
myColorMap = zeros(116,3);

myLabel = AAL_116.name(Right_Left_permutation);

ROI_name = {'Posterior Fossa', ...
            'Temporal lobe', ...
            'Central Structures', ...
            'Parietal Lobe', ...
            'Occipital Lobe', ...
            'Insula and Cingulate Gyri', ...
            'Frontal Lobe'};

for jj = 1:116
    newStr = strrep(myLabel{jj},'_',' ');
    myLabel{jj} = newStr;
    scale = 0.8;
    switch ROI_color(jj)
        case 1
            xx = [1 0 0].*scale;
        case 2
            xx = [0 1 0].*scale;
        case 3
            xx = [0 0 1].*scale;
        case 4
            xx = [1 1 0].*scale;
        case 5
            xx = [1 0 1].*scale;
        case 6
            xx = [0 1 1].*scale;
        case 7
            xx = [0.5 0.6 0.8].*scale;
    end
    myColorMap(jj, :) = xx;
    
end
figurepath = './results2plot/figures/';
figure(1)
% nexttile;
% myColorMap = lines(length(GC_TDC_toplot));
h1=circularGraph(GC_TDC_toplot(Right_Left_permutation, Right_Left_permutation)','Colormap',myColorMap,'Label',myLabel);
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
delete(h1.ShowButton)
delete(h1.HideButton)
exportgraphics(gcf,[figurepath,'reviewer_response_estim_D2Kbootstrap_circular_TDC.png'])
exportgraphics(gcf,[figurepath,'reviewer_response_estim_D2Kbootstrap_circular_TDC.eps'])

figure(2)
h2=circularGraph(GC_ADHD_toplot(Right_Left_permutation, Right_Left_permutation)','Colormap',myColorMap,'Label',myLabel);
pp = get(0, 'Screensize');
pp(3) = pp(3)*1;
set(gcf, 'Position', pp);
delete(h2.ShowButton)
delete(h2.HideButton)
exportgraphics(gcf,[figurepath,'reviewer_response_estim_D2Kbootstrap_circular_ADHD.png'])
exportgraphics(gcf,[figurepath,'reviewer_response_estim_D2Kbootstrap_circular_ADHD.eps'])

%%
% ORBsupmedL -> ACG.L (blue)
% REC.L -> PHG.R (blue)
% ORBsupmed.R -> ACG.R (red)
clf
close all
myLabel = AAL_116.name([25, 26, 27, 31, 32, 40]);
for jj = 1:length(myLabel)
    newStr = strrep(myLabel{jj},'_',' ');
    myLabel{jj} = newStr;
    
end

myColorMap=[0 0 1; ...
            1 0 0; ...
            0 0 1; ...
            0 0 0; ...
            0 0 0; ...
            0 0 0];

GC=[0 0 0 0 0 0; ...
    0 0 0 0 0 0; ...
    0 0 0 0 0 0; ...
    1 0 0 0 0 0; ...
    0 1 0 0 0 0; ...
    0 0 1 0 0 0];
 % 25 26 27 31 32 40 
permutation_index = [1,3,2,5,4,6];
 
 h2=circularGraph(GC(permutation_index,permutation_index)','Colormap',myColorMap(permutation_index, :),'Label',myLabel(permutation_index));
h2.Node(4).Visible = 0;
h2.Node(5).Visible = 0;
h2.Node(6).Visible = 0;
% pp = get(0, 'Screensize');
% pp(3) = pp(3)*1;
% set(gcf, 'Position', pp);
delete(h2.ShowButton)
delete(h2.HideButton)

figurepath = './results2plot/figures/';

exportgraphics(gcf,[figurepath,'reviewer_response_mostdiffnode.png'])
exportgraphics(gcf,[figurepath,'reviewer_response_mostdiffnode.eps'])