%% Consensus Result on Real data
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0);
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0);
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
y_total(:,:,1) = reshape(y_TDC,[116,size(y_TDC,2)*size(y_TDC,3)]);
y_total(:,:,2) = reshape(y_ADHD_C,[116,size(y_ADHD_C,2)*size(y_ADHD_C,3)]);
% map
AAL_116.name = {'PreCG_L','PreCG_R','SFGdor_L','SFGdor_R','ORBsup_L','ORBsup_R','MFG_L','MFG_R','ORBmid_L','ORBmid_R','IFGoperc_L','IFGoperc_R','IFGtriang_L','IFGtriang_R','ORBinf_L','ORBinf_R','ROL_L','ROL_R','SMA_L','SMA_R','OLF_L','OLF_R','SFGmed_L','SFGmed_R','ORBsupmed_L','ORBsupmed_R','REC_L','REC_R','INS_L','INS_R','ACG_L','ACG_R','MCG_L','MCG_R','PCG_L','PCG_R','HIP_L','HIP_R','PHG_L','PHG_R','AMYG_L','AMYG_R','CAL_L','CAL_R','CUN_L','CUN_R','LING_L','LING_R','SOG_L','SOG_R','MOG_L','MOG_R','IOG_L','IOG_R','FFG_L','FFG_R','PoCG_L','PoCG_R','SPG_L','SPG_R','IPG_L','IPG_R','SMG_L','SMG_R','ANG_L','ANG_R','PCUN_L','PCUN_R','PCL_L','PCL_R','CAU_L','CAU_R','PUT_L','PUT_R','PAL_L','PAL_R','THA_L','THA_R','HES_L','HES_R','STG_L','STG_R','TPOsup_L','TPOsup_R','MTG_L','MTG_R','TPOmid_L','TPOmid_R','ITG_L','ITG_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R','Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R','Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R','Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R','Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R','Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8','Vermis_9','Vermis_10'};
AAL_116.name_full = {'Precentral gyrus_L','Precentral gyrus_R','Superior frontal gyrus (dorsolateral)_L','Superior frontal gyrus (dorsolateral)_R','Superior frontal gyrus (orbital)_L','Superior frontal gyrus (orbital)_R','Middle frontal gyrus_L','Middle frontal gyrus_R','Middle frontal gyrus (orbital)_L','Middle frontal gyrus (orbital)_R','Inferior frontal gyrus (opercular)_L','Inferior frontal gyrus (opercular)_R','Inferior frontal gyrus (triangular)_L','Inferior frontal gyrus (triangular)_R','Inferior frontal gyrus (orbital)_L','Inferior frontal gyrus (orbital)_R','Rolandic operculum_L','Rolandic operculum_R','Supplementary motor area_L','Supplementary motor area_R','Olfactroy cortex_L','Olfactroy cortex_R','Superior frontal gyrus (medial)_L','Superior frontal gyrus (medial)_R','Superior frontal gyrus (medial orbital)_L','Superior frontal gyrus (medial orbital)_R','Rectus gyrus_L','Rectus gyrus_R','Insula_L','Insula_R','Anterior cingulate gyrus_L','Anterior cingulate gyrus_R','Median cingulate gyrus_L','Median cingulate gyrus_R','Posterior cingulate gyrus_L','Posterior cingulate gyrus_R','Hippocampus_L','Hippocampus_R','Parahippocampal gyrus_L','Parahippocampal gyrus_R','Amygdala_L','Amygdala_R','Calcarine cortex_L','Calcarine cortex_R','Cuneus_L','Cuneus_R','Lingual gyrus_L','Lingual gyrus_R','Superior occipital gyrus_L','Superior occipital gyrus_R','Middle occipital gyrus_L','Middle occipital gyrus_R','Inferior occipital gyrus_L','Inferior occipital gyrus_R','Fusiform gyrus_L','Fusiform gyrus_R','Postcentral gyrus_L','Postcentral gyrus_R','Superior parietal gyrus_L','Superior parietal gyrus_R','Inferior parietal gyrus_L','Inferior parietal gyrus_R','Supramarginal gyrus_L','Supramarginal gyrus_R','Angular gyrus_L','Angular gyrus_R','Precuneus_L','Precuneus_R','Paracentral lobule_L','Paracentral lobule_R','Caudate_L','Caudate_R','Putamen_L','Putamen_R','Pallidum_L','Pallidum_R','Thalamus_L','Thalamus_R','Heschl gyrus_L','Heschl gyrus_R','Superior temporal gyrus_L','Superior temporal gyrus_R','Temporal pole (superior)_L','Temporal pole (superior)_R','Middle temporal gyrus_L','Middle temporal gyrus_R','Temporal pole (middle)_L','Temporal pole (middle)_R','Inferior temporal gyrus_L','Inferior temporal gyrus_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R','Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R','Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R','Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R','Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R','Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8','Vermis_9','Vermis_10'};
%% K=2
% formulation D
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected')
M = fix_loglikelihood(M,y_total);
M = augment_score(M,size(y_total,2),'LLH_hetero');
tmp = [M.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];
eBIC = [tmp.eBIC];
[~,I] = sort(eBIC);

for ii=1:3
result.TDC_index.D2K{ii} = M.model(I(ii)).ind{1}{1};
result.ADHD_index.D2K{ii} = M.model(I(ii)).ind{1}{2};
end
clear M

load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected')
M = fix_loglikelihood(M,y_total);
M = augment_score(M,size(y_total,2),'LLH_hetero');
tmp = [M.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];
eBIC = [tmp.eBIC];
[~,I] = sort(eBIC);

for ii=1:3
result.TDC_index.S2K{ii} = M.model(I(ii)).ind{1}{1};
result.ADHD_index.S2K{ii} = M.model(I(ii)).ind{1}{2};
end
clear M
%% K=18
clear I
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered')
M.TDC = fix_loglikelihood(M.TDC,y_TDC);
M.TDC = augment_score(M.TDC,size(y_TDC,2),'LLH_hetero');
tmp = [M.TDC.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];
eBIC = [tmp.eBIC];
[~,I.TDC] = sort(eBIC);

M.ADHD_C = fix_loglikelihood(M.ADHD_C,y_ADHD_C);
M.ADHD_C = augment_score(M.ADHD_C,size(y_ADHD_C,2),'LLH_hetero');
tmp = [M.ADHD_C.model]; tmp = [tmp.stat]; tmp = [tmp.model_selection_score];
eBIC = [tmp.eBIC];
[~,I.ADHD_C] = sort(eBIC);

for ii=1:3
result.TDC_index.C18K{ii} = M.TDC.model(I.TDC(ii)).ind_common{1};
result.ADHD_index.C18K{ii} = M.ADHD_C.model(I.ADHD_C(ii)).ind_common{1};
end
clear M

save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\summary_real_DS2K_C18K_timecorrected','result')
%% Extract common along setting axis
clear
clc
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\summary_real_DS2K_C18K_timecorrected')

for jj=1:3
    TDC_common_index{jj} = setdiff(1:1:116^2,1:116+1:116^2);
    ADHD_common_index{jj} = setdiff(1:1:116^2,1:116+1:116^2);
% namelist = {'D2K','C18K','D18K'};
namelist = {'D2K','S2K','C18K'};
for ii=1:3
    TDC_common_index{jj} = intersect(TDC_common_index{jj},result.TDC_index.(namelist{ii}){jj});
    ADHD_common_index{jj} = intersect(ADHD_common_index{jj},result.ADHD_index.(namelist{ii}){jj});
end

extra_link{jj} = setdiff(ADHD_common_index{jj},TDC_common_index{jj});
missing_link{jj} = setdiff(TDC_common_index{jj},ADHD_common_index{jj});

missing_grid{jj} = zeros(116);
extra_grid{jj} = zeros(116);
missing_grid{jj}(missing_link{jj})=1;
extra_grid{jj}(extra_link{jj})=1;

AAL_grid{jj} = zeros(116);
AAL_grid{jj}(missing_link{jj}) = -1;
AAL_grid{jj}(extra_link{jj}) = 1;
end

%% for {'D2K','S2K','C18K'}
clf
close all

[cause.extra.links,cause.extra.ind] = sort(sum(extra_grid,1),'descend');
[effect.extra.links,effect.extra.ind] = sort(sum(extra_grid,2),'descend');

[cause.missing.links,cause.missing.ind] = sort(sum(missing_grid,1),'descend');
[effect.missing.links,effect.missing.ind] = sort(sum(missing_grid,2),'descend');

subtype_name = {'extra','missing'};
clc
RANK_NUMBER = 6;
for subtype=1:2
    AAL_index = cause.(subtype_name{subtype}).ind(1:RANK_NUMBER)';
    roi_list = {AAL_116.name_full{AAL_index}}';
    no_of_links = cause.extra.links(1:RANK_NUMBER)';
    
    cause_table = table(AAL_index,roi_list,no_of_links);
    fprintf(['Cause table:',subtype_name{subtype},'\n'])
    disp(cause_table)
    
    AAL_index = effect.(subtype_name{subtype}).ind(1:RANK_NUMBER);
    roi_list = {AAL_116.name_full{AAL_index}}';
    no_of_links = effect.extra.links(1:RANK_NUMBER);
    effect_table = table(AAL_index,roi_list,no_of_links);
    fprintf(['Effect table:',subtype_name{subtype},'\n'])
    disp(effect_table)
end


%% Check outward connectivity
clc
tmp = zeros(116);
tmp(TDC_common_index) = 1;
selected_index = 86;
fprintf(['(TDC) Starting node: ',AAL_116.name_full{selected_index},' \n'])
AAL_index = find(tmp(:,selected_index));
ending_node = {AAL_116.name_full{AAL_index}}';
pt = table(AAL_index,ending_node);
disp(pt)

tmp = zeros(116);
tmp(ADHD_common_index) = 1;
fprintf(['(ADHD) Starting node: ',AAL_116.name_full{selected_index},' \n'])
AAL_index = find(tmp(:,selected_index));
ending_node = {AAL_116.name_full{AAL_index}}';
pt = table(AAL_index,ending_node);
disp(pt)
%% Check inward connectivity
clc
tmp = zeros(116);
tmp(TDC_common_index) = 1;
G.TDC = digraph(tmp');
selected_index = 47; % 68<-50<-47
fprintf(['(TDC) Ending node: ',AAL_116.name_full{selected_index},' \n'])
AAL_index = find(tmp(selected_index,:))';
ending_node = {AAL_116.name_full{AAL_index}}';
pt = table(AAL_index,ending_node);
disp(pt)

tmp = zeros(116);
tmp(ADHD_common_index) = 1;
G.ADHD = digraph(tmp');
fprintf(['(ADHD) Ending node: ',AAL_116.name_full{selected_index},' \n'])
AAL_index = find(tmp(selected_index,:))';
ending_node = {AAL_116.name_full{AAL_index}}';
pt = table(AAL_index,ending_node);
disp(pt)
%% causality path
clc

tmp = zeros(116);
tmp(TDC_common_index) = 1;
G.TDC = digraph(tmp'); % transpose because our matrix is the transpose of adjacency matrix

tmp = zeros(116);
tmp(ADHD_common_index) = 1;
G.ADHD = digraph(tmp');


path_list = [74 68; ... %       Putamen_R to Precuneus_R: missing
             68 74; ... %       Precuneus_R to Putamen_R: missing
                        % ======= Literature =======
             60 9 ; ... % Superior parietal gyrus_R to Middle frontal gyrus (orbital)_L: path changed
             9 60 ; ... %       Superior parietal gyrus_R to Middle frontal gyrus (orbital)_L: extra
             18 93; ... % Rolandic operculum_R to Cerebelum_Crus2_L: no path
             93 18; ... %       Cerebelum_Crus2_L to Rolandic operculum_R: extra
             18 97; ... % Rolandic operculum_R to Cerebelum_4_5_L: path changed
             97 18; ... % Cerebelum_4_5_L to Rolandic operculum_R: path changed
             30 47; ... % Insula_R to Lingual gyrus_L: path changed
             47 30; ... % Lingual gyrus_L to Insula_R: path changed
             46 55; ... % Cuneus_R to Fusiform gyrus_L: path changed
             55 46; ... % Fusiform gyrus_L to Cuneus_R: path changed
             49 93; ... % Superior occipital gyrus_L to Cerebelum_Crus2_L: no path
             93 49; ... %       Cerebelum_Crus2_L to Superior occipital gyrus_L: extra
             61 63; ... % Inferior parietal gyrus_L to Supramarginal gyrus_L: path changed
             63 61; ... %       Supramarginal gyrus_L to Inferior parietal gyrus_L: missing
             67 89; ... % Precuneus_L to Inferior temporal gyrus_L: no path
             89 67; ... % Inferior temporal gyrus_L to Precuneus_L: path change
             85 115; ... % Middle temporal gyrus_L to Vermis_9: path change
             115 85; ... % Vermis_9 to Middle temporal gyrus_L : no path
             67 113; ... % Precuneus_L to Vermis_7: path change
             113 67];    % Vermis_7 to Precuneus_L: path change
% path_list = [25,11; ...
%              25,12; ...
%              25,13; ...
%              25,14; ...
%              25,15; ...
%              25,16; ...
%              25,19; ...
%              25,20];
%          
% path_list = [109,86; ...
%              110,86; ...
%              111,86; ...
%              112,86; ...
%              113,86; ...
%              114,86; ...
%              115,86; ...
%              116,86];
%          
%  path_list = [86,63];
%  
%  
%  path_list = [109,10; ...
%              110,10; ...
%              111,10; ...
%              112,10; ...
%              113,10; ...
%              114,10; ...
%              115,10; ...
%              116,10];
%          
%   path_list = [109,24; ...
%              110,24; ...
%              111,24; ...
%              112,24; ...
%              113,24; ...
%              114,24; ...
%              115,24; ...
%              116,24];
%   path_list = [31,30;32,30];
%   path_list = [55,7;55,8];
         
%  path_list = [path_list;path_list(:,2) path_list(:,1)];
for ii=1:size(path_list,1)
    fprintf('row number: %d\n',ii)
    P = shortestpath(G.TDC,path_list(ii,1),path_list(ii,2));
    fprintf(['Causal path TDC: <strong>',AAL_116.name_full{path_list(ii,1)},'</strong> to <strong>',AAL_116.name_full{path_list(ii,2)},'</strong> \n'])
    if isempty(P)
        fprintf(['No causal path from <strong>',AAL_116.name_full{path_list(ii,1)},'</strong> to <strong>',AAL_116.name_full{path_list(ii,2)},'</strong> \n'])
    else
        for jj=1:length(P)-1
            fprintf(['<strong> ',AAL_116.name_full{P(jj)},'</strong> ','->'])
        end
        fprintf(['<strong> ',AAL_116.name_full{P(jj+1)},'</strong> ','\n'])
    end
    
    
    P = shortestpath(G.ADHD,path_list(ii,1),path_list(ii,2));
    fprintf(['Causal path ADHD: <strong>',AAL_116.name_full{path_list(ii,1)},'</strong> to <strong>',AAL_116.name_full{path_list(ii,2)},'</strong> \n'])
    if isempty(P)
        fprintf(['No causal path from <strong>',AAL_116.name_full{path_list(ii,1)},'</strong> to <strong>',AAL_116.name_full{path_list(ii,2)},'</strong> \n'])
    else
        for jj=1:length(P)-1
            fprintf(['<strong> ',AAL_116.name_full{P(jj)},'</strong> ','->'])
        end
        fprintf(['<strong> ',AAL_116.name_full{P(jj+1)},'</strong> ','\n'])
        
    end
    fprintf('\n')
end



%%
figure(1);
% atlas_index = [74,84,61,68]; % shortest path from putamen_R -> Precuneus_R
% atlas_index = [68,12,77,47,71,74];

% atlas_index = [60,27,9]; % ADHD
% atlas_index = [9,87,60];

atlas_index = [60,35,21,16,9]; % TDC
M.TDC.GC = (M.model(M.index.eBIC).GC(:,:,1)).*(1-eye(116));  
M.ADHD_C.GC = (M.model(M.index.eBIC).GC(:,:,2)).*(1-eye(116));  
nametag = {'TDC','ADHD_C'};
titletag = {'TDC','ADHD'};

for ii=1:2

subplot(1,2,ii)
imagesc(M.(nametag{ii}).GC(atlas_index,atlas_index))
grid on
set(gca,'xtick',1:1:length(atlas_index),'ytick',1:1:length(atlas_index), ...
    'xticklabel',AAL_116.name(atlas_index),'yticklabel',AAL_116.name(atlas_index))
axis('square')
title(titletag{ii})
colormap(jet)
% caxis([0 1])
xtickangle(30)
end

%% K=18