%% AAL116 ATLAS DESCRIPTION
inpath = './experiment/experiment_real_data/';
outpath = './experiment/experiment_real_data/';
fid = fopen([inpath,'aal116NodeNames.txt']);
data = textscan(fid,'%s%s%s');
fclose(fid);
AAL_116.name = data{1};
% According to [1], the Default mode network indices in AAL are

DMN_list = {'Frontal_Sup_Medial_L', ...
    'Frontal_Sup_Medial_R', ...
    'Cingulum_Ant_L', ...
    'Cingulum_Ant_R', ...
    'Cingulum_Post_L', ...
    'Angular_L', ...
    'Angular_R', ...
    'Precuneus_L', ...
    'Precuneus_R'};
for ii=1:length(DMN_list)
    Index = find(contains(AAL_116.name,DMN_list{ii}));
    AAL_116.DMN(ii) = Index;
end

% save([outpath,'AAL_116'],'AAL_116')

%[1] Oliver, I.; Hlinka, J.; Kopal, J.; Davidsen, J. Quantifying the Variability in Resting-State Networks. Entropy 2019, 21, 882.

