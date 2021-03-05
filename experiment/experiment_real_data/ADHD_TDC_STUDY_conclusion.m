%% sorting important edges
clear
clc
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\summary_real_DS2K_C18K_timecorrected_LLHcorrected')
namelist = {'D2K','S2K','C18K'};
for jj=1:3
    for ii=1:3
        link.extra{ii,jj} = setdiff(result.ADHD_index.(namelist{ii}){jj},result.TDC_index.(namelist{ii}){jj});
        link.missing{ii,jj} = setdiff(result.TDC_index.(namelist{ii}){jj},result.ADHD_index.(namelist{ii}){jj});
        grid.missing{ii,jj} = zeros(116);
        grid.extra{ii,jj} = zeros(116);
        grid.missing{ii,jj}(link.missing{ii,jj})=1;
        grid.extra{ii,jj}(link.extra{ii,jj})=1;
        AAL_grid{ii,jj} = zeros(116);
        AAL_grid{ii,jj}(link.missing{ii,jj}) = -1;
        AAL_grid{ii,jj}(link.extra{ii,jj}) = 1;
    end
end
tt = tiledlayout(3,3);
tt.Padding = 'compact';
tt.TileSpacing = 'compact';
for ii=1:3
    for jj=1:3
        nexttile;
        imagesc(AAL_grid{ii,jj})
    end
end
%% print all combination
RANK_NUMBER = 6;
AAL_116.name = {'PreCG_L','PreCG_R','SFGdor_L','SFGdor_R','ORBsup_L','ORBsup_R','MFG_L','MFG_R','ORBmid_L','ORBmid_R','IFGoperc_L','IFGoperc_R','IFGtriang_L','IFGtriang_R','ORBinf_L','ORBinf_R','ROL_L','ROL_R','SMA_L','SMA_R','OLF_L','OLF_R','SFGmed_L','SFGmed_R','ORBsupmed_L','ORBsupmed_R','REC_L','REC_R','INS_L','INS_R','ACG_L','ACG_R','MCG_L','MCG_R','PCG_L','PCG_R','HIP_L','HIP_R','PHG_L','PHG_R','AMYG_L','AMYG_R','CAL_L','CAL_R','CUN_L','CUN_R','LING_L','LING_R','SOG_L','SOG_R','MOG_L','MOG_R','IOG_L','IOG_R','FFG_L','FFG_R','PoCG_L','PoCG_R','SPG_L','SPG_R','IPG_L','IPG_R','SMG_L','SMG_R','ANG_L','ANG_R','PCUN_L','PCUN_R','PCL_L','PCL_R','CAU_L','CAU_R','PUT_L','PUT_R','PAL_L','PAL_R','THA_L','THA_R','HES_L','HES_R','STG_L','STG_R','TPOsup_L','TPOsup_R','MTG_L','MTG_R','TPOmid_L','TPOmid_R','ITG_L','ITG_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R','Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R','Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R','Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R','Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R','Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8','Vermis_9','Vermis_10'};
AAL_116.name_full = {'Precentral gyrus_L','Precentral gyrus_R','Superior frontal gyrus (dorsolateral)_L','Superior frontal gyrus (dorsolateral)_R','Superior frontal gyrus (orbital)_L','Superior frontal gyrus (orbital)_R','Middle frontal gyrus_L','Middle frontal gyrus_R','Middle frontal gyrus (orbital)_L','Middle frontal gyrus (orbital)_R','Inferior frontal gyrus (opercular)_L','Inferior frontal gyrus (opercular)_R','Inferior frontal gyrus (triangular)_L','Inferior frontal gyrus (triangular)_R','Inferior frontal gyrus (orbital)_L','Inferior frontal gyrus (orbital)_R','Rolandic operculum_L','Rolandic operculum_R','Supplementary motor area_L','Supplementary motor area_R','Olfactroy cortex_L','Olfactroy cortex_R','Superior frontal gyrus (medial)_L','Superior frontal gyrus (medial)_R','Superior frontal gyrus (medial orbital)_L','Superior frontal gyrus (medial orbital)_R','Rectus gyrus_L','Rectus gyrus_R','Insula_L','Insula_R','Anterior cingulate gyrus_L','Anterior cingulate gyrus_R','Median cingulate gyrus_L','Median cingulate gyrus_R','Posterior cingulate gyrus_L','Posterior cingulate gyrus_R','Hippocampus_L','Hippocampus_R','Parahippocampal gyrus_L','Parahippocampal gyrus_R','Amygdala_L','Amygdala_R','Calcarine cortex_L','Calcarine cortex_R','Cuneus_L','Cuneus_R','Lingual gyrus_L','Lingual gyrus_R','Superior occipital gyrus_L','Superior occipital gyrus_R','Middle occipital gyrus_L','Middle occipital gyrus_R','Inferior occipital gyrus_L','Inferior occipital gyrus_R','Fusiform gyrus_L','Fusiform gyrus_R','Postcentral gyrus_L','Postcentral gyrus_R','Superior parietal gyrus_L','Superior parietal gyrus_R','Inferior parietal gyrus_L','Inferior parietal gyrus_R','Supramarginal gyrus_L','Supramarginal gyrus_R','Angular gyrus_L','Angular gyrus_R','Precuneus_L','Precuneus_R','Paracentral lobule_L','Paracentral lobule_R','Caudate_L','Caudate_R','Putamen_L','Putamen_R','Pallidum_L','Pallidum_R','Thalamus_L','Thalamus_R','Heschl gyrus_L','Heschl gyrus_R','Superior temporal gyrus_L','Superior temporal gyrus_R','Temporal pole (superior)_L','Temporal pole (superior)_R','Middle temporal gyrus_L','Middle temporal gyrus_R','Temporal pole (middle)_L','Temporal pole (middle)_R','Inferior temporal gyrus_L','Inferior temporal gyrus_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R','Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R','Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R','Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R','Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R','Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8','Vermis_9','Vermis_10'};
clc
for jj=1:3
    for ii=1:3
        [cause.extra.links,cause.extra.ind] = sort(sum(grid.extra{ii,jj},1),'descend');
        [effect.extra.links,effect.extra.ind] = sort(sum(grid.extra{ii,jj},2),'descend');
        
        [cause.missing.links,cause.missing.ind] = sort(sum(grid.missing{ii,jj},1),'descend');
        [effect.missing.links,effect.missing.ind] = sort(sum(grid.missing{ii,jj},2),'descend');
        
        subtype_name = {'extra','missing'};
        
        for subtype=1:2
            AAL_index = cause.(subtype_name{subtype}).ind(1:RANK_NUMBER)';
            roi_list = {AAL_116.name_full{AAL_index}}';
            no_of_links = cause.extra.links(1:RANK_NUMBER)';
            
            cause_table = table(AAL_index,roi_list,no_of_links);
            H(ii,jj).cause.(subtype_name{subtype})= cause_table;
            fprintf(['Cause table:',subtype_name{subtype},'\n'])
            disp(cause_table)
            
            AAL_index = effect.(subtype_name{subtype}).ind(1:RANK_NUMBER);
            roi_list = {AAL_116.name_full{AAL_index}}';
            no_of_links = effect.extra.links(1:RANK_NUMBER);
            effect_table = table(AAL_index,roi_list,no_of_links);
            H(ii,jj).effect.(subtype_name{subtype})= effect_table;
            fprintf(['Effect table:',subtype_name{subtype},'\n'])
            disp(effect_table)
        end
    end
end
clc
setting = {'D2K','S2K','C18K'};
ranking = {'1','2','3'};
cc_list = {'cause','effect'};
mm_list = {'extra','missing'};
for ii=1:3 % setting
    for jj=1:3 % ranking
        for cc=1:2
            for mm=1:2
                fprintf(['Setting ',setting{ii},', Rank: ',ranking{jj},' Causal type: ',cc_list{cc},', Abnormal type: ',mm_list{mm},'\n'])
                disp(H(ii,jj).(cc_list{cc}).(mm_list{mm}).no_of_links);
            end
        end
        
    end
end
%% D2K
clear extracted_index
fm =1;
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_D_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
type_connection = {'cause','effect'};
type_abnormal = {'extra','missing'};
for ii=1:RANK_NUMBER
    for jj=1:length(type_connection)
        for kk=1:length(type_abnormal)
            selected_index{ii,jj,kk}=H(fm,1).(type_connection{jj}).(type_abnormal{kk}).AAL_index(ii);
            switch type_connection{jj}
                case 'cause'
                    extracted_index{ii,jj,kk} = find(grid.(type_abnormal{kk}){fm,1}(:,selected_index{ii,jj,kk}));
                    
                case 'effect'
                    extracted_index{ii,jj,kk} = find(grid.(type_abnormal{kk}){fm,1}(selected_index{ii,jj,kk},:))';
            end
        end
    end
end
% sort connectivity according to extracted index
no_GC = 10;
for ii=1:RANK_NUMBER
    for jj=1:length(type_connection)
        for kk=1:length(type_abnormal)
            if jj==1 % cause
                switch type_abnormal{kk}
                    case 'extra' % sort in ADHD
                        [tmp,I] = sort(GC.total(extracted_index{ii,jj,kk},selected_index{ii,jj,kk},2),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                    case 'missing' % sort in TDC
                        [tmp,I] = sort(GC.total(extracted_index{ii,jj,kk},selected_index{ii,jj,kk},1),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                end
            else % effect
                switch type_abnormal{kk}
                    case 'extra' % sort in ADHD
                        [tmp,I] = sort(GC.total(selected_index{ii,jj,kk},extracted_index{ii,jj,kk},2),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                    case 'missing' % sort in TDC
                        [tmp,I] = sort(GC.total(selected_index{ii,jj,kk},extracted_index{ii,jj,kk},1),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                end
            end
        end
    end
end
%% Check D2K
ii=1; % Rank 1
jj=1; % cause
kk = 1; % extra
tile = tiledlayout(2,2);
for jj=1:2
    for kk=1:2
        nexttile;
        sorted_index = GC_index{ii,jj,kk};
        if jj==1
            a1 = GC.total(sorted_index,selected_index{ii,jj,kk},1); % TDC
            a2 = GC.total(sorted_index,selected_index{ii,jj,kk},2) ;% ADHD
            plot([a1 a2])
        else
            a1 = GC.total(selected_index{ii,jj,kk},sorted_index,1)'; % TDC
            a2 = GC.total(selected_index{ii,jj,kk},sorted_index,2)' ;% ADHD
            plot([a1 a2])
        end
        title([type_connection{jj},type_abnormal{kk}])
        legend('TDC','ADHD')
    end
end

% result verified the descending trend of GC strength given by GC_index{ii,jj,kk}
%% D2K: Effective connectivity differences visualization
clc

jj=1;
kk=2;
cnt=0;
for ii=1:6
    a1 = AAL_116.name_full(selected_index{ii,jj,kk});
    a2 = AAL_116.name_full(GC_index{ii,jj,kk})';
    for nn=1:no_GC
        cnt=cnt+1;
        ARR{cnt,1} = a1;
        ARR{cnt,2} = a2{nn};
        
    end
end
%% D2K: Effective connectivity differences visualization
jj=2;
kk=2;
cnt=0;
for ii=1:6
    a1 = AAL_116.name_full(selected_index{ii,jj,kk});
    a2 = AAL_116.name_full(GC_index{ii,jj,kk})';
    for nn=1:no_GC
        cnt=cnt+1;
        ARR{cnt,1} = a2{nn};
        ARR{cnt,2} = a1;
        
    end
end

%%












%% S2K
clear extracted_index
fm =2;
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_2K_S_unfiltered_timecorrected_LLHcorrected')
GC.total = M.model(M.index.eBIC).GC;
type_connection = {'cause','effect'};
type_abnormal = {'extra','missing'};
for ii=1:RANK_NUMBER
    for jj=1:length(type_connection)
        for kk=1:length(type_abnormal)
            selected_index{ii,jj,kk}=H(fm,1).(type_connection{jj}).(type_abnormal{kk}).AAL_index(ii);
            switch type_connection{jj}
                case 'cause'
                    extracted_index{ii,jj,kk} = find(grid.(type_abnormal{kk}){fm,1}(:,selected_index{ii,jj,kk}));
                    
                case 'effect'
                    extracted_index{ii,jj,kk} = find(grid.(type_abnormal{kk}){fm,1}(selected_index{ii,jj,kk},:))';
            end
        end
    end
end
% sort connectivity according to extracted index
no_GC = 10;
for ii=1:RANK_NUMBER
    for jj=1:length(type_connection)
        for kk=1:length(type_abnormal)
            if jj==1 % cause
                switch type_abnormal{kk}
                    case 'extra' % sort in ADHD
                        [tmp,I] = sort(GC.total(extracted_index{ii,jj,kk},selected_index{ii,jj,kk},2),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                    case 'missing' % sort in TDC
                        [tmp,I] = sort(GC.total(extracted_index{ii,jj,kk},selected_index{ii,jj,kk},1),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                end
            else % effect
                switch type_abnormal{kk}
                    case 'extra' % sort in ADHD
                        [tmp,I] = sort(GC.total(selected_index{ii,jj,kk},extracted_index{ii,jj,kk},2),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                    case 'missing' % sort in TDC
                        [tmp,I] = sort(GC.total(selected_index{ii,jj,kk},extracted_index{ii,jj,kk},1),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                end
            end
        end
    end
end
%% Check S2K
ii=1; % Rank 1
jj=1; % cause
kk = 1; % extra
tile = tiledlayout(2,2);
for jj=1:2
    for kk=1:2
        nexttile;
        sorted_index = GC_index{ii,jj,kk};
        if jj==1
            a1 = GC.total(sorted_index,selected_index{ii,jj,kk},1); % TDC
            a2 = GC.total(sorted_index,selected_index{ii,jj,kk},2) ;% ADHD
            plot([a1 a2])
        else
            a1 = GC.total(selected_index{ii,jj,kk},sorted_index,1)'; % TDC
            a2 = GC.total(selected_index{ii,jj,kk},sorted_index,2)' ;% ADHD
            plot([a1 a2])
        end
        title([type_connection{jj},type_abnormal{kk}])
        legend('TDC','ADHD')
    end
end

% result verified the descending trend of GC strength given by GC_index{ii,jj,kk}
%% S2K: Effective connectivity differences visualization
clc

jj=1;
kk=2;
cnt=0;
for ii=1:6
    a1 = AAL_116.name_full(selected_index{ii,jj,kk});
    a2 = AAL_116.name_full(GC_index{ii,jj,kk})';
    for nn=1:no_GC
        cnt=cnt+1;
        ARR{cnt,1} = a1;
        ARR{cnt,2} = a2{nn};
        
    end
end
%% S2K: Effective connectivity differences visualization
jj=2;
kk=2;
cnt=0;
for ii=1:6
    a1 = AAL_116.name_full(selected_index{ii,jj,kk});
    a2 = AAL_116.name_full(GC_index{ii,jj,kk})';
    for nn=1:no_GC
        cnt=cnt+1;
        ARR{cnt,1} = a2{nn};
        ARR{cnt,2} = a1;
        
    end
end
%%














%% C18K
clear extracted_index
fm =3;
load('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_18K_C_unfiltered_LLHcorrection','M')
GC.total(:,:,1) = mean(M.TDC.model(M.TDC.index.eBIC).GC,3);
GC.total(:,:,2) = mean(M.ADHD_C.model(M.ADHD_C.index.eBIC).GC,3);
type_connection = {'cause','effect'};
type_abnormal = {'extra','missing'};
for ii=1:RANK_NUMBER
    for jj=1:length(type_connection)
        for kk=1:length(type_abnormal)
            selected_index{ii,jj,kk}=H(fm,1).(type_connection{jj}).(type_abnormal{kk}).AAL_index(ii);
            switch type_connection{jj}
                case 'cause'
                    extracted_index{ii,jj,kk} = find(grid.(type_abnormal{kk}){fm,1}(:,selected_index{ii,jj,kk}));
                    
                case 'effect'
                    extracted_index{ii,jj,kk} = find(grid.(type_abnormal{kk}){fm,1}(selected_index{ii,jj,kk},:))';
            end
        end
    end
end
% sort connectivity according to extracted index
no_GC = 3;
for ii=1:RANK_NUMBER
    for jj=1:length(type_connection)
        for kk=1:length(type_abnormal)
            if jj==1 % cause
                switch type_abnormal{kk}
                    case 'extra' % sort in ADHD
                        [tmp,I] = sort(GC.total(extracted_index{ii,jj,kk},selected_index{ii,jj,kk},2),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                    case 'missing' % sort in TDC
                        [tmp,I] = sort(GC.total(extracted_index{ii,jj,kk},selected_index{ii,jj,kk},1),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                end
            else % effect
                switch type_abnormal{kk}
                    case 'extra' % sort in ADHD
                        [tmp,I] = sort(GC.total(selected_index{ii,jj,kk},extracted_index{ii,jj,kk},2),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                    case 'missing' % sort in TDC
                        [tmp,I] = sort(GC.total(selected_index{ii,jj,kk},extracted_index{ii,jj,kk},1),'descend');
                        GC_index{ii,jj,kk} = extracted_index{ii,jj,kk}(I(1:no_GC));
                end
            end
        end
    end
end
%% Check C18K
ii=1; % Rank 1
jj=1; % cause
kk = 1; % extra
tile = tiledlayout(2,2);
for jj=1:2
    for kk=1:2
        nexttile;
        sorted_index = GC_index{ii,jj,kk};
        if jj==1
            a1 = GC.total(sorted_index,selected_index{ii,jj,kk},1); % TDC
            a2 = GC.total(sorted_index,selected_index{ii,jj,kk},2) ;% ADHD
            plot([a1 a2])
        else
            a1 = GC.total(selected_index{ii,jj,kk},sorted_index,1)'; % TDC
            a2 = GC.total(selected_index{ii,jj,kk},sorted_index,2)' ;% ADHD
            plot([a1 a2])
        end
        title([type_connection{jj},type_abnormal{kk}])
        legend('TDC','ADHD')
    end
end

% result verified the descending trend of GC strength given by GC_index{ii,jj,kk}
%% C18K: Effective connectivity differences visualization
clc
clear ARR
jj=1;
kk=2;
cnt=0;
for ii=1:6
    a1 = AAL_116.name_full(selected_index{ii,jj,kk});
    a2 = AAL_116.name_full(GC_index{ii,jj,kk})';
    for nn=1:no_GC
        cnt=cnt+1;
        ARR{cnt,1} = a1;
        ARR{cnt,2} = a2{nn};
        
    end
end
%% D2K: Effective connectivity differences visualization
jj=2;
kk=2;
cnt=0;
for ii=1:6
    a1 = AAL_116.name_full(selected_index{ii,jj,kk});
    a2 = AAL_116.name_full(GC_index{ii,jj,kk})';
    for nn=1:no_GC
        cnt=cnt+1;
        ARR{cnt,1} = a2{nn};
        ARR{cnt,2} = a1;
        
    end
end