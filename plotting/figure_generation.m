%% plotting results
% This script plots multiple bar graph
% input: score, size(score) = [#variation, #realization]
% graph: splits total, common, differential as 3 subgraphs
%        each graph represents one performance score, e.g. F1 score
%        each subgraph varies the varying factors, e.g. differential
%        density of 1% and 5%
clear;clc;clf;close all
inpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
skrip_path = './experiment/Experiment_compare/skrip_code_2/data_R_formulationD/';
% 
% load([inpath,'formulation_D_result_total'])
% our.total = R;
% load([inpath,'formulation_D_result_common'])
% our.common = R;
% load([inpath,'formulation_D_result_differential'])
% our.differential = R;

load([inpath,'formulation_D_result_all'])
our =R;

load([skrip_path,'skrip_formulationD_accuracy'])
score_list = {'TPR','FPR','ACC','F1','MCC'};
for ii=1:2
    for jj=1:length(score_list)
        ss = score_list{jj};
        ref.total.(ss)(ii,:) =score(ii).total.(ss);
        ref.common.(ss)(ii,:) = score(ii).common.(ss);
        ref.differential.(ss)(ii,:) = score(ii).differential.(ss);
    end
end
plot_seq = {'total','common','differential'};
diff_density = {'1%','5%'};

for fig = 1:length(score_list)
    figure(fig)
    
    for subfig = 1:length(plot_seq)
        subplot(1,length(plot_seq),subfig)
        T = [mean(our.(plot_seq{subfig}).(score_list{fig}),2) mean(ref.(plot_seq{subfig}).(score_list{fig}),2)];
        bar(T);
        set(gca,'xticklabel',diff_density,'FontSize',14);
        title(plot_seq{subfig})
        if subfig==1
            ylabel(score_list{fig})
        elseif subfig==2
            xlabel('differential density')
        elseif subfig==3
            legend('our','ref')
        end
    end
end
%%

%% K=50 case
clear;clc;clf;close all
inpath = 'G:/My Drive/0FROM_SHARED_DRIVE/THESIS/formulation_D_result/';
skrip_path = './experiment/Experiment_compare/skrip_code_2/data_R_formulationD/';
% 

load([inpath,'formulation_D_result_K50all'])
our =R;

score_list = {'TPR','FPR','ACC','F1','MCC'};
load([skrip_path,'skrip_formulationD_accuracy_K50'])
for ii=1:1
    for jj=1:length(score_list)
        ss = score_list{jj};
        ref.total.(ss)(ii,:) =score(ii).total.(ss);
        ref.common.(ss)(ii,:) = score(ii).common.(ss);
        ref.differential.(ss)(ii,:) = score(ii).differential.(ss);
    end
end
plot_seq = {'total','common','differential'};
diff_density = {'1%','5%'};

for fig = 1:length(score_list)
    figure(fig)
    
    for subfig = 1:length(plot_seq)
        subplot(1,length(plot_seq),subfig)
        T = [mean(our.(plot_seq{subfig}).(score_list{fig}),2) [mean(ref.(plot_seq{subfig}).(score_list{fig}),2);0]];
        bar(T);
        set(gca,'xticklabel',diff_density,'FontSize',14);
        title(plot_seq{subfig})
        if ~strcmp(score_list{fig},'FPR')
            ylim([0 1])
        else
            ylim([0 0.2])
        end
        if subfig==1
            ylabel(score_list{fig})
        elseif subfig==2
            xlabel('differential density')
        elseif subfig==3
            legend('our','ref')
        end
    end
end



