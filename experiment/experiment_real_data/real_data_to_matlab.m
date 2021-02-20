x.nyu=1;
x.pku1=1;
x.pku2=1;
x.pku3=1;
list_u = {'nyu','pku1','pku2','pku3'};
%%
for ww=1:length(list_u)
subject_list.(list_u{ww}) = cell(length(x.(list_u{ww})),1);
for ii = 1:length(x.(list_u{ww}))
    kk = int2str(x.(list_u{ww})(ii));
    if length(kk)~=7
        subject_list.(list_u{ww}){ii} = ['00',kk];
    else
        subject_list.(list_u{ww}){ii} = kk;
    end
end
end
%%
folder_name = {'NYU','Peking_1','Peking_2','Peking_3'};
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\ADHD200_AAL_TCs_filtfix\matlab_format\';
for ii=2:length(folder_name)
    for jj=1:length(subject_list.(list_u{ii}))
    a = dir(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\ADHD200_AAL_TCs_filtfix\',folder_name{ii},'\',subject_list.(list_u{ii}){jj},'\sfnwmrda', ...
        subject_list.(list_u{ii}){jj},'_session_*']);
    file_name = {a.name};
    folder_name_full = {a.folder};
    y = [];
    for nn=1:length(file_name)
       [data,header,raw] = tsvread([folder_name_full{nn},'\',file_name{nn}]);
       y = [y data(2:end,3:end)'];
       
    end
    fprintf([folder_name{ii},'_sfnwmrda',subject_list.(list_u{ii}){jj},' length:\t %d \n'],size(y,2))
    save([outpath,folder_name{ii},'_sfnwmrda',subject_list.(list_u{ii}){jj}],'y')
%     disp({a.name})
    
    end
end
%% unfilter
folder_name = {'NYU'};
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\ADHD200_AAL_TCs_filtfix\matlab_format\';
for ii=1:length(folder_name)
    for jj=1:length(subject_list.(list_u{ii}))
    a = dir(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\ADHD200_AAL_TCs_filtfix\',folder_name{ii},'\',subject_list.(list_u{ii}){jj},'\snwmrda', ...
        subject_list.(list_u{ii}){jj},'_session_*']);
    file_name = {a.name};
    folder_name_full = {a.folder};
    y = [];
    for nn=1:length(file_name)
       [data,header,raw] = tsvread([folder_name_full{nn},'\',file_name{nn}]);
       y = [y data(2:end,3:end)'];
       
    end
    fprintf([folder_name{ii},'_snwmrda',subject_list.(list_u{ii}){jj},' length:\t %d \n'],size(y,2))
    save([outpath,folder_name{ii},'_snwmrda',subject_list.(list_u{ii}){jj}],'y')
%     disp({a.name})
    
    end
end