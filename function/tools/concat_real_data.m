function y_K = concat_real_data(id_list,n,toggle,is_filter)
inpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\ADHD200_AAL_TCs_filtfix\matlab_format\';
switch toggle
    case 'nyu'
        name_u = 'NYU';
        T_max = 172;
    case 'pku1'
        name_u = 'Peking_1';
        T_max = 232;
    case 'pku2'
        name_u = 'Peking_2';
        T_max = 232;
    case 'pku3'
        name_u = 'Peking_3';
        T_max = 232;
end
K = length(id_list);
y_K = zeros(n,T_max,K);
<<<<<<< HEAD
if is_filter
    for kk=1:K
        load([inpath,name_u,'_','sfnwmrda',id_list{kk}])
        y_K(:,:,kk) = y(:,1:T_max);
    end
else
    for kk=1:K
        load([inpath,name_u,'_','snwmrda',id_list{kk}])
        y_K(:,:,kk) = y(:,1:T_max);
    end
=======
for kk=1:K
    load([inpath,name_u,'_','snwmrda',id_list{kk}])
    y_K(:,:,kk) = y(:,1:T_max);
>>>>>>> ca504cf053124889a8cbabc1e20efd61b16da0e4
end
end