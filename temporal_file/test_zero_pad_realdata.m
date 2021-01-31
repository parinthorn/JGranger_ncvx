
for ii = 1:length(tmp)
    kk = int2str(tmp(ii));
    if length(kk)~=7
        id_list{ii} = ['00',kk];
    else
        id_list{ii} = kk;
    end
end