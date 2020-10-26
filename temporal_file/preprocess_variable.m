clear
clc
mname = {'1','5'};
for m=1:length(mname)
    for rr=1:20
        fname = ['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\formulation_D_result\result_formulationD_',mname{m},'percent_lag1_K5_',int2str(rr)];
        load(fname);
        M.model = rmfield(M.model,'ind_VAR');
        for ii=1:30
            for jj=1:30
                
                for kk=1:5
                    M.model(ii,jj).ind_VAR{kk} = find(M.model(ii,jj).A(:,:,:,kk));
                end
                M.model(ii,jj).ind_VAR = M.model(ii,jj).ind_VAR;
            end
        end
        M.model = orderfields(M.model,[1 2 3 4 5 6 7 9 8]);
        save(fname,'M')
    end
end