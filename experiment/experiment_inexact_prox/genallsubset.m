
function[ind] = genallsubset(n)

ind{1} = [];
ii = 1;
for k=1:n
    num_nchoosek = nchoosek(n,k);
    tmp = nchoosek(1:n,k); % general all indices of length k
    for jj=1:num_nchoosek
        ii = ii + 1;
        ind{ii} = tmp(jj,:);
    end
end