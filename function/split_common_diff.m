function [ind_common,ind_differential] = split_common_diff(nz_ind,dim)
  n= dim(1);
  p= dim(2);
  K= dim(3);
  diag_ind = 1:n+1:n^2;
  ind_common = setdiff(1:n^2,diag_ind);
  ind_differential = cell(K,1);
  for kk=1:K
    ind_common = intersect(ind_common,ind_nz{kk})
  end
  for kk=1:K
    tmp = setdiff(ind_nz{kk},diag_ind); % remove diagonal parts
    ind_differential{kk} = setdiff(tmp,ind_common); % remove common part, so that
                                                    % the remaining is differential part
  end
end
