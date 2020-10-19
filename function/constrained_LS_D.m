function [x,flag] = constrained_LS_D(G,y,nz_ind)
x = zeros(size(G,2),1);
gtmp = G(:,nz_ind);
if length(y)>length(nz_ind)
  x(nz_ind) = (gtmp\y);
  flag = 0; % normal least square
else
  x(nz_ind) = (gtmp'*gtmp+0.01*eye(length(nz_ind)))\(gtmp'*y);
  flag=-1; % L2 reg solution
end

end
