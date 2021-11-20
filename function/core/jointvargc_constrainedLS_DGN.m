function [x,flag] = jointvargc_constrainedLS_DGN(G,b,nz_ind)
% This function computes the constraint least-square on the DGN formulation
% which is to solve min_{x} ||Gx-b||_2 constrained on the sparsity x 
% obtained from DGN. If the least-square solution is undefined, we replaced
% with the L2 regularized solution with penalty parameter 0.01.
% Input
% G: vectorized matrix H of all k
% b: vectorized matrix Y of all k
% nz_ind: non-zero index of learned VAR coefficients
% Output
% x: constrained least-square version of VAR coefficients
% flag: 0 -> least-square
%      -1 -> L2 regularized solution
%
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce
x = zeros(size(G,2),1);
gtmp = G(:,nz_ind);
if length(b)>length(nz_ind)
  x(nz_ind) = (gtmp\b);
  flag = 0; % normal least square
else
  x(nz_ind) = (gtmp'*gtmp+0.01*eye(length(nz_ind)))\(gtmp'*b);
  flag=-1; % L2 reg solution
end

end
