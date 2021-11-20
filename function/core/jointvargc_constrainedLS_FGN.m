function x = jointvargc_constrainedLS_FGN(G,b,D,Dx,P,Px,varargin)
% This function computes the constraint least-square on the DGN formulation
% which is to solve min_{x} ||Gx-b||_2 constrained on Px=0, Dx=0
% obtained from FGN. 
% Input
% G: vectorized matrix H of all k
% b: vectorized matrix Y of all k
% D: difference matrix
% Dx: difference matrix applies on results obtained from FGN
% P: projection matrix
% Px: projection matrix applies on results obtained from FGN
% Output
% x: constrained least-square version of VAR coefficients
% 
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce
if ~isempty(varargin)
    options = optimset('Display',varargin{1});
else
    options = optimset('Display','notify');
end
tmpD = D;
tmpD(Dx~=0,:) = [];
tmpP = P;
tmpP(Px~=0,:) = [];
z=zeros(size(tmpD,1)+size(tmpP,1),1);
B = [tmpD;tmpP];
x = lsqlin(G,b,[],[],B,z,[],[],[],options);
end