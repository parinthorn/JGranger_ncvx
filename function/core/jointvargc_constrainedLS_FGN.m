function x = jointvargc_constrainedLS_FGN(G,y,D,Dx,P,Px,varargin)
if ~isempty(varargin)
    options = optimset('Display',varargin{1});
else
    options = optimset('Display','notify');
end
tmpD = D;
tmpD(Dx~=0,:) = [];
tmpP = P;
tmpP(Px~=0,:) = [];
b=zeros(size(tmpD,1)+size(tmpP,1),1);
% B = blkdiag(tmpD,tmpP);
% disp(size(tmpD));
% disp(size(tmpP));
B = [tmpD;tmpP];
x = lsqlin(G,y,[],[],B,b,[],[],[],options);
end