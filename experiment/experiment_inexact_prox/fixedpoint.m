function[u1,u2] = fixedpoint(x,varargin)
global p K q a1 a2 z
% u1 = G(x), 
% u2 = G(x)
% x is p x K
% indx is the index of nonzero column 

w = x; 
if isempty(varargin) % if sparsity of block in x is not specified
    u1 = x; u2 = w;
    indx = 1:size(x,2);
else % use specified column index
    indx = varargin{1};
    u1 = zeros(p,K); % sequential fixed point
    u2 = zeros(p,K); % parallel fixed point
    
    if isempty(indx)  % if indx = [], return the zero solution
        return
    end
end



    for kk = indx  % run fixed point only on specified column
        Dx = a2*q*(norm(x(:),2))^(q-2)+a1*q*(norm(x(:,kk),2))^(q-2);
        u1(:,kk) = z(:,kk)-Dx*x(:,kk);
        x(:,kk) = u1(:,kk); % update the whole vector x
        
        Dx = a2*q*(norm(w(:),2))^(q-2)+a1*q*(norm(w(:,kk),2))^(q-2);
        u2(:,kk) = z(:,kk)-Dx*w(:,kk);
    end
end
