function z = normpq_adaptive(x,p,q,gLen,a)
n = length(x);
gNo = n/gLen;
tmp = reshape(x,gLen,gNo);
tmp = (sum(abs(tmp).^(p),1)).^(1/p);
z = sum(a'.*(tmp.^q));
% z = norm(,q)^q;

end