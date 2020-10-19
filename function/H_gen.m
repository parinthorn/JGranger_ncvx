function [H,Y] = H_gen(data,p)
[m,n] = size(data);
pm = p*m;
X = zeros(pm,n-p);
Y = data(:,p+1:end);
for j =1:p
    X(m*(j-1)+1:m*j,:) = data(:,p+1-j:n-j);
end
H = X;
end