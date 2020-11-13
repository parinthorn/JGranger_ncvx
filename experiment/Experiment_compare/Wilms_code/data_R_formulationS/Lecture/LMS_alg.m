function [y,w,bias,J,step_size] = LMS_alg(d,x,filter_order,w_opt,tracking,p,R)
M = filter_order;
w = zeros(M,1);
u = 0.001;
beta = 0.00001;
e = zeros(length(d),1);

% Jmin = var(d);
% w =w_opt;
if tracking
    w =zeros(filter_order,length(d));
end
for k=M:length(d)
    xx = x(k:-1:k-M+1);

    step_size(k) = u*(1/(beta+norm(xx)^2));
    if tracking
            y(k) = w(:,k-1)'*xx;
    e(k) = d(k)-y(k);
        w(:,k) = w(:,k-1)+ u*e(k)*(1/(beta+norm(xx)^2))*xx;
        err = w(:,k)-w_opt;
        bias(k) = norm(err)/norm(w_opt);
    else
            y(k) = w'*xx;
    e(k) = d(k)-y(k);
        w = w + u*e(k)*(1/(beta+norm(xx)^2))*xx;
        err = w-w_opt;
        bias = norm(err)/norm(w_opt);
    end
    
    
    J(k) = var(d)-p'*(R\p)+ (err)'*R*err;
end

end