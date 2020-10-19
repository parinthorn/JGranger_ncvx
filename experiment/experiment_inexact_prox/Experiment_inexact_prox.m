clear
clc
clf
p=5;
K = 10;

q = 0.5;
x0 = 0.1*ones(p*K,1);
x0 = reshape(x0,[p,K]);
z = 1+1*randn(p,K);

T = 2000;
obj_equal = zeros(T,1);


%% Parinthorn fixed point Derivation [CONVERGE]
a1 = 3;
a2 = 0;
obj_prox = @(x) a2*norm(x,2)^(q)+a1*sum(norms(reshape(x,[p,K]),2,1).^(q))+ ...
    0.5*norm(x(:)-z(:),2)^2;
equality_obj = @(x,kk) norm((a2*q*sqrt(norm(x,'fro')^2)^(q-2)+ ...
    a1*q*sqrt(norm(x(:,kk),'fro')^2)^(q-2)+1)*(x(:,kk))-z(:,kk),2)^2;
PARAMETERS  = [2,q,p];
prox_x = prox_pq_eff(z(:),a1,PARAMETERS);
x = x0;
for tt=1:T
    for kk=1:K
        rcp = a2*q*(norm(x(:),2))^(q-2)+a1*q*(norm(x(:,kk),2))^(q-2);
        x(:,kk) = (z(:,kk)-rcp*x(:,kk));
    end
    
    for kk=1:K
        obj_equal(tt) = obj_equal(tt) + equality_obj(x,kk);
    end
    tmp = x(:);
    tmp(prox_x==0) = 0;
    err(tt) = abs(obj_prox(tmp)-obj_prox(prox_x));
    sol_diff(tt) = norm(tmp-prox_x)/norm(prox_x);
    if (tt>2) && (abs((obj_equal(tt) - obj_equal(tt-1))) <1e-6)
        disp('equal')
        break
    end 
end
obj_equal = obj_equal(1:tt);
figure(1)
subplot(2,3,6)
rr = zeros(p,K);
for kk=1:K
    rr(:,kk) = equality_obj(x,kk);
end
plot(rr(:))

% Find Sparse solution
reference_loss = obj_prox(x);
zero_loss = obj_prox(zeros(p*K,1));
z_ind = [];
h = x;
for kk=1:K
    tmp = h(:,kk);
    h(:,kk) = 0;
    if obj_prox(h)<reference_loss
        z_ind = [z_ind kk];
        h(:,kk) = tmp;
    else
        h(:,kk) = tmp;
    end
end
nz_ind = setdiff(1:K,z_ind);
h=x;
h(:,z_ind) = 0;

if zero_loss<obj_prox(h)
    h = zeros(p*K,1);
end

figure(1)
subplot(2,3,1)
semilogy(abs(diff(obj_equal)))
title('sum square error difference')
subplot(2,3,2)
plot(((obj_equal)))
title('sum square error')


subplot(2,3,3)
imagesc([x(:) prox_x h(:)])
h = x(:);
SSE = sum((x(prox_x~=0)-prox_x(prox_x~=0)).^2)/sum(prox_x(prox_x~=0).^2);
title(sprintf('non-zero SSE: %f',SSE))

subplot(2,3,4)
plot(err)
% semilogy(abs(diff(err)))
subplot(2,3,5)
plot(sol_diff)


% semilogy(abs(diff(sol_diff)))
 %% JSS fixed point Derivation [CONVERGE]
% obj_equal = zeros(T,1);
% x=x0;
% for tt=1:T
%     for kk=1:K
%         rcp = a2*q*(norm(x(:),2))^(q-2)+1;
%         x(:,kk) = (z(:,kk)-a1*q*((norm(x(:,kk),2))^(q-2))*x(:,kk))/rcp;
%     end
%     for kk=1:K
%         obj_equal(tt) = obj_equal(tt) + equality_obj(x,kk);
%     end
%     tmp = x(:);
%     err(tt) = obj_prox(tmp(prox_x~=0))-obj_prox(prox_x(prox_x~=0));
%     if (tt>2) && (abs((obj_equal(tt) - obj_equal(tt-1))) <1e-6)
%         disp('equal')
%         break
%     end
%     
%     
%     
% end
% obj_equal = obj_equal(1:tt);
% 
% % Find Sparse solution
% reference_loss = obj_prox(x);
% zero_loss = obj_prox(zeros(p*K,1));
% z_ind = [];
% h = x;
% for kk=1:K
%     tmp = h(:,kk);
%     h(:,kk) = 0;
%     if obj_prox(h)<reference_loss
%         z_ind = [z_ind kk];
%         h(:,kk) = tmp;
%     else
%         h(:,kk) = tmp;
%     end
% end
% nz_ind = setdiff(1:K,z_ind);
% h=x;
% h(:,z_ind) = 0;
% 
% if zero_loss<obj_prox(h)
%     h = zeros(p*K,1);
% end
% 
% figure(1)
% subplot(2,3,4)
% semilogy(abs(diff(obj_equal)))
% title('sum square error difference')
% subplot(2,3,5)
% plot(((obj_equal)))
% title('sum square error')
% subplot(2,3,6)
% imagesc([x(:) prox_x h(:)])
% h = x(:);
% SSE = sum((x(prox_x~=0)-prox_x(prox_x~=0)).^2)/sum(prox_x(prox_x~=0).^2);
% title(sprintf('non-zero SSE: %f',SSE))
% 
