clear all; clc;

global p K q a1 a2 z
% prox_obj_fun = @(x) a2*norm(x(:),'fro')^(q)+a1*sum(norms(x,2,1).^(q))+ ...
%     0.5*norm(x(:)-z(:),'fro')^2; % size(x) = [p,K]
% sse_fun = @(x,kk) norm((a2*q*sqrt(norm(x,'fro')^2)^(q-2)+ ...
%     a1*q*sqrt(norm(x(:,kk),'fro')^2)^(q-2)+1)*(x(:,kk))-z(:,kk),'fro')^2; % sum square error of optimality condition

%% Test fixed-point iterations (all subsets)
% This will generate all fixed-point given some blocks of x (varied in
% column)
% I hypothesized that for a fixed (a1,a2) the divergent solution (for a given sparsity) cannot be optimal solution

p = 5; K = 4; q = 0.5;

a1 = 5*abs(randn(1)); % can vary a scaling to see if the fixed point converge
a2 = 6*abs(randn(1)); % can vary a scaling to see if the fixed point converge
z= 1*randn(p,K); x0 = randn(p,K); 

ITER = 50;
xlog = zeros(p,K,2^K,ITER); % historical values of x 
proxcost = zeros(1,2^K); % cost objective of prox problem (1/2)|| x - z||^2 + a1*g1(x) + a2*g2(x)
sseval = zeros(K,2^K);
indx = genallsubset(K); % all indices of nonzero column in x
% 
for kk=1:2^K
    x = x0;
    for ii=1:ITER
        [x,~] = fixedpoint(x,indx{kk});         
        xlog(:,:,kk,ii) = x;
    end
    proxcost(kk) = 0.5*(norm(x-z,'fro'))^2 + a2*norm(x(:),'fro')^(q)+a1*sum(sqrt(sum(x.^2,1)).^q);
    for jj=1:K
    tmpval=norm((a2*q*sqrt(norm(x,'fro')^2)^(q-2)+a1*q*sqrt(norm(x(:,jj),'fro')^2)^(q-2)+1)*(x(:,jj))-z(:,jj),'fro')^2;
    if ~isnan(tmpval)
        
    sseval(jj,kk) = sseval(jj,kk)+tmpval;
    end
    end
end
dx = diff(xlog(:,:,:,end-10:end),1,4); % difference along iteration

oscil_ind = find(squeeze(sum(sum(sum(abs(dx),4),1),2)) >= 1e-3)% sum difference along p,K,ITER dimension
% this index correspond to the subsets that we should not compare the
% scores (because iterations do not converge)

figure(1); % plot cell K x 2^K. Each column specify the result of forcing some nonzero blocks. Each row is kth (out of K) block of x.
vecx = (1:ITER)'; 
k = 0;
for jj=1:2
    for kk=1:2^K
        k = k+1;
        subplot(K,2^K,k); x2plot = squeeze(xlog(:,jj,kk,:)); % p x ITER
        plot(vecx,x2plot'); % plot x versus ITER
        if jj == 1
            title(num2str(indx{kk})); % specify the non zero block index of x
        end
        if kk == 1
            ylabel(['block #',num2str(jj)]);
        end
    end
end
[val,index_prox] = min(proxcost)
sd_xlog = std(xlog(:,:,:,end-5:end),0,4);
figure(2)
hh = squeeze(sum(sd_xlog,1));
imagesc(max(hh))
%% Test fixed-point iterations (compare updating schemes)
% 
% a1 = 1*abs(randn(1)); 
% a2 = abs(randn(1));
% z= randn(p,K); x0 = randn(p,K); y0 = randn(p,K);
% x = x0; y = x0; % y has same initial as x but updated in different ways
% xlog = []; ylog = [];
% ITER = 50;
% for ii=1:ITER
% %     [x,y] = fixedpoint(x);
%     [x,y] = fixedpoint(x,[1 3 4]); 
% %     xlog = [xlog x(:)];
% %     ylog = [ylog y(:)];
%     xlog(:,:,ii) = x;
%     ylog(:,:,ii) = y;
% end
% 
% vecx = (1:ITER)'; 
% figure(1);
% k = 0;
% for ii=1:p
%     for jj=1:K
%         k = k+1;
%         subplot(p,K,k); plot(vecx,squeeze(xlog(ii,jj,:)),vecx,squeeze(ylog(ii,jj,:)));
%     end
%     legend('sequential','parallel');
% end
% result is that 
% 1) when a1, a2 = randn, the last entry of x and y  to the same value (so updating scheme may not
% matter) but the iteration may not converge for all entries
% 2) when a1 , a2 are large (100*randn), x and y updates may not be the
% same and osciilate

%% Test contractive constant
% z= 100*randn(p,K); 
% cnt = 0;
% fprintf('|| G(x) - G(y) ||  || x - y ||   constant \n');
% for ii=1:1000
% %     x = randn(p,K); y = randn(p,K); % contractive 
%     x = 0.001*randn(p,K); y = randn(p,K); % not contractive
% %     x = 100*randn(p,K); y = randn(p,K); % contractive
% %     x = 100*randn(p,K); y = 1000*randn(p,K); % contractive
%     gx = fixedpoint1(x); gy = fixedpoint1(y);
%     normgxy = norm(gx-gy,'fro'); normxy = norm(x-y,'fro');
%     cnt = cnt+ (normgxy > normxy);
%     fprintf('%2.3f %2.3f %2.5f \n',[normgxy normxy normgxy/normxy]); 
% end
% cnt

% the result depends on the value of x, y, so G might be contractive only in
% some region


