clear
clc
R = [0 0 0]; %[forward backward greedy]
for tt=1:1
    p=5;
    K = 5;
    q = 0.5;
    % a1 = 3;a2 = 4;
    z = 1+1*randn(p,K);
    normz_vect_pK = ((((sum((reshape(z,[p*K,1])).^2,1))).^(0.5)')./1.5).^(3/2);
    normz_vect_p = ((((sum((z).^2,1))).^(0.5))'./1.5).^(3/2);
    lambda1 = linspace(0,1.01,50)*max(normz_vect_p);
    lambda2 = linspace(0,1.01,50)*normz_vect_pK;
    NZ_MAP = zeros(50,50);
    % a1 = lambda1(end);
    % a2 = 0;
    %     fprintf('%3s\t%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n','a1','a2', ...
    %         'obj forward','error forward','obj backward','error backward', 'obj greedy','error greedy');
    %%
    for ii=1:length(lambda1)
        for jj=1:length(lambda2)
            a1 = lambda1(ii);
            a2 = lambda2(jj);
            PARAMETERS  = [2,q,p];
            tic;
            prox_x = prox_pq_eff(z(:),a1,PARAMETERS);
            toc;
            prox_x = reshape(prox_x,[p,K]);
            %
            prox_obj_fun = @(x) a2*norm(x(:),'fro')^(q)+a1*sum(norms(x,2,1).^(q))+ ...
                0.5*norm(x(:)-z(:),'fro')^2;
            tic
            [x_forward,obj_prox_forward,obj_sse_forward] = best_subset_selection(z,a1,a2,q,'forward');
            toc;
            tic
            [x_backward,obj_prox_backward,obj_sse_backward] = best_subset_selection(z,a1,a2,q,'backward');
            toc;
            tic
            [x_optimal,obj_prox_optimal,obj_sse_optimal] = best_subset_selection(z,a1,a2,q,'optimal');
            toc
            tic
            [x_greedy,obj_prox_greedy,obj_sse_greedy] = best_subset_selection(z,a1,a2,q,'greedy');
            toc;
            yf=prox_obj_fun(x_forward);
            yb=prox_obj_fun(x_backward);
            yo=prox_obj_fun(x_optimal);
            yg=prox_obj_fun(x_greedy);
            % figure(1);
            % subplot(1,4,1)
            % imagesc(x_optimal)
            % subplot(1,4,2)
            % imagesc(x_forward)
            % subplot(1,4,3)
            % imagesc(x_backward)
            % subplot(1,4,4)
            % imagesc(x_greedy)
            % pause(0.1)
            if any(x_optimal~=0)
                NZ_MAP(ii,jj) = 1;
            end
            ef = norm(x_forward-x_optimal,'fro')/norm(x_optimal,'fro');
            eb = norm(x_backward-x_optimal,'fro')/norm(x_optimal,'fro');
            eg = norm(x_greedy-x_optimal,'fro')/norm(x_optimal,'fro');
            
            
            if (ef >=1e-3) &&(norm(x_forward-x_optimal,'fro')~=0)
                R(1) = R(1)+1;
            end
            if (eb >=1e-3)&&(norm(x_backward-x_optimal,'fro')~=0)
                R(2) = R(2)+1;
            end
            if (eg >=1e-3)&&(norm(x_greedy-x_optimal,'fro')~=0)
                R(3) = R(3)+1;
            end
            disp([ef eb eg;yf yb yg])
            disp(yo)
            %     fprintf('%3s\t%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n',ii,jj, ...
            %         yf,ef,yb,eb, yg,eg);
        end
    end
    
end