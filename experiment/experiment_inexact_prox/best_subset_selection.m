function [x,obj_prox,obj_sse] = best_subset_selection(z,a1,a2,qq,selection_type)
prox_obj_fun = @(x) a2*norm(x(:),'fro')^(qq)+a1*sum(sqrt(sum(x.^2,1)).^(qq))+ ...
    0.5*norm(x(:)-z(:),'fro')^2;
[p,K] = size(z);
switch selection_type
    case 'forward' % begin with null model
        reference_loss = prox_obj_fun(zeros(p,K));
        sol_found = 1;
        nz_ind = [];
        z_ind = 1:K;
        while sol_found
            obj_prox_val = ones(K,1)*inf;
            h = zeros(p,K,length(z_ind));
            for kk=z_ind % find nonzero candidate from zero entries
                nz_ind_tmp = [nz_ind kk];
                [h(:,:,kk),prox_obj,~] = fix_point_prox(z,nz_ind_tmp,a1,a2,qq);
                obj_prox_val(kk) =prox_obj(end);
                nz_ind_tmp = [];
            end
            [minval,minval_ind] = min(obj_prox_val);
            if minval<reference_loss
                reference_loss = minval;
                nz_ind = [nz_ind minval_ind];
                z_ind(z_ind==minval_ind) = [];
                if isempty(z_ind)
                    sol_found=0;
                end
            else
                
                sol_found = 0;
            end
        end
        [x,prox_obj,obj_sse] = fix_point_prox(z,nz_ind,a1,a2,qq);
        obj_prox = prox_obj(end);
        
    case 'backward'
        
        [~,prox_obj,~] = fix_point_prox(z,1:K,a1,a2,qq);
        reference_loss =prox_obj(end);
        sol_found = 1;
        z_ind = [];
        nz_ind = 1:K;
        while sol_found
            obj_prox_val = ones(K,1)*inf;
            h = zeros(p,K,length(nz_ind));
            for kk=nz_ind % find zero candidate from nonzero entries
                nz_ind_tmp = nz_ind;
                nz_ind_tmp(nz_ind_tmp==kk) = [];
                [h(:,:,kk),prox_obj,~] = fix_point_prox(z,nz_ind_tmp,a1,a2,qq);
                obj_prox_val(kk) =prox_obj(end);
            end
            [minval,minval_ind] = min(obj_prox_val);
            if minval<=reference_loss
                
                reference_loss = minval;
                nz_ind = [nz_ind minval_ind];
                z_ind(z_ind==minval_ind) = [];
                if isempty(z_ind)
                    sol_found=0;
                end
            else
                sol_found = 0;
            end
        end
        [x,prox_obj,obj_sse] = fix_point_prox(z,nz_ind,a1,a2,qq);
        obj_prox = prox_obj(end);
    case 'optimal'
        reference_loss = prox_obj_fun(zeros(p,K));
        h = zeros(p,K);
        for kk=1:K % choose kk from K
            ind = combnk(1:K,kk); % all combination
            %     x_tmp = h;
            obj_prox_val = zeros(size(ind,1),1);
            x = zeros(p,K,size(ind,1));
            for jj = 1:size(ind,1)
                nz_ind = ind(jj,:);
                [x(:,:,jj),prox_obj,obj_sse] = fix_point_prox(z,nz_ind,a1,a2,qq);
                obj_prox_val(jj) =prox_obj(end);
            end
            [~,choose_ind] = min(obj_prox_val);
            if obj_prox_val(choose_ind)<reference_loss
                h = x(:,:,choose_ind);
                reference_loss = prox_obj_fun(h);
            end
        end
        x = h;
        obj_prox = prox_obj(end);
    case 'greedy'
        [x_greedy,prox_obj,obj_sse] = fix_point_prox(z,1:1:K,a1,a2,qq);
        obj_prox =prox_obj(end);
        % Find Sparse solution
        reference_loss = prox_obj(end);
        zero_loss = prox_obj_fun(zeros(p,K));
        z_ind = [];
        h = x_greedy;
        for kk=1:K
            tmp = h(:,kk);
            h(:,kk) = 0;
            if prox_obj_fun(h)<reference_loss
                z_ind = [z_ind kk];
                h(:,kk) = tmp;
            else
                h(:,kk) = tmp;
            end
        end
        nz_ind = setdiff(1:K,z_ind);
        h=x_greedy;
        h(:,z_ind) = 0;
        
        if zero_loss<prox_obj_fun(h)
            h = zeros(p,K);
            x = h;
        else
            [x,prox_obj,obj_sse] = fix_point_prox(z,setdiff(1:K,z_ind),a1,a2,qq);
            obj_prox = prox_obj(end);
        end
        
        
end

end