function [Lambdacrit_1,Lambdacrit_2, L1, L2] = gen_critical_lambdas(G,b, xLS,PARAMETER,qnorm,penalty_weight,formulation,varargin)
switch qnorm
    case 'cvx'
        q=1;
    case 'ncvx'
        q=1/2;
end
if isempty(varargin)
    gamma=q;
else
    gamma = varargin{1};
end
n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);
m1 = PARAMETER(4);
m2 = PARAMETER(5);

if ~strcmp(formulation,'fgn')
    switch penalty_weight
        case 'LS'
            [~,~,P1_m1,P1_m2,~] = weighted_projection(xLS,q,[n,p,K],gamma);
            
            Lambda0_1 = lambdamax_grouplasso_v2(G,b,m1,[n ,p ,K],P1_m1);
            Lambda0_2 = lambdamax_grouplasso_v2(G,b,m2,[n ,p ,K],P1_m2);
            
            weight_1 = max(P1_m1,[],2);
            weight_2 = max(P1_m2,[],2);
            
            Lambdacrit_1 = weight_1.*Lambda0_1; % sparsest weight vector
            Lambdacrit_2 = weight_2.*Lambda0_2; % sparsest weight vector
            
            
        case 'uniform'
            Lambdacrit_1 = lambdamax_grouplasso_v2(G,b,m1,[n ,p ,K]);
            Lambdacrit_2 = lambdamax_grouplasso_v2(G,b,m2,[n ,p ,K]);
            
        otherwise
            error('Please provide penalty weight option, (LS, uniform)')
    end
else
    switch penalty_weight
        case 'LS'
            [~,~,P1,~,weight_1,weight_2] = weighted_projection_S(xLS,q,[n,p,K],gamma);
            Lambda0 = lambdamax_grouplasso_v2(G,b,m1,[n ,p ,K],P1);
            Lambdacrit_1 = weight_1.*Lambda0;
            Lambdacrit_2 = weight_2.*Lambda0;
            
        case 'uniform'
            Lambdacrit_1 = lambdamax_grouplasso_v2(G,b,m1,[n ,p ,K]);
            Lambdacrit_2 = Lambdacrit_1;
        otherwise
            error('Please provide penalty weight option, (LS, uniform)')
    end
end

end