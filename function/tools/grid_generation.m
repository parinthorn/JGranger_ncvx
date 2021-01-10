function [Lambda_1,Lambda_2,opt] = grid_generation(G,b,GridSize,PARAMETER,q,toggle)
Lambda = logspace(-6,0,GridSize);
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);
m1 = PARAMETER.dim(4);
m2 = PARAMETER.dim(5);
opt.detail = toggle;
switch toggle
    case 'static'
        [P,~] = offdiagJSS(n,p,K);
        Lmax_1 = lambdamax_grouplasso_v2(G,b,m1,[n ,p ,K]);
        Lmax_2 = lambdamax_grouplasso_v2(G,b,m2,[n ,p ,K]);
        Lambda_1 = Lambda*Lmax_1;
        Lambda_2 = Lambda*Lmax_2;
        opt.L1 = P;
        opt.L2 = P;
    case 'adaptive_P'
        xLS = G\b;
        [P_m1,P_m2,P1_m1,P1_m2] = weighted_projection(xLS,q,PARAMETER);
        Lmax_1 = lambdamax_grouplasso_v2(G,b,m1,[n ,p ,K],P1_m1);
        Lmax_2 = lambdamax_grouplasso_v2(G,b,m2,[n ,p ,K],P1_m2);
        Lambda_1 = Lambda*Lmax_1;
        Lambda_2 = Lambda*Lmax_2;
        opt.L1 = P_m1;
        opt.L2 = P_m2;
    case 'adaptive_L'
        xLS = G\b;
        [~,~,P1_m1,P1_m2,P] = weighted_projection(xLS,q,PARAMETER);
        Lambda0_1 = lambdamax_grouplasso_v2(G,b,m1,[n ,p ,K],P1_m1);
        Lambda0_2 = lambdamax_grouplasso_v2(G,b,m2,[n ,p ,K],P1_m2);
        
        weight_1 = max(P1_m1,[],2);
        weight_2 = max(P1_m2,[],2);
        
        Lambda_1 = weight_1.*Lambda.*Lambda0_1;
        Lambda_2 = weight_2.*Lambda.*Lambda0_2;
        
        opt.L1 = P;
        opt.L2 = P;
end
end