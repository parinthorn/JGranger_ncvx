function ALG_PARAMETER = gen_alg_params(qnorm, formulation)

switch qnorm
    case 'cvx'
        if strcmp(formulation,'cgn')
            ALG_PARAMETER.MAXITERS=50000;
            ALG_PARAMETER.ABSTOL=1e-7; %1e-7
            ALG_PARAMETER.RELTOL=1e-5; %1e-5
            ALG_PARAMETER.PRINT_RESULT=0;
            ALG_PARAMETER.IS_ADAPTIVE =1;
            ALG_PARAMETER.rho_init = 1;
            ALG_PARAMETER.epscor = 0.5;
            ALG_PARAMETER.Ts = 2;
            ALG_PARAMETER.is_chol = 1;
            ALG_PARAMETER.multiplier = 2;
        elseif strcmp(formulation,'dgn')
            ALG_PARAMETER.MAXITERS=50000;
            ALG_PARAMETER.ABSTOL=1e-7; %1e-7
            ALG_PARAMETER.RELTOL=1e-5; %1e-5
            ALG_PARAMETER.PRINT_RESULT=0;
            ALG_PARAMETER.IS_ADAPTIVE =1;
            ALG_PARAMETER.rho_init = 1;
            ALG_PARAMETER.epscor = 0.5;
            ALG_PARAMETER.Ts = 2;
            ALG_PARAMETER.is_chol = 1;
            ALG_PARAMETER.multiplier = 2;
            
        elseif strcmp(formulation,'fgn')
            ALG_PARAMETER.MAXITERS=50000;
            ALG_PARAMETER.ABSTOL=1e-7; %1e-7
            ALG_PARAMETER.RELTOL=1e-5; %1e-5
            ALG_PARAMETER.PRINT_RESULT=0;
            ALG_PARAMETER.IS_ADAPTIVE =1;
            ALG_PARAMETER.rho_init = 1;
            ALG_PARAMETER.epscor = 0.5;
            ALG_PARAMETER.Ts = 2;
            ALG_PARAMETER.is_chol = 1;
            ALG_PARAMETER.multiplier = 2;
        end
        
    case 'ncvx'
        if strcmp(formulation,'cgn')
            ALG_PARAMETER.MAXITERS=50000;
            ALG_PARAMETER.ABSTOL=1e-7; %1e-7
            ALG_PARAMETER.RELTOL=1e-5; %1e-5
            ALG_PARAMETER.PRINT_RESULT=0;
            ALG_PARAMETER.IS_ADAPTIVE =1;
            ALG_PARAMETER.rho_init = 1;
            ALG_PARAMETER.epscor = 0.1;
            ALG_PARAMETER.Ts = 50;
            ALG_PARAMETER.is_chol = 1;
            ALG_PARAMETER.multiplier = 2;
        elseif strcmp(formulation,'dgn')
            ALG_PARAMETER.MAXITERS=50000;
            ALG_PARAMETER.ABSTOL=1e-7; %1e-7
            ALG_PARAMETER.RELTOL=1e-5; %1e-5
            ALG_PARAMETER.PRINT_RESULT=0;
            ALG_PARAMETER.IS_ADAPTIVE =1;
            ALG_PARAMETER.rho_init = 1;
            ALG_PARAMETER.epscor = 0.1;
            ALG_PARAMETER.Ts = 200;
            ALG_PARAMETER.is_chol = 1;
            ALG_PARAMETER.multiplier = 2;
        elseif strcmp(formulation,'fgn')
            ALG_PARAMETER.MAXITERS=50000;
            ALG_PARAMETER.ABSTOL=1e-7; %1e-7
            ALG_PARAMETER.RELTOL=1e-5; %1e-5
            ALG_PARAMETER.PRINT_RESULT=0;
            ALG_PARAMETER.IS_ADAPTIVE =1;
            ALG_PARAMETER.rho_init = 1;
            ALG_PARAMETER.epscor = 0.1;
            ALG_PARAMETER.Ts = 100;
            ALG_PARAMETER.is_chol = 1;
            ALG_PARAMETER.multiplier = 2;
        end
end
ALG_PARAMETER.FREQ_PRINT = 100;
end