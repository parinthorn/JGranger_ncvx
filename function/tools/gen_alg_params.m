function ALG_PARAMETER = gen_alg_params(qnorm, formulation)
% This function can be used to generate algorithm parameters for
% jointvargc.m
% 
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce
switch qnorm
    case 'cvx'
        if strcmp(formulation,'cgn')
            ALG_PARAMETER.MAXITERS=50000;  % max iterations
            ALG_PARAMETER.ABSTOL=1e-7; %1e-7  % absolute tolerance for ADMM convergence detection
            ALG_PARAMETER.RELTOL=1e-5; %1e-5  % relative tolerance
            ALG_PARAMETER.PRINT_RESULT=0;  % print to console [0, 1]
            ALG_PARAMETER.IS_ADAPTIVE =1;  % adaptive rho [0, 1]
            ALG_PARAMETER.rho_init = 1;  % intial rho
            ALG_PARAMETER.epscor = 0.5;  % threshold for spectral rule
            ALG_PARAMETER.Ts = 2;  % update rho every Ts iterations
            ALG_PARAMETER.is_chol = 1;  % Cholesky precomputation or direct solve
            ALG_PARAMETER.multiplier = 2;  % rho multiplier
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