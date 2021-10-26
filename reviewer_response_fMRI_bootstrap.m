%% Data concatenation [Run this first]
clear
clc
selected_TDC = {'0010023';'0010024';'0010070';'0010088';'0010092';'0010112';'0010122';'0010123';'0010128';'1000804';'1435954';'3163200';'3243657';'3845761';'4079254';'8692452';'8834383';'9750701'};
selected_ADHD_C = {'0010013';'0010017';'0010019';'0010022';'0010025';'0010037';'0010042';'0010048';'0010050';'0010086';'0010109';'1187766';'1283494';'1918630';'1992284';'2054438';'2107638';'2497695'};
y_TDC = concat_real_data(selected_TDC,116,'nyu',0);
y_ADHD_C = concat_real_data(selected_ADHD_C,116,'nyu',0); % load data in format (n,T,K)
y_TDC = y_TDC-mean(y_TDC,2);
y_ADHD_C = y_ADHD_C-mean(y_ADHD_C,2);
K = size(y_TDC,3);
y_total = cat(3,y_TDC,y_ADHD_C);
y_total = y_total-mean(y_total,2); % detrend in time axis


parameter.cvx.varorder = 1;
parameter.cvx.formulation = 'dgn'; % cgn, dgn, fgn
parameter.cvx.penalty_weight = 'LS'; % LS, uniform
parameter.cvx.GridSize = 30;
parameter.cvx.data_concat = 0;
parameter.cvx.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.cvx.qnorm = 'cvx';  % cvx, ncvx
ALG_PARAMETER.cvx = gen_alg_params(parameter.cvx.qnorm, parameter.cvx.formulation);

patient_sample = 18;
T_per_patient = 172;
n = 116;
p_est = parameter.cvx.varorder;


eff_T = (T_per_patient-p_est)*patient_sample;
H_TDC = zeros(n*p_est,eff_T);
Y_TDC = zeros(n,eff_T);
H_ADHD = zeros(n*p_est,eff_T);
Y_ADHD = zeros(n,eff_T);


for kk=[1:18]
    % TDC
    [H_TDC(:,(T_per_patient-p_est)*(kk-1)+1:(T_per_patient-p_est)*(kk)), ...
     Y_TDC(:,(T_per_patient-p_est)*(kk-1)+1:(T_per_patient-p_est)*(kk))] = ...
        H_gen(y_total(:,:,kk),p_est);
end
for kk=[1:18]
    % ADHD_C
        [H_ADHD(:,(T_per_patient-p_est)*(kk-1)+1:(T_per_patient-p_est)*(kk)), ...
         Y_ADHD(:,(T_per_patient-p_est)*(kk-1)+1:(T_per_patient-p_est)*(kk))] = ...
            H_gen(y_total(:,:,kk+18),p_est);
end

H(:,:,1) = H_TDC;
H(:,:,2) = H_ADHD;

Y(:,:,1) = Y_TDC;
Y(:,:,2) = Y_ADHD;

% bootstrapping

block_length = 162;
T_eff = 3078;

num_blocks = 19;
repetitions = 5;

for ii=0:T_eff-block_length-1
    bootstrap_sample{ii+1} = (1:block_length) + ii;
%     bootstrap_index = randi([1 100],1,5);
end
s = RandStream('mlfg6331_64');

for tt=1:repetitions

bootstrap_index = randsample(s, length(bootstrap_sample),num_blocks,false);

Y_bootstrap = zeros(n, T_eff, 2);
H_bootstrap = zeros(n*p_est, T_eff, 2);
for bb=1:length(bootstrap_index)
    
    
    disp('Concatenating H, Y matrix')
    H_bootstrap(:,(block_length)*(bb-1)+1:(block_length)*(bb), :) = H(:,bootstrap_sample{bootstrap_index(bb)},:);
    
    Y_bootstrap(:,(block_length)*(bb-1)+1:(block_length)*(bb), :) = Y(:,bootstrap_sample{bootstrap_index(bb)},:);
    
end


parameter.cvx.is_YH = 1;
parameter.cvx.Y= Y_bootstrap;
parameter.cvx.H=H_bootstrap;
parameter.cvx.eff_T=eff_T;

M_bootstrap = jointvargc(Y_bootstrap,parameter.cvx,ALG_PARAMETER.cvx);

save(['G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_D2K_bootstrap_', int2str(tt)],'M') % verified for reproduce
end

% 
% 
% ALG_PARAMETER = gen_alg_params(parameter.qnorm, parameter.formulation);
% M = jointvargc(y_total,parameter,ALG_PARAMETER); % data with dimension (n,T*K,2), K is # subjects in each TDC, ADHD
% save('G:\My Drive\0FROM_SHARED_DRIVE\THESIS\Real_data\experiment_real_data_result\estim_D2K','M') % verified for reproduce