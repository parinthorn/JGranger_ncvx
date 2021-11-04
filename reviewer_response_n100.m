
%% Configuration
clear
clc
model_parameter_path = './experiment/model_parameters/';
n = 20; % time-series dimension
p = 2;  % ground-truth VAR order
K = 3; % number of models
realization = 1; % number of model's realizations
common_density = [0.1]; % for p=1, common GC network density [We actually used only 0.1, 0.2]
differential_density = [0.05]; % differential GC network density
model = {'common','differential','similar'}; % type of similarity
mname = {'C','D','S'};
cnt = 0;
E=cell(length(model),length(common_density),length(differential_density),realization);
for m=1:length(model)
    for d=1:length(common_density)
        opts.common_density = common_density(d);
        for diff_d =1:length(differential_density)
            opts.differential_density = differential_density(diff_d);
            opts.type = model{m};

            for b=1:realization %number of [C,S,D] VAR model generated
                if strcmp(mname{m},'D')
                    E{m,d,diff_d,b} = gen_multi_VAR([n,p,K],opts,E{1,d,diff_d,b}.A); % look for C type model [code 1] to generate D type
                else
                    E{m,d,diff_d,b} = gen_multi_VAR([n,p,K],opts);
                end
            end
        end
    end
end
% save([model_parameter_path,['reviewer_response_model_K' int2str(K),'_p',int2str(p)]],'E')
%% plot ground-truth
% model_parameter_path = './experiment/model_parameters/';
K=3;
p=2;
% load([model_parameter_path,['model_K' int2str(5),'_p',int2str(p)]],'E')
sample_model = 1;
tt = tiledlayout(K,3);
tt.TileSpacing = 'None';
tt.Padding = 'None';
mname = {'common','differential','fused'};
figurepath = './results2plot/figures/';
for ii=1:3
for kk =1:K

    model = E{ii,1,1,sample_model};
    GC = model.GC(:,:,1:K);
    A = model.A;
    A = squeeze(sqrt(sum(A.^2,3)));

    
    
    [n,~,K] = size(A);
    commonNZ = ones(n,n)-eye(n);
    diag_ind = 1:n+1:n^2;
    for zz=1:K
        tmp = A(:,:,zz);
        commonNZ = commonNZ & (tmp~=0);
    end
    s = [];
    
        tmp = A(:,:,kk);
        nexttile;
%         spy(tmp, 10, 'k')
%         spyc(tmp)
%         sA = tmp;
%         indx=find(sA);
%         [Nx Ny]=size(sA);
%         sA=full(sA(indx));
%         ns = length(indx);
%         [ix iy]=ind2sub([Nx Ny],indx);
%         xx = GC(:);
%         imap = sA;
%         scatter(iy,ix,[],imap,'Marker','s', 'MarkerFaceColor', 'flat')
imagesc(tmp)
set(gca,'box','on', 'BoxStyle','full')

% H = gca;
% H.LineWidth = 1.5;
% set(gca,'box','off')
        
        
        hold on
        spy(tmp.*(1-commonNZ),30,'r')
        hold off
%         axis('square')
%         xlim([-5 105])
%         ylim([-5 105])
        set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
        if ii==1
                        title(sprintf('model #%d', kk))

%             title()
        end
        if kk==1
%             title()
ylabel(mname{ii})
        end
    end
    hold off
    colormap((1-gray))
%             colorbar;

    
end
pp = get(0, 'Screensize');
pp(3) = pp(3)*1/3;
pp(4) = pp(4)*1/3;
pp = pp+200;
set(gcf, 'Position', pp);
set(findall(gcf,'-property','FontSize'),'FontSize',28)
print([figurepath,'reviewer_response_sampleGC'],'-painters','-depsc','-r300')
print([figurepath,'reviewer_response_sampleGC'],'-painters','-dpng','-r300')

%% Estimation
clear
clc
inpath = './experiment/model_parameters/';
outpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx\';
mkdir(outpath)

T = 2500;
p_true = 2;
p_est = 2;
K = 2;
load([inpath,'reviewer_response_model_K',int2str(K),'_p',int2str(p_true)]) % struct E
m= size(E,4);
realz = m;
GridSize = 30;
mname = {'1','5'};
ii=2;
parameter.cvx.varorder = p_est;
parameter.cvx.formulation = 'cgn'; % cgn, dgn, fgn
parameter.cvx.penalty_weight = 'LS'; % LS, uniform
parameter.cvx.GridSize = GridSize;
parameter.cvx.data_concat = 0;
parameter.cvx.noisecov = 'full'; % 'full', 'diag', 'identity'
parameter.cvx.qnorm = 'cvx';  % cvx, ncvx



ALG_PARAMETER.cvx = gen_alg_params(parameter.cvx.qnorm, parameter.cvx.formulation);
ALG_PARAMETER.cvx.IS_ADAPTIVE = 1;
ALG_PARAMETER.cvx.PRINT_RESULT = 0;
% ALG_PARAMETER.cvx.Ts = 1;

parameter.ncvx = parameter.cvx;
parameter.ncvx.qnorm = 'ncvx';
ALG_PARAMETER.ncvx = gen_alg_params(parameter.ncvx.qnorm, parameter.ncvx.formulation);
% ALG_PARAMETER.ncvx.PRINT_RESULT = 1;
%% Non-CVX
for jj=10:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
    M = jointvargc(y,parameter.ncvx,ALG_PARAMETER.ncvx);
    
    save([outpath,'reviewer_response_estim_CGN_T2500_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')
end
%% CVX
inference_time = zeros(realz,1);
for jj=1:realz
    % generate data from given seed
    model = E{2,jj}; % type D
    y = sim_VAR(model.A,T,1,model.seed,0);
    
    t1 = tic;
    M = jointvargc(y,parameter.cvx,ALG_PARAMETER.cvx);
    inference_time(jj) = toc(t1);
    fprintf("model: %d, time:%.2f", jj, inference_time(jj))
    save([outpath,'reviewer_response_estim_CGN_cvx_T2500_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

end

%% cvx-CGN evaluation
clear
% clc
inpath = './experiment/model_parameters/';
resultpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx/';
performance_path = './results2plot/';
mname = {'5'};
name_list = {'bic_lasso','bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
realization = 10;
load([inpath,'reviewer_response_model_K2_p2'])
toggle_list = {'common'};
p_est = 2;
K=2;
for ii=1:length(mname)
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{2,1,1,jj};
        fname = [resultpath,'reviewer_response_estim_CGN_cvx_T2500_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)];
        load(fname)
        model_acc = performance_eval(M,GTmodel);
        ALL_RESULT(ii,jj).model_acc = model_acc;
        for kk=1:length(name_list)
            R.index(ii,jj).(name_list{kk}) = M.index.(name_list{kk});
        end
        for tt = 1:length(toggle_list)
            toggle = toggle_list{tt};
            R.(toggle).F1(ii,jj) =model_acc(M.index.eBIC).(toggle).F1;
            R.(toggle).MCC(ii,jj) =model_acc(M.index.eBIC).(toggle).MCC;
            R.(toggle).TPR(ii,jj) =model_acc(M.index.eBIC).(toggle).TPR;
            R.(toggle).FPR(ii,jj) =model_acc(M.index.eBIC).(toggle).FPR;
            R.(toggle).ACC(ii,jj) =model_acc(M.index.eBIC).(toggle).ACC;
        end
        R.bias(ii,jj) =model_acc(M.index.eBIC).bias;
        
    end
end
% save([performance_path,'reviewer_response_CVX_CGN_result'],'R')
% save([performance_path,'reviewer_response_CVX_CGN_ALL_RESULT'],'ALL_RESULT')

%% CGN evaluation
clear
clc
inpath = './experiment/model_parameters/';
resultpath = 'G:\My Drive\0FROM_SHARED_DRIVE\THESIS\experiment_compare_cvx/';
performance_path = './results2plot/';
mname = {'5'};
name_list = {'bic_lasso','bic','aic','aicc','eBIC','GIC_2','GIC_3','GIC_4','GIC_5','GIC_6'};
realization = 10;
load([inpath,'reviewer_response_model_K2_p2'])
toggle_list = {'common'};
p_est = 2;
K=2;
for ii=1:length(mname)
    for jj=1:realization
        fprintf('(%d,%d)\n',ii,jj)
        GTmodel = E{2,1,1,jj};
        fname = [resultpath,'reviewer_response_estim_CGN_T2500_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)];
        load(fname)
        model_acc = performance_eval(M,GTmodel);
        ALL_RESULT(ii,jj).model_acc = model_acc;
        for kk=1:length(name_list)
            R.index(ii,jj).(name_list{kk}) = M.index.(name_list{kk});
        end
        for tt = 1:length(toggle_list)
            toggle = toggle_list{tt};
            R.(toggle).F1(ii,jj) =model_acc(M.index.eBIC).(toggle).F1;
            R.(toggle).MCC(ii,jj) =model_acc(M.index.eBIC).(toggle).MCC;
            R.(toggle).TPR(ii,jj) =model_acc(M.index.eBIC).(toggle).TPR;
            R.(toggle).FPR(ii,jj) =model_acc(M.index.eBIC).(toggle).FPR;
            R.(toggle).ACC(ii,jj) =model_acc(M.index.eBIC).(toggle).ACC;
        end
        R.bias(ii,jj) =model_acc(M.index.eBIC).bias;
        
    end
end
% save([performance_path,'reviewer_response_CGN_result'],'R')
% save([performance_path,'reviewer_response_CGN_ALL_RESULT'],'ALL_RESULT')
%%
load([performance_path,'reviewer_response_CVX_CGN_result'],'R')
load([performance_path,'reviewer_response_CVX_CGN_ALL_RESULT'],'ALL_RESULT')
%%
load([performance_path,'reviewer_response_CGN_result'],'R')
load([performance_path,'reviewer_response_CGN_ALL_RESULT'],'ALL_RESULT')
