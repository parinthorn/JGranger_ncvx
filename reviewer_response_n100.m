
%% Configuration
clear
clc
model_parameter_path = './experiment/model_parameters/';
n = 100; % time-series dimension
p = 2;  % ground-truth VAR order
K = 2; % number of models
realization = 100; % number of model's realizations
common_density = [0.01]; % for p=1, common GC network density [We actually used only 0.1, 0.2]
differential_density = [0.005]; % differential GC network density
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
save([model_parameter_path,['reviewer_response_model_K' int2str(K),'_p',int2str(p)]],'E')
%% plot ground-truth
model_parameter_path = './experiment/model_parameters/';
K=2;
p=2;
load([model_parameter_path,['reviewer_response_model_K' int2str(K),'_p',int2str(p)]],'E')
sample_model = 60;
tt = tiledlayout(2,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
mname = {'common','differential','fused'};
figurepath = './results2plot/figures/';

for kk =1:K
for ii=1:3
    model = E{ii,1,1,sample_model};
    GC = model.GC;
    
    
    [n,~,K] = size(GC);
    commonNZ = ones(n,n)-eye(n);
    diag_ind = 1:n+1:n^2;
    for zz=1:K
        tmp = GC(:,:,zz);
        tmp(diag_ind) = 0;
        commonNZ = commonNZ & (tmp~=0);
    end
    s = [];
    
        tmp = GC(:,:,kk);
        tmp(diag_ind) = 0;
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
        spy(tmp.*(1-commonNZ),4,'sr')
        hold off
%         axis('square')
        xlim([-5 105])
        ylim([-5 105])
        set(gca,'xticklabel',[],'yticklabel',[],'xlabel',[])
        if ii==1
            ylabel(sprintf('model #%d', kk))
        end
        if kk==1
            title(mname{ii})
        end
    end
    hold off
    colormap((1-gray).^1.8)
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
p_est = 1;
K = 2;
load([inpath,'reviewer_response_model_K',int2str(K),'_p',int2str(p_true)]) % struct E
m= size(E,4);
realz = m;
GridSize = 30;
mname = {'1','5'};
ii=2;
parameter.cvx.varorder = p_est;
parameter.cvx.formulation = 'dgn'; % cgn, dgn, fgn
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
for jj=1:realz
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
    save([outpath,'reviewer_response_estim_DGN_cvx_T2500_',mname{ii},'percent','_lag',int2str(p_est),'_K',int2str(K),'_',int2str(jj)],'M')

end

