%% reviewer response n=100

%% Generate model

%% Configuration
n = 100; % time-series dimension
p = 2;  % ground-truth VAR order
K = 2; % number of models
realization = 100; % number of model's realizations
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
save([outpath,['reviewer_response_model_K',int2str(K),'_p',int2str(p)]],'E')