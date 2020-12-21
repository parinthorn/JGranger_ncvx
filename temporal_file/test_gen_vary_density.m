clear
% clc
clf
close all
outpath = './data_compare/';

n = 20;
p = 1;
K = 5;
realization = 20;
% common_density = [0.02;0.04;0.05;0.06;0.08]; % for p=1
common_density = [0.05,0.15,0.25];
% common_density = [0.01;0.05]; % for p=3
% differential_density = 0.1-common_density;
differential_density = 0.3-common_density;
model = {'common','differential','similar'};
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
save([outpath,['vary_CD_model',int2str(K),'_p',int2str(p)]],'E')
