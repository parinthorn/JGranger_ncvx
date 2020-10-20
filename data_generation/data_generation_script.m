clear
clc
clf
close all
outpath = './data_compare/';

n = 20;
p = 1;
K = 50;
realization = 100;

common_density = [0.01;0.05;0.1;0.2];
differential_density = [0.01;0.05];
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
save([outpath,['model_allcase_K',int2str(K)]],'E')
