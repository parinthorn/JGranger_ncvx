clear
clf
close all
clc
T=100000;
N = 0.1;
seed = 154921;
rng(seed)
wn = randn(T,1);
wn = detrend(sqrt(N)*wn/std(wn));

pd1 = makedist('Stable','alpha',1.8,'beta',0,'gam',1,'delta',0);
an = random(pd1,T,1);
dom = -10:0.1:10;
true_pdf = pdf(pd1,dom);
histogram(an,'Normalization','pdf')
hold on
plot(dom,true_pdf,'--r','linewidth',1.4)
hold off
xlim([-10,10])
