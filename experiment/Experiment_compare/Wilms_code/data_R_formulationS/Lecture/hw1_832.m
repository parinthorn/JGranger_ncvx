%% 1
clear
clc

seed = 154921; % 6270154921 is my ID
rng(seed)
T = 100;
tmp =randn(T,1);
v = detrend(tmp)/(std(tmp)); % zero mean unit variance white noise
% MA process
x = zeros(T,1);
M_all = [2,5,10];
figure(1)
for m=1:length(M_all)
ax(m) = subplot(length(M_all),1,m);
M = M_all(m);
syms s
f = 1/(1+0.75*s+0.25*s^2);
Taylor_approx = taylor(f, 'Order', M+10);
a = fliplr(sym2poly(Taylor_approx));
b=a;
b(1) = [];
u = zeros(T,1);
for tt=M+1:T
    x(tt) = v(tt)+0.75*v(tt-1)+0.25*v(tt-2);
    u(tt) = -b(1:M)*flipud(u(tt-M:tt-1))+v(tt);
end
plot([M+1:T],[v(M+1:end) x(M+1:end) u(M+1:end)],'linewidth',1.2)
rel_err(m) = norm(x(M+1:end)-u(M+1:end))/norm(x(M+1:end));
title(sprintf('relative error:%f',rel_err(m)))
legend('white noise','MA',sprintf('AR order:%d',M),'Location','southeast')
end
for ii=2:length(ax)
linkaxes([ax(ii) ax(ii-1)],'x')
end
% figure(2)
% plot(M_all,rel_err) % plot to see performance at each M
%%
M_set=10;
f = 1/(1+0.75*s+0.25*s^2);
Taylor_approx = taylor(f, 'Order', M_set+10);
a = fliplr(sym2poly(Taylor_approx));

MA_num = [1 0.75 0.25];
AR_den = a(1:M_set+1);
MA_filter = tf(MA_num,1,1);
AR_filter = tf(1,AR_den,1);
% figure(3)
figure(3)
pzmap(MA_filter)
hold on
pzmap(AR_filter)
hold off
axis('square')
fvtool(1,MA_num,AR_den,1)
% freqz(1,MA_num)
% hold on
% freqz(AR_den,1)
% hold off
%% 2
clear
clc
clf
close all
T = 10000;
N=0.1;
pd1 = makedist('Stable','alpha',1.8,'beta',0,'gam',1,'delta',0);
an = random(pd1,T,1);
wn = randn(T,1);
wn = sqrt(N)*(wn-mean(wn))/std(wn);