clear
clc
clf
close all
T = 100000;
ADD = 1000;
NN = [0.1,0.3,1.5];
WW = [100,500,1000];
cnt = 0;
for N = NN
    for WSIZE = WW
        cnt = cnt+1;
        figure(cnt)
% N= 0.1 ;
seed = 154921;
rng(seed)
wn = randn(T,1);
wn = detrend(sqrt(N)*wn/std(wn));
f_ord = [1:30];
for filter_order=f_ord
pd1 = makedist('Stable','alpha',1.8,'beta',0,'gam',1,'delta',0);
an = random(pd1,T+ADD,1);
y = zeros(T+ADD,1);
% y(1:T+ADD-1) = 20;
for tt=11:T+ADD
    y(tt) = y(tt-1)+an(tt)-an(tt-10);
% for kk=0:9
%     y(tt) = y(tt)+an(tt-kk);
% end
end
y = y(ADD+1:end)+wn;
an = an(ADD+1:end);
% WSIZE = 100;

yb = reshape(y,[WSIZE,T/WSIZE]);
anb = reshape(an,[WSIZE,T/WSIZE]);

% filter_order = 5;
Ryy_avg = zeros(2*filter_order-1,1);
Rxy_avg = zeros(2*filter_order-1,1);
for ii=1:T/WSIZE
    [Ryy,lagyy] = xcorr(anb(:,ii),anb(:,ii),filter_order-1,'unbiased');
    [Rxy,lagxy] = xcorr(anb(:,ii),yb(:,ii),filter_order-1,'unbiased');
    Ryy_avg = Ryy_avg+Ryy/(T/WSIZE);
    Rxy_avg = Rxy_avg+Rxy/(T/WSIZE);
end
toep_Ryy_avg = toeplitz(Ryy_avg(find(lagyy==0):end));

% w = toep_Ryy_avg\(Ryx_avg(find(lagyx==0)-filter_order+1:find(lagyx==0)));
% w = toep_Ryy_avg\(Ryx_avg(find(lagxy==0):end));
w = toep_Ryy_avg\flipud(Rxy_avg(find(lagxy==0)-filter_order+1:find(lagxy==0)));
% plot(w)
if filter_order>9
w_opt = [ones(9,1);zeros(filter_order-9,1)];
else
    w_opt = ones(9,1);
    w = [w;zeros(9-filter_order,1)];
end
err = w_opt-w;
WSNR(filter_order) = 10*log10((w_opt'*w_opt)/(err'*err));
end
subplot(2,1,1)
plot(WSNR,'linewidth',1.2)
xlabel('filter order')
ylabel('WSNR (dB)')
title(sprintf('window size:%d, Noise power:%.2f',WSIZE,N))
% subplot(3,1,2)
% plot(w)
% subplot(3,1,3)
% plot(flipud(Ryx_avg(find(lagyx==0)-filter_order+1:find(lagyx==0))))
% plot(Ryx_avg)
subplot(2,1,2)
plot(lagxy,Rxy_avg,'linewidth',1.2)
xlabel('lags')
ylabel('Cross-correlation')
print(figure(cnt),['p2_',int2str(cnt)],'-dpng')
    end
end