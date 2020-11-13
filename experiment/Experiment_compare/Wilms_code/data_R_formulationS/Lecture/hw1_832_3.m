clear
clc
load('data.mat')
filter_order = 2;
% reference=[reference reference reference reference];
% primary = [primary primary primary primary];
u = reference;
% r = [r r r r r r];
y = primary;
% p = [p p p p p p];
T = length(u);

WSIZE = 1000;
y = reshape(y,[WSIZE,T/WSIZE]);
u = reshape(u,[WSIZE,T/WSIZE]);

% filter_order = 5;
Ryy_avg = zeros(2*filter_order-1,1);
Rxy_avg = zeros(2*filter_order-1,1);
for ii=1:T/WSIZE
    [Ryy,lagyy] = xcorr(u(:,ii),u(:,ii),filter_order-1,'unbiased');
    [Rxy,lagxy] = xcorr(u(:,ii),y(:,ii),filter_order-1,'unbiased');
    Ryy_avg = Ryy_avg+Ryy/(T/WSIZE);
    Rxy_avg = Rxy_avg+Rxy/(T/WSIZE);
end
toep_Ryy_avg = toeplitz(Ryy_avg(find(lagyy==0):end));

% w = toep_Ryy_avg\(Ryx_avg(find(lagyx==0)-filter_order+1:find(lagyx==0)));
% w = toep_Ryy_avg\(Ryx_avg(find(lagxy==0):end));
w_theory = toep_Ryy_avg\flipud(Rxy_avg(find(lagxy==0)-filter_order+1:find(lagxy==0)));

[y,w,bias,J,step_size] = LMS_alg(primary',reference',filter_order,w_theory,1,flipud(Rxy_avg(find(lagxy==0)-filter_order+1:find(lagxy==0))),toep_Ryy_avg);
% plot(bias)
figure(1)
plot(bias)
figure(2)
plot([primary' y'])