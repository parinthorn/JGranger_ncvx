clear
clc
t = 0:0.05:1;
f=1;
y = sin(2*pi*f*t);

stem(t,y)


fy = fft(y);

plot(abs(fy))

pad_fy = [fy(1:11) zeros(1,20) fy(12:end)];

yhat = ifft(pad_fy)*length(pad_fy)/(length(fy));

stem(t,y)
hold on
stem(linspace(0,1,length(pad_fy)),yhat)
y_I = interp(y,10);
stem(linspace(0,1,length(y_I)),y_I)


hold off

