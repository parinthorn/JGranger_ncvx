function [snr, trSy, trSe] = compute_ar_SNR(Ak, sigmasq, PARAMETER)
n = PARAMETER(1);
p = PARAMETER(2);

A = [reshape(Ak,n,n*p) ; [eye(n*(p-1)) zeros(n*(p-1),n)]];
B = [eye(n); zeros(n*(p-1),n)]; C = [eye(n) zeros(n,n*(p-1))];
Se = sigmasq*eye(n); Sx = dlyap(A,B*Se*B'); Sy = C*Sx*C'; % sigmasq is the noise variance you generated in VAR system
trSy = trace(Sy);
trSe = trace(Se);
snr = 10*log10(trSy/trSe);

end