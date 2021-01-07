function x = log_stirling_approx(n)
% based on Ramanujan approximation
x = n*log(n)-n+log(n*(1+4*n*(1+2*n)))/6+log(pi)/2;
end