

function [ind_nz,A,y] = gen_sparseAR_similar(B,ind_nz,noise_var,Num,density)
% gen_sparseAR generates a sparse vector autoregressive model
% 
% [ind_nz,phis,y] =  gen_sparseAR_similar(A,ind_nz,noise_var,Num,density)
%
% y(t) = A1*y(t-1) + A2*y(t-2) + ... + Ap*y(t-p) + u(t)
%
% 'A' represents AR coefficients A1,A2,...,Ap and is stored as a p-dimensional array
% 
% 'A' has SIMILAR structure to B except a few diagonal elements controlled by 'density'
% 
% The input arguments are 
% 'B': the original AR coefficients
% 'n': dimension of AR coefficient matrices
% 'p': AR model order
% 'noise_var': variance of u(t) (noise)
% 'density': the fraction of different entries in A from B
% 'Num': number of data points in time series
% 
% The AR coefficients are sparse with a common sparsity pattern. The
% indices of nonzero entries are saved in 'ind_nz'.
% 
% 'y' is a time series generated from the model and has size n x Num
%  y = [y(1) y(2) ... y(Num)]
% 
% if p = 0, 'y' is simply a random variable. In this case, A is the
% covariance matrix of u with sparse inverse.
% 'S' is the given sparse matrix and the AR model will have the same
% sparsity pattern as S

%% Static case
n = size(B,1); p = size(B,3);

if (p==0),
%     S = sparse(2*eye(n)+sign(sprandsym(n,density)));  
    [i,j]=find(S);
    S = S+sparse(ceil(max(0,-min(eig(S))))*eye(n));
    A = S\eye(n); % covariance matrix with sparse inverse
    R = chol(phi);   
    y = R'*randn(n,Num); % y reduces to a random variable with covariance 'phi'
    ind_nz = sub2ind([n n],i,j);
    figure;plot_spy(ind_nz,n,'image');
    title('correct sparsity');
    return;
end
    
%% Randomize AR coefficients   


MAX_EIG = 1;    
diag_ind = (1:n+1:n^2)';
% k = length(diag_ind);
% diag_ind3D = kron(n^2*(0:p-1)',ones(k,1))+kron(ones(p,1),diag_ind);

A = zeros(n,n,p);
S = sprand(n,n,density);
S(diag_ind) = 0;
S(ind_nz) = 0; 
ii = 0;
while MAX_EIG,
    ii = ii+1;
    for k=1:p,
        A(:,:,k) = B(:,:,k)+ 0.2*sprandn(S); % add some off-diagonal entries
    end
    AA =[];
    for k=1:p,
        AA = [AA A(:,:,k)];
    end
    AA = sparse([AA ; [speye(n*(p-1)) zeros(n*(p-1),n)]]);
    
    if max(abs(eigs(AA))) < 1
        MAX_EIG = 0;
    end
end
abs(eigs(AA))
ii

%% the sparsity pattern of A1,A2,...,Ap
S1 = S; S1(ind_nz) = 1; 
ind_nz = find(S1);
figure;
subplot(1,2,1);
plot_spy(ind_nz,n,'image');title('sparsity of AR coefficients');

%% Generate time series
% noise_var = 0;
u = sqrt(noise_var)*randn(n,Num);

y = zeros(n,Num);
y(:,1:p) = randn(n,p);
matA = reshape(A,n,n*p);
for k=p+1:Num,
    y(:,k) = matA*reshape(y(:,k-1:-1:k-p),n*p,1);
    y(:,k) = y(:,k) + u(:,k);
end

aux = max(abs(y),[],2); % find the largest peak in y
[tmp,I] = max(aux);
subplot(1,2,2);
plot(y(I,:)); % plot the entry of y that has the largest peak
title(['component ',num2str(I),'of time series']);xlabel('t');
y = detrend(y','constant')';

close all; % plot and close
