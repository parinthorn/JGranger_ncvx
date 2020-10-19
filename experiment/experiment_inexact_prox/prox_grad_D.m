function [x,history] = prox_grad_D(G,b,PInd,a1,a2,pp,qq,PARAMETERS,varargin)
IT_MAX = 200000;
TOL = 1e-6;
objTOL = 1e-6;
ALLPRINT = 1;FREQ = 1;
n = PARAMETERS(1);
p = PARAMETERS(2);
K = PARAMETERS(3);
obj = @(u) 0.5*(norm(G*u-b))^2+a1*normpq(u(PInd),pp,qq,p)+a2*normpq(u(PInd),pp,qq,p*K);
if isempty(varargin)
    x0 = zeros(n^2*p*K,1);
else
    x0 = varargin{1};
end
xk = x0;
GtG = G'*G;
Gtb = G'*b;
GtGxk_Gtb = GtG*x0-Gtb;
STEP_SIZE = 1/eigs(GtG,1);
kk=1;
while kk<IT_MAX
    if a2==0
        x_Px = prox_pq_Px(xk-STEP_SIZE*(GtGxk_Gtb),PInd,STEP_SIZE,a1,pp,qq,PARAMETERS,p);
        x_Px_eff = prox_pq_Px_efficient(xk-STEP_SIZE*(GtGxk_Gtb),PInd,STEP_SIZE,a1,pp,qq,PARAMETERS,p);
        x_num =  prox_Px_numerical(xk-STEP_SIZE*(GtGxk_Gtb),PInd,STEP_SIZE*a1,STEP_SIZE*a2,qq,PARAMETERS,'optimal');
        xkp1 = x_num;
        fprintf('exact-numerical:%d, numerical-exact:%d \n',length(setdiff(find(x_Px~=0),find(x_num~=0))),length(setdiff(find(x_num~=0),find(x_Px~=0))))
    else
        xkp1 =  prox_Px_numerical(xk-STEP_SIZE*(GtGxk_Gtb),PInd,a1,a2,qq,PARAMETERS,'optimal');
    end
    
    history.obj(kk,1) = obj(xkp1);
    if ((mod(kk,FREQ)==0) && ALLPRINT)
        fprintf('%3d\t%10.4f\t%10.4f\n',kk,history.obj(kk,1),norm(xkp1-xk)/norm(xk))
    end
    
    if norm(xkp1-xk) < TOL*norm(xk)
        x = xkp1;
        history.FLAG=1;
        break
    end
    
    xk = xkp1;
    GtGxk_Gtb = GtG*xk-Gtb;
    kk = kk+1;
end


end