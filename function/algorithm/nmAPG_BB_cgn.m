function [x,history]= nmAPG_BB_cgn(G,b,lambda,PARAMETER,varargin)
IT_MAX = 100000;
TOL = 1e-6;
objTOL = 1e-8;
ALLPRINT = 0;FREQ = 1;

PInd = PARAMETER.proj_ind;
pp = PARAMETER.p;
qq = PARAMETER.q;
n = PARAMETER.dim(1);
p = PARAMETER.dim(2);
K = PARAMETER.dim(3);
d = PARAMETER.delta;
eta = PARAMETER.eta;
rho = PARAMETER.rho;

IS_LINESEARCH = PARAMETER.IS_LINESEARCH;
if isempty(varargin)
    x0 = zeros(n^2*p*K,1);
else
    x0 = varargin{1};
end
% obj = @(u) 0.5*(norm(G*u-b))^2+a1*normpq(u(PInd),pp,qq,p)+a2*normpq(u(PInd),pp,qq,p*K);
obj = @(u) 0.5*(norm(G*u-b))^2+normpq_adaptive(u(PInd),pp,qq,p*K,lambda);
objx0 = obj(x0);
% APG INITIALIZATION
k=1;
xkm1 = x0;
zk = xkm1;
xk = xkm1;
tk = 1;
tkm1 = 0;
ck = obj(xk);
qk = 1;
ykm1 = zeros(n^2*p*K,1);
GtG = G'*G;
STEP_SIZE = 1/eigs(GtG,1);
Gtb = G'*b;
axITR = 0;
ayITR = 0;
if ALLPRINT
    fprintf('%3s\t%10s\t%10s\n','iter', 'objective','step size');
end
t_iteration = tic;
while k<IT_MAX
    yk = xk+(tkm1/tk)*(zk-xk)+((tkm1-1)/tk)*(xk-xkm1);
    GtGyk_Gtb = GtG*yk-Gtb;
    
    if (IS_LINESEARCH) && (k>2)
        sk = yk-ykm1;
        rk = GtG*(sk);
        ay = norm(sk,2)^2/(sk'*rk);

%         zkp1 = prox_Px_numerical(yk-ay*(GtGyk_Gtb),a1*ay,a2*ay,PARAMETER);
        zkp1 = prox_pq_Px_weighted(yk-ay*(GtGyk_Gtb),PInd,lambda,ay,pp,qq,[n,p,K],p*K);
        while (obj(zkp1)>=obj(yk)-d*norm(zkp1-yk)^2 )&& (obj(zkp1)>=ck-d*norm(zkp1-yk)^2 )
            ayITR = ayITR+1;
            ay = ay*rho;
            if ay<=STEP_SIZE
                ay = STEP_SIZE;

%                 zkp1 = prox_Px_numerical(yk-ay*(GtGyk_Gtb),a1*ay,a2*ay,PARAMETER);
                zkp1 = prox_pq_Px_weighted(yk-ay*(GtGyk_Gtb),PInd,lambda,ay,pp,qq,[n,p,K],p*K);
                break
            end

%             zkp1 = prox_Px_numerical(yk-ay*(GtGyk_Gtb),a1*ay,a2*ay,PARAMETER);
            zkp1 = prox_pq_Px_weighted(yk-ay*(GtGyk_Gtb),PInd,lambda,ay,pp,qq,[n,p,K],p*K);
        end
    else
        ay = STEP_SIZE;

%         zkp1 = prox_Px_numerical(yk-ay*(GtGyk_Gtb),a1*ay,a2*ay,PARAMETER);
        zkp1 = prox_pq_Px_weighted(yk-ay*(GtGyk_Gtb),PInd,lambda,ay,pp,qq,[n,p,K],p*K);
    end
    
    
    if (obj(zkp1)<=ck-d*norm(zkp1-yk)^2 )
        xkp1 = zkp1;
    else
        GtGxk_Gtb = GtG*xk-G'*b;
        if (IS_LINESEARCH) && (k>2)
            sk = xk-ykm1;
            rk = GtG*(xk-ykm1);
            ax = norm(sk,2)^2/(sk'*rk);
%             vkp1 = prox_Px_numerical(xk-ax*(GtGxk_Gtb),a1*ax,a2*ax,PARAMETER);
            vkp1 = prox_pq_Px_weighted(xk-ax*(GtGxk_Gtb),PInd,lambda,ax,pp,qq,[n,p,K],p*K);
            while obj(vkp1) >= ck-d*norm(vkp1-xk,2)^2
                axITR = axITR+1;
                ax = ax*rho;
                if ax<= STEP_SIZE
                    ax = STEP_SIZE;
%                     vkp1 = prox_Px_numerical(xk-ax*(GtGxk_Gtb),a1*ax,a2*ax,PARAMETER);
                    vkp1 = prox_pq_Px_weighted(xk-ax*(GtGxk_Gtb),PInd,lambda,ax,pp,qq,[n,p,K],p*K);
                    break
                end
%                 vkp1 = prox_Px_numerical(xk-ax*(GtGxk_Gtb),a1*ax,a2*ax,PARAMETER);
                vkp1 = prox_pq_Px_weighted(xk-ax*(GtGxk_Gtb),PInd,lambda,ax,pp,qq,[n,p,K],p*K);
            end
        else
            ax = STEP_SIZE;
%             vkp1 = prox_Px_numerical(xk-ax*(GtGxk_Gtb),a1*ax,a2*ax,PARAMETER);
            vkp1 = prox_pq_Px_weighted(xk-ax*(GtGxk_Gtb),PInd,lambda,ax,pp,qq,[n,p,K],p*K);
        end
        if obj(zkp1)<=obj(vkp1)
            xkp1 = zkp1;
        else
            xkp1 = vkp1;
        end
    end
    
    
    objval = obj(xkp1);
    tkp1 = 0.5*(sqrt(4*tk^2+1)+1);
    qkp1 = eta*qk+1;
    ckp1 = (eta*qk*ck+objval)/qkp1;
    ykm1 = yk;
    history.t(k,1) = toc(t_iteration);
    history.objval(k,1) = objval;
    history.c(k,1) = ckp1;
    if k==1
        history.reldiff_norm(k)  = 1;
    else
        history.reldiff_norm(k)  = norm(xkp1-xk)/norm(xk);
    end
    if (norm(xkp1-xk)<(TOL*norm(xk))) % || ((k>1) && (abs(objval-history.objval(k-1))<objTOL*objval))
        x = sparse(xkp1);
        history.flag = 0;
        history.objval = [objx0;history.objval];
        break
    else
        x = sparse(xkp1);
        history.flag = -1; % not converge yet
    end
    
    
    k = k+1;
    if ((mod(k,FREQ)==0) && ALLPRINT)
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\n',k,objval,ckp1,norm(xkp1-xk)/norm(xk))
        %         toc(ta);
    end
    tk = tkp1;
    tkm1 = tk;
    zk = zkp1;
    xkm1 = xk;
    xk = xkp1;
    qk = qkp1;
    ck = ckp1;
    
end
end
function z = normpq(x,p,q,gLen)
n = length(x);
gNo = n/gLen;
tmp = reshape(x,gLen,gNo);
z = norm((sum(abs(tmp).^(p),1)).^(1/p),q)^q;
end