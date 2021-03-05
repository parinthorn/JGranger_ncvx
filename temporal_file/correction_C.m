function M_out = correction_C(y,M,varargin)
%% This program estimates Ground truth with their original model
% Input are
%             y : 3D array [n,Num,K] which is dimension of timeseries,
%       timepoints, # models respectively
%             P : Off-diagonal projection matrix
% (optional)  p : lag-order
%      GridSize : dimension of a regularization grid
% Output is
%       M : structure containing
%           M.x_est
M_out = M;
[n,T,K] = size(y);
len_varargin = length(varargin);
if isempty(varargin)
  p=1;
  GridSize = 30;
  toggle = 'LLH_full';
  data_concat=0;
elseif len_varargin==1
  p = varargin{1};
  GridSize = 30;
  toggle = 'LLH_full';
  data_concat = 0;
elseif len_varargin ==2
  p = varargin{1};
  GridSize = varargin{2};
  toggle = 'LLH_full';
  data_concat = 0;
elseif len_varargin ==3
  p = varargin{1};
  GridSize = varargin{2};
  toggle = varargin{3};
  data_concat=0;
elseif len_varargin ==4
  p = varargin{1};
  GridSize = varargin{2};
  toggle = varargin{3};
  data_concat = varargin{4};
else
    error('must be atmost 6 input')
end
if (data_concat)
    Ktmp = K/2;K=2; % divides to 2 groups and set K=2
    eff_T = Ktmp*(T-p);
    y1 = y(:,:,1:Ktmp);
    y2 = y(:,:,(Ktmp+1):end);
    H = zeros(n*p,eff_T,2);
    Y = zeros(n,eff_T,2);
    disp('Concatenating H, Y matrix')
    for kk=1:Ktmp
        [H(:,(T-p)*(kk-1)+1:(T-p)*(kk),1),Y(:,(T-p)*(kk-1)+1:(T-p)*(kk),1)] = H_gen(y1(:,:,kk),p);
        [H(:,(T-p)*(kk-1)+1:(T-p)*(kk),2),Y(:,(T-p)*(kk-1)+1:(T-p)*(kk),2)] = H_gen(y2(:,:,kk),p);
    end
else
    eff_T = T-p;
    H = zeros(n*p,eff_T,K);
    Y = zeros(n,eff_T,K);
    disp('Generating H, Y matrix')
    for kk=1:K
        [H(:,:,kk),Y(:,:,kk)] = H_gen(y(:,:,kk),p);
    end
end
idx = efficient_vect([n,p,K]);
for ii=1:GridSize % test 20
    
    x_cls = M.model(ii).A(idx);
    A_cls =devect(full(x_cls),n,p,K);
    df = M.model(ii).stat.model_selection_score.df;
    score(1,ii) = model_selection_S(Y,H,A_cls,df,toggle);
    
end
GIC_LIST = fieldnames(score);
for nn=1:length(GIC_LIST)
    [~,M_out.index.(GIC_LIST{nn})] = min([score.(GIC_LIST{nn})]);
end
for ii=1:GridSize
    M_out.model(ii).stat.model_selection_score =  score(ii);
end
M_out.LLH_type = toggle;
end
