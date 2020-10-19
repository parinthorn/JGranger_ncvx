%% This experiment estimate VAR with formulation S by ADMM
clear
clc
SFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS\FormulationS_Model20Sep19';
DataFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\THESIS_DATA_Sparse';
ModelFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\MODEL';
Mtype = {'C','D','S'}; % DATA TYPE
dp = [50,300];
[P,~] = offdiagJSS(15,2,4);
D = diffmat(15,2,4);
DP = sparse(D*P);
for b = 5:10 %C12
    s = load(strcat(ModelFolderName,'\MODEL20Sep19.mat'));
    fprintf('Databank number: %d \n',b)
    for tp = 1:3 % TEST
        Mname = Mtype{tp};
        for f=1:2 % TEST
            disp(strcat('Model name:',Mname,int2str(f)))
            r = load(strcat(DataFolderName,'\DATA_BANK_',int2str(b),'_',Mname,int2str(f),'.mat'));
%             trials = size(r.data.y,4);
            trials = 10;
            for dpx = 1:2
                E(trials).M =struct();
                Num = dp(dpx);
                fprintf('Data point : %d \n',Num)
                for t=1:trials
                    fprintf('Trial number : %d\n',t)
                    M = est_Formulation_S(r.data.y(:,1:Num,:,t),s.E(b).M.(strcat(Mname,int2str(f))).A,P,DP);
                    E(t).M = M;
                    E(t).M.DataInfo = r.data.info;
                    E(t).M.DataSize = Num;
                end
                save(strcat(SFolderName,'\ESTIMATED_BANK_',int2str(b),'_',Mname,int2str(f),'_',int2str(Num),'_adaptive_6trial.mat'),'E')
                clear E
            end
        end
    end
end