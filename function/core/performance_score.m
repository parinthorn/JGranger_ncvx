function score = performance_score(stat)
% This function compute performance index from confution matrix
% Input
% stat: confusion matrix
% Output
% score: performance score [F1, ACC, TPR, MCC]
%
% Originally written by Parinthorn Manomaisaowapak
% Please email to parinthorn@gmail.com before reuse, reproduce

if size(stat,3) ==4
    TP=stat(:,:,1);
    TN=stat(:,:,2);
    FP=stat(:,:,3);
    FN=stat(:,:,4);
else
    TP=stat(1);
    TN=stat(2);
    FP=stat(3);
    FN=stat(4);
end
score.FPR = FP./(FP+TN);
if score.FPR>1
    error('something wrong')
end
score.TPR = TP./(TP+FN);
score.ACC = (TP+TN)./(TP+FN+TN+FP);
score.F1 = 2*TP./(2*TP+FP+FN);
tmp = sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));
tmp(tmp==0) = inf;
score.MCC = (TP.*TN-FP.*FN)./tmp;
end