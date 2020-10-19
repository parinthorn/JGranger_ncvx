function score = performance_score(stat)
TP=stat(:,:,1);
TN=stat(:,:,2);
FP=stat(:,:,3);
FN=stat(:,:,4);

score.FPR = FP./(FP+TN);
score.TPR = TP./(TP+FN);
score.ACC = (TP+TN)./(TP+FN+TN+FP);
score.F1 = 2*TP./(2*TP+FP+FN);
end