function [pctCorr,isCorrect] = getPctCorr_trainData(glmResult,X,Y,dv)

yhat = glmval(glmResult.beta,X,'logit') > dv;
isCorrect = (Y-yhat) == 0;
pctCorr = mean(isCorrect);


end