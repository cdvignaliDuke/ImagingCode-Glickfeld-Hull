function [pctCorr,isCorrect] = getPctCorr_weightsOnly(beta,X,Y,dv)

yhat = glmval(beta,X,'logit') > dv;
isCorrect = (Y-yhat) == 0;
pctCorr = mean(isCorrect);


end