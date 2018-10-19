function pctCorr = getPctCorr_trainData(glmResult,X,Y,dv)

yhat = glmval(glmResult.beta,X,'logit') > dv;
pctCorr = mean((Y-yhat) == 0);

end