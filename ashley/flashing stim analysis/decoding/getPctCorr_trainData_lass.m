function pctCorr = getPctCorr_trainData_lass(beta,X,Y,dv)

yhat = glmval(beta,X,'logit') > dv;
pctCorr = mean((Y-yhat) == 0);

end