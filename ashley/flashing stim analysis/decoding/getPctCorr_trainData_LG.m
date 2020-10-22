function [pctCorr MI] = getPctCorr_trainData_LG(weights,X,Y,dv)

prob_out = glmval(weights,X,'logit');
yhat = prob_out>dv;
pctCorr = mean((Y-yhat) == 0);
MI = sum(prob_out.*log(prob_out) + ((1-prob_out)./(log(1-prob_out))));

end