function pctCorr = getPctCorr_hoData_lass(X,Y,dv)

nt = length(Y);
yhat = nan(nt,1);
for i = 1:nt
    X_holdout = X(i,:);
    othersInd = [1:(i-1),(i+1):nt];
    
    [B,othersGLM] = lassoglm(X(othersInd,:),Y(othersInd),'binomial');
    yhat(i) = glmval([othersGLM.Intercept(1);B(:,1)],X_holdout,'logit') > dv;
end
pctCorr = mean((Y-yhat) == 0);

end