function pctCorr = getPctCorr_hoData(X,Y,dv)

nt = length(Y);
yhat = nan(nt,1);
for i = 1:nt
    X_holdout = X(i,:);
    othersInd = [1:(i-1),(i+1):nt];
    
    [~,~,othersGLM] = glmfit(X(othersInd,:),Y(othersInd),'binomial');
    yhat(i) = glmval(othersGLM.beta,X_holdout,'logit') > dv;
end
pctCorr = mean((Y-yhat) == 0);

end