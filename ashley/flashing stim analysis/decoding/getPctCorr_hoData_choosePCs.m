function pctCorr = getPctCorr_hoData_choosePCs(X,Y,dv,pc_ind)

nt = length(Y);
yhat = nan(nt,1);
for i = 1:nt
    X_holdout = X(i,:);
    othersInd = [1:(i-1),(i+1):nt];
    
    [~,~,othersGLM] = glmfit(X(othersInd,:),Y(othersInd),'binomial');
    yhat(i) = glmval(othersGLM.beta([1,(pc_ind+1)]),X_holdout(:,pc_ind),'logit') > dv;
end
pctCorr = mean((Y-yhat) == 0);

end