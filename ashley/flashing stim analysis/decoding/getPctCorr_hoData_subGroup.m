function pctCorr = getPctCorr_hoData_subGroup(X_all,Y_all,groupIndNumber,dv)

X_subgroup = X_all(groupIndNumber,:);
Y_subgroup = Y_all(groupIndNumber);
nt_group = length(groupIndNumber);
nt_all = length(Y_all);

yhat = nan(nt_group,1);
for i = 1:nt_group
    X_holdout = X_subgroup(i,:);
    
    groupTrialInd = groupIndNumber(i);
    othersInd = [1:(groupTrialInd-1),(groupTrialInd+1):nt_all];
    
    [~,~,othersGLM] = glmfit(X_all(othersInd,:),Y_all(othersInd),'binomial');
    yhat(i) = glmval(othersGLM.beta,X_holdout,'logit') > dv;
end
pctCorr = mean((Y_subgroup-yhat) == 0);

end