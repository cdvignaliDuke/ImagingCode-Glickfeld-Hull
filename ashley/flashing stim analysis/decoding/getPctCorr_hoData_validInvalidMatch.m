function [pctCorr_val,pctCorr_inv] = getPctCorr_hoData_validInvalidMatch(...
    X_valid_all,Y_valid_all,groupIndNumber,dv,X_invalid,Y_invalid)

X_valid_subgroup = X_valid_all(groupIndNumber,:);
Y_subgroup = Y_valid_all(groupIndNumber);
nt_group = length(groupIndNumber);
nt_all = length(Y_valid_all);
nt_inv = length(Y_invalid);

yhat_val = nan(nt_group,1);
yhat_inv = nan(nt_group,nt_inv);
for i = 1:nt_group
    X_holdout = X_valid_subgroup(i,:);
    
    groupTrialInd = groupIndNumber(i);
    othersInd = [1:(groupTrialInd-1),(groupTrialInd+1):nt_all];

    [~,~,othersGLM] = glmfit(X_valid_all(othersInd,:),Y_valid_all(othersInd),'binomial');
    yhat_val(i) = glmval(othersGLM.beta,X_holdout,'logit') > dv;
    yhat_inv(i,:) = glmval(othersGLM.beta,X_invalid,'logit') > dv;
end
pctCorr_val = mean((Y_subgroup-yhat_val) == 0);
pctCorr_inv = mean((yhat_inv-Y_invalid') == 0);

end