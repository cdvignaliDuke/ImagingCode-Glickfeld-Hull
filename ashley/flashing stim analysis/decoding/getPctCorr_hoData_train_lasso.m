function [pctCorr] = getPctCorr_hoData_train_lasso(X,Y,ind_sub,ind_train,dv)
%X- responses for all trials all cells
%Y- category for all trials
%ind_sub- set of trials to test
%ind_train- set of trials to train;
nt = length(ind_sub);
yhat = nan(nt,1);
C = eye(size(X,2));
[B FitInfo] = lassoglm(X(ind_train,:)*C,Y(ind_train),'binomial');
B0 = FitInfo.Intercept(75);
targetWeight = [B0; B(:,75)];
for i = 1:nt
    ind_test = ind_sub(i);
    X_holdout = X(ind_test,:);
    if find(ind_train==ind_test)
        othersInd = setdiff(ind_train,ind_test);

        [B FitInfo] = lassoglm(X(othersInd,:),Y(othersInd),'binomial');
        B0 = FitInfo.Intercept(75);
        othersWeight = [B0; B(:,75)];
        yhat(i) = glmval(othersWeight,X_holdout,'logit') > dv;
    else
        yhat(i) = glmval(targetWeight,X_holdout,'logit') > dv; 
    end
end
pctCorr = mean((Y(ind_sub)-yhat) == 0);

end