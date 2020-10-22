function [pctCorr,isCorrect] = getPctCorr_trainData(glmResult,X,Y,dv,varargin)

if isempty(varargin)
    pc_ind = 1:size(X,2);
else
    pc_ind = varargin{1};
end

yhat = glmval(glmResult.beta([1,(pc_ind+1)]),X(:,pc_ind),'logit') > dv;
isCorrect = (Y-yhat) == 0;
pctCorr = mean(isCorrect);


end