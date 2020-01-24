function [yFit,yFit_err,rsq,slope, pVal,corrVal] = fitLine(x,y)
    mdl = fitlm(x,y);
    rsq = mdl.Rsquared.Ordinary;
    [yFit,yFit_ci] = predict(mdl,x);
    yFit_err = yFit-yFit_ci(:,1);
    slope = mdl.Coefficients.Estimate(2);
    [corrmat,pmat] = corrcoef(x,y);
    pVal = pmat(1,2);
    corrVal = corrmat(1,2);
end