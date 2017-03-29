function dFoverF = getDFoverF(F_trials_cells,Find)
    
F = mean(F_trials_cells(Find,:,:),1);
dF = bsxfun(@minus,F_trials_cells,F);
dFoverF = bsxfun(@rdivide,dF,F);

end