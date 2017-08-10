function pref = fitPref(dataFits)
% data should be orienations x cells
fit_ind = ~isnan(mean(dataFits,1)) & mean(dataFits,1) ~= 0;
[~, fit_pref_ind] = max(dataFits,[],1);
fit_pref_ind(~fit_ind) = NaN;
fit_pref_ind = fit_pref_ind-1;
fit_pref_ind(fit_pref_ind == 180) = 0;
pref = fit_pref_ind;
end