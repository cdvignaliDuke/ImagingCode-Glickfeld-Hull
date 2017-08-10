function tuning_shifted = alignPeakOriTuning(tuning_mat,oris)
% Align orientation tuning curves to peak
% tuning mat should be responses of size n oris x n cells
[~, center_ind] = min(abs(oris-90));
center_ind = center_ind-1;
[~, max_ind] = max(tuning_mat,[],1);
nc = size(tuning_mat,2);

tuning_shifted = nan(length(oris),nc);
for icell = 1:nc
   peak = max_ind(icell);
   tc = tuning_mat(:,icell);   
   nshift = center_ind-peak;
   tuning_shifted(:,icell) = circshift(tc,nshift);
end

end