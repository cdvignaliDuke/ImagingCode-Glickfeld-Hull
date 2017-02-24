F = mean(data_TC,1);
dF = bsxfun(@minus,data_TC,F);
dFoverF = bsxfun(@rdivide,dF,F);

dFoverF_diff = diff(dFoverF);
dF_mean = mean(dFoverF,1);
dF_std = std(dFoverF,[],1);

thr2std = 2*dF_std;

dFdiffThr = bsxfun(@gt,dFoverF_diff,dF_std);
dFdiffThr = cat(1,dFdiffThr,zeros(1,size(dFoverF,2)));

bin_siz = 10;
dFdiffThr_binned = reshape(dFdiffThr,bin_siz,size(dFdiffThr,1)/bin_siz,size(dFdiffThr,2));
dFoverF_binsum = reshape(squeeze(mean(dFdiffThr_binned,1)),size(dFdiffThr,1)/bin_siz,size(dFdiffThr,2));

diffCorr = corrcoef(dFdiffThr);
diffCorr = tril(diffCorr,-1);

diffBinnedCorr = corrcoef(dFoverF_binsum);
diffBinnedCorr = tril(diffBinnedCorr,-1);

bar_bins = -1:0.05:1;
diffCorr_hist = hist(diffCorr(diffCorr ~=0),bar_bins);
diffBinnedCorr_hist = hist(diffBinnedCorr(diffBinnedCorr ~=0),bar_bins);
diffMean = mean(diffCorr(diffCorr ~= 0));

figure;
suptitle('thresholded time-courses')
colormap(brewermap([],'*RdBu'))
subplot(2,2,1)
bar(bar_bins,diffCorr_hist)
hold on
vline(diffMean,'k--')
xlim([-1 1])
xlabel('r')
ylabel('n_pairs')
axis square
title({['iti; mean corr' num2str(chop(diffMean,2))]})
subplot(2,2,2)
imagesc(diffCorr)
axis square
colorbar
caxis([-1 1])
title({'iti; corr of events'})
print([fnpath '\' exptName 'thresholdedItiCorrs.pdf'],'-dpdf','-fillpage')
