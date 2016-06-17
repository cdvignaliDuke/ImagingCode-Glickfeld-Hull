ttTarEndInd = find(chop(ttMs,3) == 1000);
ttMsCrop = ttMs(1:ttTarEndInd);

%%
set(0,'defaultfigurepaperorientation','portrait'); 
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
%sort all trials
[peakR_all_sorted peakR_all_sort_ind] = sort(resp_all_all);
%%
respHeatmap_allTrials = figure;
suptitle('avg all trials - sort across all trials')
colormap(brewermap(1001,'*RdBu'))
subplot(2,2,1)
allTrials_sort = imagesc(tc_all(20:ttTarEndInd,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('all trials')
axis square

subplot(2,2,2)
visTrials_allSort = imagesc(tc_val_all(20:ttTarEndInd,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('valid')
axis square
subplot(2,2,3)
audTrials_allSort = imagesc(tc_inv_all(20:ttTarEndInd,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('invalid')
axis square

% concatenate DG TCs
oriTuningResp_tc_crop = oriTuningResp_tc(:,6:end,:);
sz = size(oriTuningResp_tc_crop);
oriTuningResp_catTC = reshape(permute(oriTuningResp_tc_crop,[2 1 3]),sz(1)*sz(2),sz(3));
oriBorders = sz(2):sz(2):sz(2)*4;

subplot(2,2,4)
oriTC_visSort = imagesc(oriTuningResp_catTC(:,fliplr(peakR_all_sort_ind))');
hold on
vline(oriBorders,'k')
xticks = [4 (4 +sz(2)) (4+(sz(2)*2)) (4+(sz(2)*3))];
xticklabels = {'0deg';'45deg';'90deg';'135deg'};
ax = gca;
ax.XTick = xticks;
ax.XTickLabel = xticklabels;
xlabel('time (s)')
ylabel('cells')
caxis([-0.3 0.3])
colorbar
title(['DG stim tc,tick at 0s'])
axis square

figure(respHeatmap_allTrials)
print([fnout 'catch_align_TCheatmatp_allTrSortPlustuningData' datasetStr ], '-dpdf')

%% last plot vis-aud, sorted by all

respHeatmap_allTrials_VsubA = figure;
suptitle('avg all trials - sort across all trials')
colormap(brewermap(1001,'*RdBu'))
subplot(2,2,1)
allTrials_sort = imagesc(tc_all(20:ttTarEndInd,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('all trials')
axis square

subplot(2,2,2)
visTrials_allSort = imagesc(tc_val_all(20:ttTarEndInd,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('valid')
axis square

subplot(2,2,3)
audTrials_allSort = imagesc(tc_inv_all(20:ttTarEndInd,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('invalid')
axis square

VsubA_tc_all = bsxfun(@minus, tc_val_all(20:ttTarEndInd,fliplr(peakR_all_sort_ind))', tc_inv_all(20:ttTarEndInd,fliplr(peakR_all_sort_ind))');

subplot(2,2,4)
VsubA = imagesc(VsubA_tc_all);
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('val - inv')
axis square

figure(respHeatmap_allTrials_VsubA)
print([fnout 'catch_align_TCheatmatp_allTrSortPlusValsubInv' datasetStr ], '-dpdf')


