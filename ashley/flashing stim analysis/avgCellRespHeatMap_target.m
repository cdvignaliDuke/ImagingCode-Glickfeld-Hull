ttTarEndInd = find(chop(ttMs,3) == 1000);
ttMsCrop = ttMs(1:ttTarEndInd);

respHeatmap_vis = figure;
colormap('hot')
suptitle('visual trials')
respHeatmap_aud = figure;
colormap('hot')
suptitle('auditory trials')

i = 1;
tc_vis_all = [];
tc_aud_all = [];
tc_all = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(13).ind;
        peakR_vis = squeeze(mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(pre_event_frames:pre_event_frames+mouse(imouse).expt(iexp).info.cyc_time,:,:),3),1));
        [peakR_vis_sorted peakR_vis_sort_ind] = sort(peakR_vis);
        cell_ind_sort = flipud(intersect(peakR_vis_sort_ind,cell_ind,'stable'));
        tc_vis_all = cat(2, tc_vis_all, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(20:end,cell_ind,:),3));
        tc_aud_all = cat(2, tc_aud_all, mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(20:end,cell_ind,:),3));
        tc_all = cat(2,tc_all, mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(20:end,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(20:end,cell_ind,:)),3));
%         figure(respHeatmap_vis);
%         subplot(n,n2,i)
%         imagesc(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(20:end,cell_ind_sort,:),3)')
% %         axis square
%         ax = gca;
%         ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)/2 find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
%         ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
%         caxis([-0.01 0.1])
%         figure(respHeatmap_aud);
%         subplot(n,n2,i)
%         imagesc(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(20:end,cell_ind_sort,:),3)')
% %         axis square
%         ax = gca;
%         ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)/2 find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
%         ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
%         caxis([-0.01 0.1])
        i = i+1;
    end 
end

%% drifting grating responses
oriTuningResp_avg = [];
oriTuningResp_tc = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(13).ind;
        mName = ['AW' mouse(imouse).expt(iexp).mouse_name(2:end)];
        dirtuning = expt(intersect( find(strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)) ,find(strcmp({expt.date},mouse(imouse).expt(iexp).date)) ) ).dirtuning;
        load(fullfile(rc.ashleyAnalysis,mName,'two-photon imaging',mouse(imouse).expt(iexp).date,dirtuning,'cellsSelect.mat'))
        oriTuningResp_avg = cat(2,oriTuningResp_avg,dFoverF_OriResp_avg_rect(:,cell_ind));
        oriTuningResp_tc = cat(3,oriTuningResp_tc,dFoverF_OriResp_TC(:,:,cell_ind));
    end
end
%%
set(0,'defaultfigurepaperorientation','portrait'); 
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
%sort all trials
% % peakR_all = squeeze(mean(tc_all(pre_event_frames-20:end,:),1));
peakR_all = squeeze(mean(tc_all(pre_event_frames-20:pre_event_frames-20+10,:),1));
[peakR_all_sorted peakR_all_sort_ind] = sort(peakR_all);
%%
respHeatmap_allTrials = figure;
suptitle('avg all trials - sort across all trials')
colormap(brewermap(1001,'*RdBu'))
subplot(2,2,1)
allTrials_sort = imagesc(tc_all(1:ttTarEndInd-20,fliplr(peakR_all_sort_ind))');
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
visTrials_allSort = imagesc(tc_vis_all(1:ttTarEndInd-20,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('visual')
axis square
subplot(2,2,3)
audTrials_allSort = imagesc(tc_aud_all(1:ttTarEndInd-20,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('auditory')
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
print([fnout 'target_align_TCheatmatp_allTrSortPlustuningData' datasetStr ], '-dpdf')

%%
%% last plot vis-aud, sorted by all

respHeatmap_allTrials_VsubA = figure;
suptitle('avg all trials - sort across all trials')
colormap(brewermap(1001,'*RdBu'))
subplot(2,2,1)
allTrials_sort = imagesc(tc_all(1:ttTarEndInd-20,fliplr(peakR_all_sort_ind))');
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
visTrials_allSort = imagesc(tc_vis_all(1:ttTarEndInd-20,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('visual')
axis square

subplot(2,2,3)
audTrials_allSort = imagesc(tc_aud_all(1:ttTarEndInd-20,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('auditory')
axis square

VsubA_tc_all = bsxfun(@minus, tc_vis_all(1:ttTarEndInd-20,fliplr(peakR_all_sort_ind))', tc_aud_all(1:ttTarEndInd-20,fliplr(peakR_all_sort_ind))');

subplot(2,2,4)
VsubA = imagesc(VsubA_tc_all);
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('vis - aud')
axis square

figure(respHeatmap_allTrials_VsubA)
print([fnout 'target_align_TCheatmatp_allTrSortPlustVsubA' datasetStr ], '-dpdf')
%%

peakR_vis = squeeze(mean(tc_vis_all(pre_event_frames-20:pre_event_frames-20+mouse(imouse).expt(iexp).info.cyc_time,:),1));
[peakR_vis_sorted peakR_vis_sort_ind] = sort(peakR_vis);
respHeatmap_all = figure;
colormap('hot')
subplot(1,3,1)
imagesc(tc_vis_all(:,fliplr(peakR_vis_sort_ind))')
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([-0.01 0.1])
colorbar
title('visual')
subplot(1,3,2)
imagesc(tc_aud_all(:,fliplr(peakR_vis_sort_ind))')
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([-0.01 0.1])
colorbar
title('auditory')
subplot(1,3,3)
imagesc(tc_vis_all(:,fliplr(peakR_vis_sort_ind))'-tc_aud_all(:,fliplr(peakR_vis_sort_ind))')
ax = gca;
ax.XTick = [find(ceil(ttMsCrop(20:end)) == 0) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000/2) find(ttMsCrop(20:end) == floor(max(ttMsCrop(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMsCrop(20:end)/1000))/2 floor(max(ttMsCrop(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([-0.01 0.1])
colorbar
title('subtract')

figure(respHeatmap_vis);
print([fnout 'target_align_TCheatmatp_visexpt' datasetStr '.pdf'], '-dpdf')
figure(respHeatmap_aud);
print([fnout 'target_align_TCheatmatp_audexpt' datasetStr '.pdf'], '-dpdf')
figure(respHeatmap_all);
print([fnout 'target_align_TCheatmatp_all' datasetStr '.pdf'], '-dpdf')

%%

