
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
%%
respHeatmap_vis = figure;
% colormap('hot')
colormap(brewermap([],'*RdBu'))
suptitle('visual trials')
respHeatmap_aud = figure;
% colormap('hot')
colormap(brewermap([],'*RdBu'))
suptitle('auditory trials')

cell_ind_ori = cell(4,1);
L = 0;

i = 1;
tc_vis_all = [];
tc_aud_all = [];
tc_all = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(1).ind,cell_ind);
        for iOri = 1:4
            if i == 1
                [int c] = intersect(cell_ind,mouse(imouse).expt(iexp).cells(iOri+1).ind);
                if ~isempty(c)
                cell_ind_ori{iOri,1} = cat(1,cell_ind_ori{iOri},c);
                end
            else
                [int c] = intersect(cell_ind,mouse(imouse).expt(iexp).cells(iOri+1).ind);
                if ~isempty(c)
                cell_ind_ori{iOri,1} = cat(1,cell_ind_ori{iOri},c+L);
                end
            end
        end
        L = length(cell_ind)+L;
        peakR_vis = squeeze(mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_event_frames:end,:,:),3),1));
        [peakR_vis_sorted peakR_vis_sort_ind] = sort(peakR_vis);
        cell_ind_sort = flipud(intersect(peakR_vis_sort_ind,cell_ind,'stable'));
        tc_vis_all = cat(2, tc_vis_all, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(20:end,cell_ind,:),3));
        tc_aud_all = cat(2, tc_aud_all, mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(20:end,cell_ind,:),3));
        tc_all = cat(2,tc_all, mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(20:end,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(20:end,cell_ind,:)),3));
        figure(respHeatmap_vis);
        subplot(n,n2,i)
        imagesc(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(20:end,cell_ind_sort,:),3)')
%         axis square
        ax = gca;
        ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
        ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
        caxis([-0.01 0.1])
        figure(respHeatmap_aud);
        subplot(n,n2,i)
        imagesc(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(20:end,cell_ind_sort,:),3)')
%         axis square
        ax = gca;
        ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
        ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
        caxis([-0.01 0.1])
        i = i+1;
    end 
end
%%

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

peakR_vis = squeeze(mean(tc_vis_all(pre_event_frames-20:end,:),1));
[peakR_vis_sorted peakR_vis_sort_ind] = sort(peakR_vis);
% sort auditory trials
peakR_aud = squeeze(mean(tc_aud_all(pre_event_frames-20:end,:),1));
[peakR_aud_sorted peakR_aud_sort_ind] = sort(peakR_aud);
%sort all trials
peakR_all = squeeze(mean(tc_all(pre_event_frames-20:end,:),1));
[peakR_all_sorted peakR_all_sort_ind] = sort(peakR_all);
%%
% plot groups of cells
respHeatmap_all = figure;
suptitle('avg all trials - sort vis, subtr vis-aud')
% colormap('hot')
colormap(brewermap([],'*RdBu'))
subplot(5,3,1)
visTrials_visSort = imagesc(tc_vis_all(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
% caxis([-0.05 0.3])
% ax = gca;
% clim_vis = [min(min(tc_vis_all)) 0.1];
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('visual')
axis square
subplot(5,3,2)
audTrials_visSort = imagesc(tc_aud_all(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
% caxis([-0.05 0.3])
colorbar
title('auditory')
axis square
subplot(5,3,3)
VisAudSubtr_visSort = imagesc(tc_vis_all(:,fliplr(peakR_vis_sort_ind))'-tc_aud_all(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([-0.08 0.08])
colorbar
title('subtract')
axis square

plt = 4:3:15;
for iOri = 1:4

    subplot(5,3,plt(iOri))
    imagesc(tc_vis_all(:,flipud(intersect(peakR_vis_sort_ind,cell_ind_ori{iOri},'stable')))');
    ax = gca;
    ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
    ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
    xlabel('time (s)')
    ylabel('cells')
    caxis(clim_vis)
    colorbar
    title([mouse(1).expt(1).cells(iOri+1).name 'deg slctv'])
    axis square
    
    subplot(5,3,plt(iOri)+1)
    imagesc(tc_aud_all(:,flipud(intersect(peakR_vis_sort_ind,cell_ind_ori{iOri},'stable')))');
    ax = gca;
    ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
    ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
    xlabel('time (s)')
    ylabel('cells')
    caxis(clim_vis)
    colorbar
    title(['n=' num2str(length(cell_ind_ori{iOri}))])
    axis square
    
    subplot(5,3,plt(iOri)+2)
    imagesc(tc_vis_all(:,flipud(intersect(peakR_vis_sort_ind,cell_ind_ori{iOri},'stable')))'-tc_aud_all(:,flipud(intersect(peakR_vis_sort_ind,cell_ind_ori{iOri},'stable')))');
    ax = gca;
    ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
    ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
    xlabel('time (s)')
    ylabel('cells')
    caxis(clim_vis)
    colorbar
    title('subtract')
    axis square
    
end

%%
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);


F2 = figure;
suptitle('avg all trials - sort vis, sort aud')
% colormap('hot')
colormap(brewermap([],'*RdBu'))
subplot(2,3,1)
imagesc(tc_vis_all(:,fliplr(peakR_vis_sort_ind))');
title('visual-sort')
subplot(2,3,2)
imagesc(tc_aud_all(:,fliplr(peakR_vis_sort_ind))');
title('auditory')
subplot(2,3,3)
imagesc(tc_vis_all(:,fliplr(peakR_vis_sort_ind))'-tc_aud_all(:,fliplr(peakR_vis_sort_ind))');
title('subtract')
subplot(2,3,4)
imagesc(tc_vis_all(:,fliplr(peakR_aud_sort_ind))');
title('visual')
subplot(2,3,5)
imagesc(tc_aud_all(:,fliplr(peakR_aud_sort_ind))');
title('auditory-sort')
subplot(2,3,6)
imagesc(tc_vis_all(:,fliplr(peakR_aud_sort_ind))'-tc_aud_all(:,fliplr(peakR_aud_sort_ind))');
title('subtract')

for iplot = 1:6
   figure(F2)
   subplot(2,3,iplot)
   hold on
   ax = gca;
    ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
    ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
    xlabel('time (s)')
    ylabel('cells')
    caxis([-0.1 0.1])
% clim_vis = [min(min(tc_vis_all)) 0.1];
% caxis(clim_vis)
    colorbar
    axis square
end
%% plot with aud-vis subtraction
F1 = figure;
suptitle('avg all trials - sort vis, subtr aud-vis')
% colormap('hot')
colormap(brewermap([],'*RdBu'))
subplot(1,3,2)
visTrials_visSort = imagesc(tc_vis_all(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
% caxis([-0.05 0.3])
% ax = gca;
clim_vis = [-0.1 0.1];
caxis(clim_vis)
colorbar
title('visual')
axis square
subplot(1,3,1)
audTrials_visSort = imagesc(tc_aud_all(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
% caxis([-0.05 0.3])
colorbar
title('auditory')
axis square
subplot(1,3,3)
AudVisSubtr_visSort = imagesc(tc_aud_all(:,fliplr(peakR_vis_sort_ind))'-tc_vis_all(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([-0.08 0.08])
colorbar
title('subtract')
axis square

%%
figure(respHeatmap_vis);
print([fnout 'press_align_TCheatmatp_visexpt' datasetStr '.pdf'], '-dpdf')
figure(respHeatmap_aud);
print([fnout 'press_align_TCheatmatp_audexpt' datasetStr '.pdf'], '-dpdf')
figure(respHeatmap_all);
print([fnout 'press_align_TCheatmatp_all_aubV-A' datasetStr '.pdf'], '-dpdf')
figure(F1);
print([fnout 'press_align_TCheatmatp_all_subA-V' datasetStr '.pdf'], '-dpdf')
figure(F2);
print([fnout 'press_align_TCheatmatp_sortVisandAud' datasetStr '.pdf'], '-dpdf')
%% plot ori tuning data alongside
oriTuningResp_avg = [];
oriTuningResp_tc = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(1).ind,cell_ind);
        mName = ['AW' mouse(imouse).expt(iexp).mouse_name(2:end)];
        dirtuning = expt(intersect( find(strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)) ,find(strcmp({expt.date},mouse(imouse).expt(iexp).date)) ) ).dirtuning;
        load(fullfile(rc.ashleyAnalysis,mName,'two-photon imaging',mouse(imouse).expt(iexp).date,dirtuning,'cellsSelect.mat'))
        oriTuningResp_avg = cat(2,oriTuningResp_avg,dFoverF_OriResp_avg_rect(:,cell_ind));
        oriTuningResp_tc = cat(3,oriTuningResp_tc,dFoverF_OriResp_TC(:,:,cell_ind));
    end
end
% plot groups of cells
respHeatmapPlusOri = figure;
suptitle('avg all trials - sort vis')
% colormap('hot')
colormap(brewermap([],'*RdBu'))
subplot(2,4,1)
visTrials_visSort = imagesc(tc_vis_all(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
% caxis([-0.05 0.3])
% ax = gca;
% clim_vis = [min(min(tc_vis_all)) 0.1];
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('visual')
axis square
subplot(2,4,2)
audTrials_visSort = imagesc(tc_aud_all(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
% caxis([-0.05 0.3])
colorbar
title('auditory')
axis square

oris = [0 45 90 135];

subplot(2,4,3)
oriResp_visSort = imagesc(oriTuningResp_avg(:,fliplr(peakR_vis_sort_ind))');
ax = gca;
ax.XTick = 1:4;
ax.XTickLabel = oris;
xlabel('ori')
ylabel('cells')
% caxis(clim_vis)
caxis([-max(max(oriTuningResp_avg)) max(max(oriTuningResp_avg))])
colorbar
title('DG resp')
axis square

plt = 5:8;
for iOri = 1:4
subplot(2,4,plt(iOri))
oriTC_visSort = imagesc(squeeze(oriTuningResp_tc(iOri,:,fliplr(peakR_vis_sort_ind)))');
ax = gca;
ax.XTick = [ 7 10 13];
ax.XTickLabel = [-1 0 1];
xlabel('time (s)')
ylabel('cells')
% caxis(clim_vis)
caxis([-max(max(max(oriTuningResp_tc))) max(max(max(oriTuningResp_tc)))])
colorbar
title([num2str(oris(iOri)) ' deg DG stim'])
axis square
end

figure(respHeatmapPlusOri)
print([fnout 'press_align_TCheatmatp_visSortPlustuningData' datasetStr '.pdf'], '-dpdf')

%% sort by avg across all trials

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);

respHeatmap_allTrials = figure;
suptitle('avg all trials - sort across all trials')
colormap(brewermap(1001,'*RdBu'))
subplot(2,2,1)
allTrials_sort = imagesc(tc_all(:,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('all trials')
axis square

subplot(2,2,2)
visTrials_allSort = imagesc(tc_vis_all(:,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('visual')
axis square
subplot(2,2,3)
audTrials_allSort = imagesc(tc_aud_all(:,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
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
print([fnout 'press_align_TCheatmatp_allTrSortPlustuningData' datasetStr ], '-dpdf')

%% last plot vis-aud, sorted by all

respHeatmap_allTrials_VsubA = figure;
suptitle('avg all trials - sort across all trials')
colormap(brewermap(1001,'*RdBu'))
subplot(2,2,1)
allTrials_sort = imagesc(tc_all(:,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('all trials')
axis square

subplot(2,2,2)
visTrials_allSort = imagesc(tc_vis_all(:,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.1 0.1];
caxis(clim_vis)
colorbar
title('visual')
axis square

subplot(2,2,3)
audTrials_allSort = imagesc(tc_aud_all(:,fliplr(peakR_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('auditory')
axis square

VsubA_tc_all = bsxfun(@minus, tc_vis_all(:,fliplr(peakR_all_sort_ind))', tc_aud_all(:,fliplr(peakR_all_sort_ind))');

subplot(2,2,4)
VsubA = imagesc(VsubA_tc_all);
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('vis - aud')
axis square

figure(respHeatmap_allTrials_VsubA)
print([fnout 'press_align_TCheatmatp_allTrSortPlustVsubA' datasetStr ], '-dpdf')