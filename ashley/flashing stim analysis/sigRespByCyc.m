endTrMI_VsubA = cellfun(@(x,y) bsxfun(@minus,x,y),endTrResp_vis,endTrResp_aud,'unif',false);

%% sort by modulation index, then plot resp to all trials, vis trials, aud trials, and vis-aud subtractions
% sorting
HMcyc = 8;
[MI_all_sorted MI_all_sort_ind] = sort(endTrMI_VsubA{HMcyc});
% plotting
respHeatmap_allTrials_VsubA = figure;
suptitle('avg all trials - sort by mod ix vis-aud endTrResp')
colormap(brewermap(1001,'*RdBu'))
subplot(2,2,1)
allTrials_sort = imagesc(tc_all(:,fliplr(MI_all_sort_ind))');
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
visTrials_allSort = imagesc(tc_vis_all(:,fliplr(MI_all_sort_ind))');
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
audTrials_allSort = imagesc(tc_aud_all(:,fliplr(MI_all_sort_ind))');
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
caxis([clim_vis])
colorbar
title('auditory')
axis square

VsubA_tc_all = bsxfun(@minus, tc_vis_all(:,fliplr(MI_all_sort_ind))', tc_aud_all(:,fliplr(MI_all_sort_ind))');

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
print([fnout 'press_align_TCheatmatp_MIsort_allTrPlusVsubA' datasetStr ], '-dpdf')

%% plot 'most modulated cells'
nc = 25;
sp1 = 5;
sp2 = 5;
cell2plot = MI_all_sort_ind(end-nc+1:end);%fliplr(MI_all_sort_ind(end-nc+1:end));
figure;
for i = 1:nc
    subplot(sp1,sp2,i)
    plot(ttMs(20:end),tc_vis_all(:,cell2plot(i)),'g')
    hold on
    plot(ttMs(20:end),tc_aud_all(:,cell2plot(i)),'k')
    hold on
    xlim([ttMs(20) ttMs(end)])
end

print([fnout 'press_align_MIallsub_exCellsMostMod' datasetStr ], '-dpdf')

%% find cells signif resp to each cycle compared to baseline

% cell_ind_ori = cell(4,1);
% L = 0;

ncyc = ceil(minTrialLengthFrames/cycTime);

tc_vis_allcells = [];
tc_aud_allcells = [];
tc_allcells = [];
mi_VsubA_late = [];
mi_VsubA_early = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        tc_vis_allcells = cat(2, tc_vis_allcells, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,:,:),3));
        tc_aud_allcells = cat(2, tc_aud_allcells, mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,:,:),3));
        tc_allcells = cat(2,tc_allcells, mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,:,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,:,:)),3));
        mi_VsubA_late = cat(2,mi_VsubA_late,bsxfun(@minus,squeeze(mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(end-cycTime*2:end,:,:),1),3)),squeeze(mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(end-cycTime*2:end,:,:),1),3))));
        mi_VsubA_early = cat(2,mi_VsubA_early,bsxfun(@minus,squeeze(mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(cycTime:cycTime*3,:,:),1),3)),squeeze(mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(cycTime:cycTime*3,:,:),1),3))));

    end 
end

pCyc_all = cell(1,ncyc);
hCyc_all = cell(1,ncyc);
p = [];
h = [];
% need to change alpha from 0.05???

for icyc = 1:ncyc
    for imouse = 1:size(mouse,2)
        for iexp = 1:size(mouse(imouse).expt,2)
        pw = mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_win,:,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(pre_win,:,:)),1);
        T = trans_win + (cycTime*(icyc-1));
        tw =  mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(T,:,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(T,:,:)),1);
        [htemp ptemp] = ttest(pw,tw,'dim', 3, 'tail', 'left');
        p = [p ptemp];
        h = [h htemp];
        end
    end
    pCyc_all{icyc} = p;
    hCyc_all{icyc} = h;
    p = [];
    h = [];
end

respCellCycInd = cellfun(@find,hCyc_all,'unif',false);
nRespCellsCyc = cell2mat(cellfun(@sum,hCyc_all,'unif',false));

% heatmap of resp cells for each cycle
c = ceil(sqrt(ncyc+1));
if (c^2)-c > ncyc+1
   c2 = c-1;
else
    c2 = c;
end

sigRespCellsEaCycFig_all = figure;
suptitle('resp to each cyc - all trials')
colormap(brewermap([],'*RdBu'))
sigRespCellsEaCycFig_vis = figure;
suptitle('resp to each cyc - visual trials')
colormap(brewermap([],'*RdBu'))
sigRespCellsEaCycFig_aud = figure;
suptitle('resp to each cyc - auditory trials')
colormap(brewermap([],'*RdBu'))
sigRespCellsEaCycFig_VsubA = figure;
suptitle('resp to each cyc - vis-aud')
colormap(brewermap([],'*RdBu'))

sigRespCellsEaCycFig_allmean = figure;
suptitle('resp to each cyc - all trials')
sigRespCellsEaCycFig_AVmean = figure;
suptitle('resp to each cyc - vis vs aud trials')


figure(sigRespCellsEaCycFig_all)
print([fnout 'press_align_cellsResp2EaCyc_alltrials'], '-dpdf')
figure(sigRespCellsEaCycFig_vis)
print([fnout 'press_align_cellsResp2EaCyc_vistrials'], '-dpdf')
figure(sigRespCellsEaCycFig_aud)
print([fnout 'press_align_cellsResp2EaCyc_audtrials'], '-dpdf')
figure(sigRespCellsEaCycFig_VsubA)
print([fnout 'press_align_cellsResp2EaCyc_VsubAtrials'], '-dpdf')
figure(sigRespCellsEaCycFig_allmean)
print([fnout 'press_align_cellsResp2EaCyc_allTrialsTC'], '-dpdf')
figure(sigRespCellsEaCycFig_AVmean)
print([fnout 'press_align_cellsResp2EaCyc_VvsATrialsTC'], '-dpdf')


for icyc = 1:ncyc
    if icyc == 1;
        C = respCellCycInd{icyc};
    else
        C = setdiff(respCellCycInd{icyc},unique(cell2mat(respCellCycInd(1:icyc-1))));
    end
    [sort_order sort_ind] = sort(squeeze(mean(tc_allcells(pre_event_frames:end,C),1)));
    C = C(sort_ind);
    
    %plot heatmap
    figure(sigRespCellsEaCycFig_all)
    subplot(c,c2,icyc)
    imagesc(tc_allcells(20:end,fliplr(C))');
    ax = gca;
    ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
    ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
    xlabel('time (s)')
    ylabel('cells')
    clim_vis = [-.1 0.1];
    caxis(clim_vis)
    colorbar
    title(['resp to cyc ' num2str(icyc)])
    axis square
    
    figure(sigRespCellsEaCycFig_vis)
    subplot(c,c2,icyc)
    imagesc(tc_vis_allcells(20:end,fliplr(C))');
    ax = gca;
    ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
    ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
    xlabel('time (s)')
    ylabel('cells')
    clim_vis = [-.1 0.1];
    caxis(clim_vis)
    colorbar
    title(['resp to cyc ' num2str(icyc)])
    axis square
    
    figure(sigRespCellsEaCycFig_aud)
    subplot(c,c2,icyc)
    imagesc(tc_aud_allcells(20:end,fliplr(C))');
    ax = gca;
    ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
    ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
    xlabel('time (s)')
    ylabel('cells')
    clim_vis = [-.1 0.1];
    caxis(clim_vis)
    colorbar
    title(['resp to cyc ' num2str(icyc)])
    axis square
    
    figure(sigRespCellsEaCycFig_VsubA)
    subplot(c,c2,icyc)
    imagesc(bsxfun(@minus,tc_vis_allcells(20:end,fliplr(C))',tc_aud_allcells(20:end,fliplr(C))'));
    ax = gca;
    ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
    ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
    xlabel('time (s)')
    ylabel('cells')
    clim_vis = [-.1 0.1];
    caxis(clim_vis)
    colorbar
    title(['resp to cyc ' num2str(icyc)])
    axis square
    
    %plot tc
    figure(sigRespCellsEaCycFig_allmean)
    subplot(c,c2,icyc)
    ste_all = std(tc_allcells(20:end,fliplr(C)),[],2)/sqrt(length(C));
    shadedErrorBar(ttMs(20:end),mean(tc_allcells(20:end,fliplr(C))',1),ste_all,'k');
    hold on
    vline(baseStimFrames/(cycTime/cycTimeMs),'k:')    
    ax = gca;
    ax.XTick = 0:500:2000;
    ax.XTickLabel = 0:500:2000;
    xlim([min(ttMs(20:end)) max(ttMs(20:end))]);
    ylim([-0.01 0.03])
    xlabel('time (ms)')
    ylabel('dF/F')
    title(['resp to cyc ' num2str(icyc) '; ' num2str(length(C)) ' cells'])
    axis square
    
    figure(sigRespCellsEaCycFig_AVmean)
    subplot(c,c2,icyc)
    ste_v = std(tc_vis_allcells(20:end,fliplr(C)),[],2)/sqrt(length(C));
    ste_a = std(tc_aud_allcells(20:end,fliplr(C)),[],2)/sqrt(length(C));
    shadedErrorBar(ttMs(20:end),mean(tc_vis_allcells(20:end,fliplr(C))',1),ste_v,'g');
    hold on
    shadedErrorBar(ttMs(20:end),mean(tc_aud_allcells(20:end,fliplr(C))',1),ste_a,'k');
    hold on
    vline(baseStimFrames/(cycTime/cycTimeMs),'k:') 
    ax = gca;
    ax.XTick = 0:500:2000;
    ax.XTickLabel = 0:500:2000;
    xlim([min(ttMs(20:end)) max(ttMs(20:end))]);
    ylim([-0.01 0.03])
    xlabel('time (ms)')
    ylabel('dF/F')
    title(['resp to cyc ' num2str(icyc) '; ' num2str(length(C)) ' cells'])
    axis square
end

figure(sigRespCellsEaCycFig_all)
print([fnout 'press_align_cellsResp2EaCyc_alltrials'], '-dpdf')
figure(sigRespCellsEaCycFig_vis)
print([fnout 'press_align_cellsResp2EaCyc_vistrials'], '-dpdf')
figure(sigRespCellsEaCycFig_aud)
print([fnout 'press_align_cellsResp2EaCyc_audtrials'], '-dpdf')
figure(sigRespCellsEaCycFig_VsubA)
print([fnout 'press_align_cellsResp2EaCyc_VsubAtrials'], '-dpdf')
figure(sigRespCellsEaCycFig_allmean)
print([fnout 'press_align_cellsResp2EaCyc_allTrialsTC'], '-dpdf')
figure(sigRespCellsEaCycFig_AVmean)
print([fnout 'press_align_cellsResp2EaCyc_VvsATrialsTC'], '-dpdf')
%% plot A vs V for rand cyc-1 responsive cells
c = 1;
if c == 1;
    C = respCellCycInd{c};
else
    C = setdiff(respCellCycInd{c},unique(cell2mat(respCellCycInd(1:c-1))));
end
nc = 25;
sp1 = 5;
sp2 = 5;

% sort by avg anti resp
% [sort_order sort_ind] = sort(squeeze(mean(tc_allcells(pre_event_frames:end,C),1)));

% sort by modulation index
[sort_order sort_ind] = sort(mi_VsubA_early(C));
% [sort_order sort_ind] = sort(mi_VsubA_late(C));

if length(C) >= nc
    % biggest vis effect
% sortCell2Plot = C(fliplr(sort_ind(end-nc+1:end)));
    % biggest aud effect
sortCell2Plot = C(fliplr(sort_ind(1:nc)));  

exCellsEarlyAudMod = figure;
suptitle({[num2str(c) ' cyc resp cells'];'early win V sub A, most modulated by aud trials'})
for i = 1:nc
    subplot(sp1,sp2,i)
    plot(ttMs(20:end),tc_vis_allcells(20:end,sortCell2Plot(i)),'g')
    hold on
    plot(ttMs(20:end),tc_aud_allcells(20:end,sortCell2Plot(i)),'k')
    hold on
    xlim([ttMs(20) ttMs(end)])
    ylim([-0.05 max(tc_vis_allcells(20:end,sortCell2Plot(i)))+0.05]);
    title(['cell ' num2str(sortCell2Plot(i))]);
end
else
    nc = length(C);
    % biggest vis effect
%    sortCell2Plot = C(fliplr(sort_ind(end-nc+1:end)));
    % biggest aud effect
sortCell2Plot = C(fliplr(sort_ind(1:nc)));  
exCellsEarlyAudMod = figure;
suptitle({[num2str(c) ' cyc resp cells'];'early win V sub A, most modulated by aud trials'})
for i = 1:nc
    subplot(sp1,sp2,i)
    plot(ttMs(20:end),tc_vis_allcells(20:end,sortCell2Plot(i)),'g')
    hold on
    plot(ttMs(20:end),tc_aud_allcells(20:end,sortCell2Plot(i)),'k')
    hold on
    xlim([ttMs(20) ttMs(end)])
    ylim([-0.05 0.1]);
    title(['cell ' num2str(sortCell2Plot(i))]);
end 
end

figure(exCellsEarlyAudMod)
print([fnout 'press_align_exCellsEarlyModAud'], '-dpdf')

% sort by modulation index
[sort_order sort_ind] = sort(mi_VsubA_late(C));
if length(C) >= nc
%     biggest vis effect
sortCell2Plot = C(fliplr(sort_ind(end-nc+1:end)));
    % biggest aud effect
% sortCell2Plot = C(fliplr(sort_ind(1:nc)));  

exCellsLateVisMod = figure;
suptitle({[num2str(c) ' cyc resp cells'];'late win V sub A, most modulated by vis trials'})
for i = 1:nc
    subplot(sp1,sp2,i)
    plot(ttMs(20:end),tc_vis_allcells(20:end,sortCell2Plot(i)),'g')
    hold on
    plot(ttMs(20:end),tc_aud_allcells(20:end,sortCell2Plot(i)),'k')
    hold on
    xlim([ttMs(20) ttMs(end)])
    ylim([-0.05 max(tc_vis_allcells(20:end,sortCell2Plot(i)))+0.05]);
    title(['cell ' num2str(sortCell2Plot(i))]);
end
else
    nc = length(C);
    % biggest vis effect
%    sortCell2Plot = C(fliplr(sort_ind(end-nc+1:end)));
    % biggest aud effect
sortCell2Plot = C(fliplr(sort_ind(1:nc)));  
exCellsLateVisMod = figure;
suptitle({[num2str(c) ' cyc resp cells'];'late win V sub A, most modulated by vis trials'})
for i = 1:nc
    subplot(sp1,sp2,i)
    plot(ttMs(20:end),tc_vis_allcells(20:end,sortCell2Plot(i)),'g')
    hold on
    plot(ttMs(20:end),tc_aud_allcells(20:end,sortCell2Plot(i)),'k')
    hold on
    xlim([ttMs(20) ttMs(end)])
    ylim([-0.05 0.1]);
    title(['cell ' num2str(sortCell2Plot(i))]);
end 
end
figure(exCellsLateVisMod)
print([fnout 'press_align_exCellsLateModVis'], '-dpdf')