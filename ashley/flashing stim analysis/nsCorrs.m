nsCorr_visB = figure;
colormap(brewermap([],'*RdBu'))
suptitle('visual trials - baseline')
nsCorr_audB = figure;
colormap(brewermap([],'*RdBu'))
suptitle('auditory trials - baseline')
nsCorr_vis1 = figure;
colormap(brewermap([],'*RdBu'))
suptitle('visual trials - 1st 3 cycles')
nsCorr_aud1 = figure;
colormap(brewermap([],'*RdBu'))
suptitle('auditory trials - 1st 3 cycles')
nsCorr_vis2 = figure;
colormap(brewermap([],'*RdBu'))
suptitle('visual trials - last 3+ cycles')
nsCorr_aud2 = figure;
colormap(brewermap([],'*RdBu'))
suptitle('auditory trials - last 3+ cycles')

i = 1;

meanNsCorr_vis = nan(3,size(expt,2));
meanNsCorr_aud = nan(3,size(expt,2));
errNsCorr_vis = nan(3,size(expt,2));
errNsCorr_aud = nan(3,size(expt,2));

absNsCorr_vis = nan(3,size(expt,2));
absNsCorr_aud = nan(3,size(expt,2));
errAbsNsCorr_vis = nan(3,size(expt,2));
errAbsNsCorr_aud = nan(3,size(expt,2));

nsc_vis_dir = nan(size(expt,2),3,4);
errnsc_vis_dir = nan(size(expt,2),3,4);
nsc_aud_dir = nan(size(expt,2),3,4);
errnsc_aud_dir = nan(size(expt,2),3,4);

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
    select_ind = cell(4,1);
    cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
    cell_ind = intersect(mouse(imouse).expt(iexp).cells(1).ind,cell_ind);
    for iOri = 1:4
        [C select_ind{iOri}] = intersect(cell_ind,mouse(imouse).expt(iexp).cells(iOri+1).ind);
    end
    trL_vis = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).tcyc;
    trL_aud = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).tcyc;
    cycs = unique(trL_vis);
    rVis = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(:,cell_ind,:);
    rAud = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(:,cell_ind,:);
    cycTime = mouse(imouse).expt(iexp).info.cyc_time;
    maxTrialLength = mouse(imouse).expt(iexp).post_event_frames;
    intVisB = [];
    intAudB = [];
    intVis1 = [];
    intAud1 = [];
    intVis2 = [];
    intAud2 = [];
    if length(cell_ind)>1
    for icyc = 1:length(cycs)
        if cycs(icyc)*cycTime <= maxTrialLength;
        ind = find(trL_vis == cycs(icyc));
        intVisB = cat(1,intVisB,squeeze(trapz(rVis(1:pre_event_frames,:,ind)))');
        intVis1 = cat(1,intVis1,squeeze(trapz(rVis(pre_event_frames:pre_event_frames+(3*cycTime),:,ind)))');
        intVis2 = cat(1,intVis2,squeeze(trapz(rVis(pre_event_frames+(3*cycTime):pre_event_frames+(cycs(icyc)*cycTime),:,ind)))');
        ind = find(trL_aud == cycs(icyc));
        intAudB = cat(1,intAudB,squeeze(trapz(rAud(1:pre_event_frames,:,ind)))');
        intAud1 = cat(1,intAud1,squeeze(trapz(rAud(pre_event_frames:pre_event_frames+(3*cycTime),:,ind)))');
        intAud2 = cat(1,intAud2,squeeze(trapz(rAud(pre_event_frames+(3*cycTime):pre_event_frames+(cycs(icyc)*cycTime),:,ind)))');
        end
    end
    end
    nsc_visB = tril(corrcoef(intVisB),-1);
    nsc_vis1 = tril(corrcoef(intVis1),-1);
    nsc_vis2 = tril(corrcoef(intVis2),-1);
    nsc_audB = tril(corrcoef(intAudB),-1);
    nsc_aud1 = tril(corrcoef(intAud1),-1);
    nsc_aud2 = tril(corrcoef(intAud2),-1);
    
    if imouse == 1 & iexp == 2
        nsCorr_exData = figure;
        colormap(brewermap([],'*RdBu'))
        suptitle([mouse(imouse).expt(iexp).mouse_name '-' mouse(imouse).expt(iexp).date])
        subplot(3,2,1)
        imagesc(nsc_visB);
        axis square
        caxis([-1 1])
        title('visual trials')
        ylabel('baseline')
        subplot(3,2,3)
        imagesc(nsc_vis1);
        axis square
        caxis([-1 1])
        ylabel('< 3 cyc')
        subplot(3,2,5)
        imagesc(nsc_vis2);
        axis square
        caxis([-1 1])
        ylabel('> 3 cyc')
        subplot(3,2,2)
        imagesc(nsc_audB);
        axis square
        caxis([-1 1])
        title('auditory trials')
        subplot(3,2,4)
        imagesc(nsc_aud1);
        axis square
        caxis([-1 1])
        subplot(3,2,6)
        imagesc(nsc_aud2);
        axis square
        caxis([-1 1])
                
        
    end
    
    for iOri = 1:4;
        ns_temp = nsc_visB(select_ind{iOri},:);
        nsc_vis_dir(i,1,iOri) = mean(ns_temp(ns_temp ~=0));        
        errnsc_vis_dir(i,1,iOri) = std(ns_temp(ns_temp ~=0))/sqrt(length(ns_temp(ns_temp ~=0)));
        ns_temp = nsc_vis1(select_ind{iOri},:);
        nsc_vis_dir(i,2,iOri) = mean(ns_temp(ns_temp ~=0));        
        errnsc_vis_dir(i,2,iOri) = std(ns_temp(ns_temp ~=0))/sqrt(length(ns_temp(ns_temp ~=0)));
        ns_temp = nsc_vis2(select_ind{iOri},:);
        nsc_vis_dir(i,3,iOri) = mean(ns_temp(ns_temp ~=0));        
        errnsc_vis_dir(i,3,iOri) = std(ns_temp(ns_temp ~=0))/sqrt(length(ns_temp(ns_temp ~=0)));
        
        ns_temp = nsc_audB(select_ind{iOri},:);
        nsc_aud_dir(i,1,iOri) = mean(ns_temp(ns_temp ~=0));        
        errnsc_aud_dir(i,1,iOri) = std(ns_temp(ns_temp ~=0))/sqrt(length(ns_temp(ns_temp ~=0)));
        ns_temp = nsc_aud1(select_ind{iOri},:);
        nsc_aud_dir(i,2,iOri) = mean(ns_temp(ns_temp ~=0));        
        errnsc_aud_dir(i,2,iOri) = std(ns_temp(ns_temp ~=0))/sqrt(length(ns_temp(ns_temp ~=0)));
        ns_temp = nsc_aud2(select_ind{iOri},:);
        nsc_aud_dir(i,3,iOri) = mean(ns_temp(ns_temp ~=0));        
        errnsc_aud_dir(i,3,iOri) = std(ns_temp(ns_temp ~=0))/sqrt(length(ns_temp(ns_temp ~=0)));
    end
    
    if length(cell_ind) > 10
    meanNsCorr_vis(1,i) = mean(nsc_visB(nsc_visB ~= 0));
    errNsCorr_vis(1,i) = std(nsc_visB(nsc_visB ~= 0))/sqrt(length(nsc_visB(nsc_visB ~= 0)));
    meanNsCorr_vis(2,i) = mean(nsc_vis1(nsc_vis1 ~= 0));
    errNsCorr_vis(2,i) = std(nsc_vis1(nsc_vis1 ~= 0))/sqrt(length(nsc_vis1(nsc_vis1 ~= 0)));
    meanNsCorr_vis(3,i) = mean(nsc_vis2(nsc_vis2 ~= 0));
    errNsCorr_vis(3,i) = std(nsc_vis2(nsc_vis2 ~= 0))/sqrt(length(nsc_vis2(nsc_vis2 ~= 0)));
    
    meanNsCorr_aud(1,i) = mean(nsc_audB(nsc_audB ~= 0));
    errNsCorr_aud(1,i) = std(nsc_audB(nsc_audB ~= 0))/sqrt(length(nsc_audB(nsc_audB ~= 0)));
    meanNsCorr_aud(2,i) = mean(nsc_aud1(nsc_aud1 ~= 0));
    errNsCorr_aud(2,i) = std(nsc_aud1(nsc_aud1 ~= 0))/sqrt(length(nsc_aud1(nsc_aud1 ~= 0)));
    meanNsCorr_aud(3,i) = mean(nsc_aud2(nsc_aud2 ~= 0)); 
    errNsCorr_aud(3,i) = std(nsc_aud2(nsc_aud2 ~= 0))/sqrt(length(nsc_aud2(nsc_aud2 ~= 0)));
    
    absNsCorr_vis(1,i) = mean(abs(nsc_visB(nsc_visB ~= 0)));
    errAbsNsCorr_vis(1,i) = std(abs(nsc_visB(nsc_visB ~= 0)))/sqrt(length(nsc_visB(nsc_visB ~= 0)));
    absNsCorr_vis(2,i) = mean(abs(nsc_vis1(nsc_vis1 ~= 0)));
    errAbsNsCorr_vis(2,i) = std(abs(nsc_vis1(nsc_vis1 ~= 0)))/sqrt(length(nsc_vis1(nsc_vis1 ~= 0)));
    absNsCorr_vis(3,i) = mean(abs(nsc_vis2(nsc_vis2 ~= 0)));
    errAbsNsCorr_vis(3,i) = std(abs(nsc_vis2(nsc_vis2 ~= 0)))/sqrt(length(nsc_vis2(nsc_vis2 ~= 0)));
    
    absNsCorr_aud(1,i) = mean(abs(nsc_audB(nsc_audB ~= 0)));
    errAbsNsCorr_aud(1,i) = std(abs(nsc_audB(nsc_audB ~= 0)))/sqrt(length(nsc_audB(nsc_audB ~= 0)));
    absNsCorr_aud(2,i) = mean(abs(nsc_aud1(nsc_aud1 ~= 0)));
    errAbsNsCorr_aud(2,i) = std(abs(nsc_aud1(nsc_aud1 ~= 0)))/sqrt(length(nsc_aud1(nsc_aud1 ~= 0)));
    absNsCorr_aud(3,i) = mean(abs(nsc_aud2(nsc_aud2 ~= 0))); 
    errAbsNsCorr_aud(3,i) = std(abs(nsc_aud2(nsc_aud2 ~= 0)))/sqrt(length(nsc_aud2(nsc_aud2 ~= 0)));
    end
    
    % heatmaps
    figure(nsCorr_visB);
    subplot(n,n2,i)
    imagesc(nsc_visB);
    axis square
    caxis([-1 1])
    figure(nsCorr_vis1);
    subplot(n,n2,i)
    imagesc(nsc_vis1);
    axis square
    caxis([-1 1])
    figure(nsCorr_vis2);
    subplot(n,n2,i)
    imagesc(nsc_vis2);
    axis square
    caxis([-1 1])
    figure(nsCorr_audB);
    subplot(n,n2,i)
    imagesc(nsc_audB);
    axis square
    caxis([-1 1])
    figure(nsCorr_aud1);
    subplot(n,n2,i)
    imagesc(nsc_aud1);
    axis square
    caxis([-1 1])
    figure(nsCorr_aud2);
    subplot(n,n2,i)
    imagesc(nsc_aud2);
    axis square
    caxis([-1 1])
    
    i = i+1;
    end
end
    
meanNsCorrPlot = figure;
subplot(3,2,1)
errorbar(nanmean(meanNsCorr_vis,2),nanmean(errNsCorr_vis,2),'go-');
hold on
errorbar(nanmean(meanNsCorr_aud,2),nanmean(errNsCorr_aud,2),'ko-');
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'base'; '< 3cyc'; '> 3cyc'};
xlabel('time period')
ylabel('mean ns corr')

subplot(3,2,2)
errorbar(nanmean(absNsCorr_vis,2),nanmean(errAbsNsCorr_vis,2),'go-');
hold on
errorbar(nanmean(absNsCorr_aud,2),nanmean(errAbsNsCorr_aud,2),'ko-');
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'base'; '< 3cyc'; '> 3cyc'};
xlabel('time period')
ylabel('mean abs val ns corr')

for iOri = 1:4
    subplot(3,2,iOri+2)
    errorbar(nanmean(squeeze(nsc_vis_dir(:,:,iOri)),1),nanmean(squeeze(errnsc_vis_dir(:,:,iOri)),1),'go-');
    hold on
    errorbar(nanmean(squeeze(nsc_aud_dir(:,:,iOri)),1),nanmean(squeeze(errnsc_aud_dir(:,:,iOri)),1),'ko-');
    ax = gca;
    ax.XTick = [1 2 3];
    ax.XTickLabel = {'base'; '< 3cyc'; '> 3cyc'};
    xlabel('time period')
    ylabel('mean ns corr')
    title([mouse(imouse).expt(iexp).cells(iOri+1).name ' deg selective'])
end
for iplot = 1:6
    figure(meanNsCorrPlot);
    subplot(3,2,iplot)
    hold on
    ylim([-0.05 0.3])
end
    

figtitles = {'baseline';'1st 3 cycs';'end of trial'};
meanNsCorrScatter = figure;
for iplot = 1:3
subplot(2,3,iplot)
errbar(meanNsCorr_vis(iplot,:),meanNsCorr_aud(iplot,:),errNsCorr_aud(iplot,:),'k');
hold on
errbar(meanNsCorr_vis(iplot,:),meanNsCorr_aud(iplot,:),errNsCorr_vis(iplot,:),'horiz','k');
hold on
plot(meanNsCorr_vis(iplot,:),meanNsCorr_aud(iplot,:),'ko');
hold on
plot(nanmean(meanNsCorr_vis(iplot,:),2),nanmean(meanNsCorr_aud(iplot,:),2),'r.','markersize', 20)
hold on
axis square
plot([-1:0.1:1],[-1:0.1:1],'k:')
xlim([-0.2 0.5])
ylim([-0.2 0.5])
title(figtitles{iplot})
if iplot ==  1
    ylabel('mean ns corr')
else
    ylabel('auditory')
end
xlabel('visual')

subplot(2,3,iplot+3)
errbar(absNsCorr_vis(iplot,:),absNsCorr_aud(iplot,:),errAbsNsCorr_aud(iplot,:),'k');
hold on
errbar(absNsCorr_vis(iplot,:),absNsCorr_aud(iplot,:),errAbsNsCorr_vis(iplot,:),'horiz','k');
hold on
plot(absNsCorr_vis(iplot,:),absNsCorr_aud(iplot,:),'ko');
hold on
plot(nanmean(absNsCorr_vis(iplot,:),2),nanmean(absNsCorr_aud(iplot,:),2),'r.','markersize', 20)
hold on
axis square
plot([-1:0.1:1],[-1:0.1:1],'k:')
xlim([-0.2 0.5])
ylim([-0.2 0.5])
title(figtitles{iplot})
if iplot ==  1
    ylabel('abs val ns corr')
else
    ylabel('auditory')
end
xlabel('visual')
end    

%% save
figure(meanNsCorrPlot)
print([fnout 'meanNsCorrPlot' datasetStr '.pdf'], '-dpdf')
figure(meanNsCorrScatter)
print([fnout 'meanNsCorrScatter' datasetStr '.pdf'], '-dpdf')
figure(nsCorr_exData);
print([fnout 'nsCorrExData' datasetStr '.pdf'], '-dpdf')