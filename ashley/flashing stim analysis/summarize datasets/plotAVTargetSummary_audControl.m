function plotAVCatchSummary(datasetStr,cellsInd,cellsOnly)
close all
respCellsInd = 14;

eval(['awFSAVdatasets' datasetStr]);

titleStr = datasetStr;
titleStr = titleStr(2:end);

rc = behavConstsAV;
if strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' datasetStr];
else
    dataGroup = [];
end
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
% mouse_ind = find(intersect(cell2mat({av.mouse}),values));
if  cellsOnly
load(fullfile(rc.caOutputDir,dataGroup, 'cells only',[mouse_str '_CaSummary' datasetStr '.mat']));
titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name]; 
fnout = fullfile(rc.caOutputDir,dataGroup, 'cells only',[titleStr '_' mouse_str]);
else
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' datasetStr '.mat']));
titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name]; 
fnout = fullfile(rc.caOutputDir,dataGroup, [titleStr '_' mouse_str '_gr' num2str(respCellsInd)]);
end

pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
cycTime = mouse(1).expt(1).info(1).cyc_time;
cycTimeMs = mouse(1).expt(1).info(1).cyc_time_ms;

targetAlign = 2;
hits = 1;
vis = 1;
aud = 2;

nexp = 0;
for imouse  = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
n = ceil(sqrt(nexp+1));
if (n^2)-n > nexp+1
    n2 = n-1;
else
    n2 = n;
end

%% set params for figures
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
% set(0,'DefaultaxesFontSize', 16)

tt = -pre_event_frames:post_event_frames-1;
ttMs = tt/(cycTime/cycTimeMs);
baseStimFrames = -(floor(pre_event_frames/cycTime)*cycTime):cycTime:0;

%% figure handles
heatmapFig = figure;
suptitle('avg all trials - sort across all trials')
colormap(brewermap(1001,'*RdBu'))
tcFig = figure;

%% get target response for each cell, each expt
trialLengthInd = 1:length(tt);

tc_vis = [];
tc_aud = [];
tc_all = [];

tc_vis_trials = [];
tc_aud_trials = [];
tc_all_trials = [];

tcyc_vis = [];
tcyc_aud = [];
tcyc_all = [];

for imouse  = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2) 
        cell_ind = mouse(imouse).expt(iexp).cells(respCellsInd).ind;
        
        vistemp = mouse(imouse).expt(iexp).align(targetAlign).av(vis).outcome(hits).resp(trialLengthInd,cell_ind,:);
        audtemp = mouse(imouse).expt(iexp).align(targetAlign).av(aud).outcome(hits).resp(trialLengthInd,cell_ind,:);
        alltemp = cat(3,vistemp,audtemp);
        visTtemp = mouse(imouse).expt(iexp).align(targetAlign).av(vis).outcome(hits).tcyc;
        audTtemp = mouse(imouse).expt(iexp).align(targetAlign).av(aud).outcome(hits).tcyc;
        allTtemp = cat(2,visTtemp,audTtemp);
        
        tc_vis_trials = cat(2,tc_vis,squeeze(mean(vistemp,2)));
        tc_aud_trials = cat(2,tc_aud,squeeze(mean(audtemp,2)));
        tc_all_trials = cat(2,tc_all,squeeze(mean(alltemp,2)));
        
        tc_vis = cat(2,tc_vis,mean(vistemp,3));
        tc_aud = cat(2,tc_aud,mean(audtemp,3));
        tc_all  = cat(2,tc_all,mean(alltemp,3));
        
        tcyc_vis = cat(2,tcyc_vis,visTtemp);
        tcyc_aud = cat(2,tcyc_aud,audTtemp);
        tcyc_all = cat(2,tcyc_all,allTtemp);
    end
end

nCells = size(tc_all,2);

%% zero all time-courses
tc_vis = bsxfun(@minus,tc_vis,mean(tc_vis(pre_win,:)));
tc_aud = bsxfun(@minus,tc_aud,mean(tc_aud(pre_win,:)));
tc_all = bsxfun(@minus,tc_all,mean(tc_all(pre_win,:)));

%% heatmap
[tcSort tcSortInd] = sort(mean(tc_all(pre_event_frames:end,:),1));
tcSortInd = fliplr(tcSortInd);

figure(heatmapFig)
subplot(2,2,1)
allTrials_sort = imagesc(tc_all(20:end,tcSortInd)');
title('all trials')
subplot(2,2,2)
visTrials_sort = imagesc(tc_vis(20:end,tcSortInd)');
title('visual trials')
subplot(2,2,3)
audTrials_sort = imagesc(tc_aud(20:end,tcSortInd)');
title('auditory trials')
subplot(2,2,4)
subtract_sort = imagesc(tc_vis(20:end,tcSortInd)'-tc_aud(20:end,tcSortInd)');
title('vis - aud')

for iplot = 1:4
    subplot(2,2,iplot)
ax = gca;
ax.XTick = [find(ceil(ttMs(20:end)) == 0) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000/2) find(ttMs(20:end) == floor(max(ttMs(20:end)/1000))*1000)];
ax.XTickLabel = [0 floor(max(ttMs(20:end)/1000))/2 floor(max(ttMs(20:end)/1000))];
xlabel('time (s)')
ylabel('cells')
clim_vis = [-.2 0.2];
caxis(clim_vis)
colorbar
axis square
end

print([fnout '_tarRespHeatmap.pdf'],'-dpdf','-fillpage')
%% average time-course
tc_all_mean = mean(tc_all,2);
tc_all_ste = std(tc_all,[],2)/sqrt(nCells);
tc_vis_mean = mean(tc_vis,2);
tc_vis_ste = std(tc_vis,[],2)/sqrt(nCells);
tc_aud_mean = mean(tc_aud,2);
tc_aud_ste = std(tc_aud,[],2)/sqrt(nCells);


figure(tcFig)
subplot(2,2,1)
shadedErrorBar(ttMs,tc_all_mean,tc_all_ste,'k')
title('all trials, align to target')

subplot(2,2,2)
shadedErrorBar(ttMs,tc_vis_mean,tc_vis_ste,'g')
hold on
shadedErrorBar(ttMs,tc_aud_mean,tc_aud_ste,'k')
title('vis(g) and aud(blk) trials, align to target')
for iplot = 1:2
    subplot(2,2,iplot)
    hold on
    xlim([-10 20]/(cycTime/cycTimeMs))
    ylim([-0.01 0.05])
    xlabel('time(ms)')
    ylabel('df/f')
    vline(baseStimFrames,':k')
    hold on
    vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
    hold on
    vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
    axis square
end

%% scatter response
tc_vis_resp = mean(tc_vis(trans_win,:),1)-mean(tc_vis(pre_win,:),1);
tc_aud_resp = mean(tc_aud(trans_win,:),1)-mean(tc_aud(pre_win,:),1);

[h,p] = ttest(tc_vis_resp,tc_aud_resp,'alpha',0.05/(nCells));

figure(tcFig)
subplot(2,2,3)
scatter(tc_vis_resp,tc_aud_resp,50,'k.')
hold on
errorbarxy(mean(tc_vis_resp),mean(tc_aud_resp),std(tc_vis_resp)/sqrt(nCells),std(tc_aud_resp)/sqrt(nCells),{'ro','r','r'})
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('visual')
ylabel('auditory')
title(['p=' num2str(p) ', ' num2str(nCells) '  cells'])

%% cdf response
figure(tcFig)
subplot(2,2,4)
h = cdfplot(tc_vis_resp);
h.Color = 'g';
hold on
h = cdfplot(tc_aud_resp);
h.Color = 'k';
[h p] = kstest2(tc_vis_resp,tc_aud_resp);
xlim([-0.05 0.05])
title(['p= ' num2str(p)])
axis square

print([fnout '_tarRespTC'],'-dpdf','-fillpage')
%% response by trial length

end