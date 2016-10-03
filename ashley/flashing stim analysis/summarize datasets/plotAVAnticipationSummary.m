function plotAVAnticipationSummary(datasetStr,cellsInd,cellsOnly)
% takes mouse FSAV Ca structure and plots responses to anticipation phase
% of task
close all
respCellsInd = 14;
av = behavParamsAV_naive;
eval(['awFSAVdatasets' datasetStr])
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
else
    titleStr = titleStr(2:end);
end

rc = behavConstsAV;
if strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' datasetStr];
else
    dataGroup = [];
end
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
mouse_ind = find(intersect(cell2mat({av.mouse}),values));
if cellsOnly
load(fullfile(rc.caOutputDir,dataGroup, 'cells only',[mouse_str '_CaSummary' datasetStr '.mat']));
titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name];
fnout = fullfile(rc.caOutputDir, dataGroup, 'cells only', [titleStr '_' mouse_str]); %% maybe lose mouse_str
else
load(fullfile(rc.caOutputDir,dataGroup, 'cells only',[mouse_str '_CaSummary' datasetStr '.mat']));
titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name];
fnout = fullfile(rc.caOutputDir, dataGroup, [titleStr '_' mouse_str]); %% maybe lose mouse_str
end


pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
minTrialLengthFrames = mouse(1).expt(1).info.minTrialLengthFrames;
ialign = 1;
visual = 1;
auditory = 2;
hits = 1;
cycTime = mouse(1).expt(1).info(1).cyc_time;
cycTimeMs = mouse(1).expt(1).info(1).cyc_time_ms;
nexp = 0;
for imouse = 1:size(mouse,2)
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

tt = -pre_event_frames:minTrialLengthFrames-1;
ttMs = tt/(cycTime/cycTimeMs);
baseStimFrames = 0:cycTime:minTrialLengthFrames-1;

%%
antiRespAllCyc
if ~strcmp(datasetStr,'_V1')
antiRespEarlyLateQuant
end
%%
depOnPrevTrialType

%%


%% plot avg trace for each expt and avg across all expt - success and responsive cells only; plot std and variance over time
respAvsVFig = figure;
suptitle({titleStr, 'avg response - Visual: Green; Auditory: black'})
stdAvsVFig = figure;
suptitle({titleStr, 'stdev response - Visual: Green; Auditory: black'})
i = 1;
resp_vis = [];
resp_aud = [];
resp_all = [];
std_vis = [];
std_aud = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        
        figure(respAvsVFig);
        subplot(n,n2,i)
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),2),3), 'g');
        hold on
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:),2),3), 'k');
        hold on
        xlim([-10 minTrialLengthFrames])
        ylim([-0.01 0.1]);       
%         vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
%         vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
        vline(baseStimFrames,':k')
        resp_vis = cat(2, resp_vis, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),3));
        resp_aud = cat(2, resp_aud, mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:),3));
        resp_all = cat(2, resp_all, mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:)),3));
        
        title({mouse(imouse).expt(iexp).date; [' n = ' num2str(length(cell_ind)) ' cells']})
        
        figure(stdAvsVFig);
        subplot(n,n2,i)
        plot(tt,std(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),2),[],3), 'g');
        hold on
        plot(tt,std(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:),2),[],3), 'k');
        hold on
        xlim([-10 minTrialLengthFrames])
        ylim([-0.01 0.1]);
%         vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
%         vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
        vline(baseStimFrames,':k')
        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
        std_vis = cat(2, std_vis, std(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),[],3));
        std_aud = cat(2, std_aud, std(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:),[],3));
        
       i = i+1;
    end
end
figure(respAvsVFig);
subplot(n,n2,i)
plot(tt, nanmean(resp_vis,2), 'g')
hold on
plot(tt, nanmean(resp_aud,2), 'k')
% vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
% vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
xlim([-10 minTrialLengthFrames])
ylim([-0.01 0.03]);
vline(baseStimFrames,':k')
ylim([-0.01 0.05]);
title(['All cells; n = ' num2str(size(resp_vis,2))])

figure(stdAvsVFig);
subplot(n,n2,i)
plot(tt,std(resp_vis,[],2), 'g');
hold on
plot(tt,std(resp_aud,[],2), 'k');
hold on
%         vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
%         vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
vline(baseStimFrames,':k')
xlim([-10 minTrialLengthFrames])
ylim([-0.01 0.05]);
title(['All cells; n = ' num2str(size(resp_vis,2))])        

figure(respAvsVFig);
print([fnout 'press_align_TCs' datasetStr '.pdf'], '-dpdf')
figure(stdAvsVFig);
print([fnout 'press_align_stdTCs' datasetStr '.pdf'], '-dpdf')


%% sort all neurons by avg resp across both trial types. 

nBins = 10;
ncells = size(resp_all,2);
peakR_all = squeeze(mean(resp_all(pre_event_frames:end,:),1));
[peakR_all_sorted peakR_all_sort_ind] = sort(peakR_all);

bin_edges = floor(linspace(0,ncells,nBins+1));

resp_vis_bin = zeros(size(resp_all,1),nBins);
resp_aud_bin = zeros(size(resp_all,1),nBins);
errV_bin = zeros(size(resp_all,1),nBins);
errA_bin = zeros(size(resp_all,1),nBins);
for ibin = 1:nBins
    ind = bin_edges(ibin)+1:bin_edges(ibin+1);
    resp_vis_bin(:,ibin) = mean(resp_vis(:,peakR_all_sort_ind(ind)),2);
    resp_aud_bin(:,ibin) = mean(resp_aud(:,peakR_all_sort_ind(ind)),2);
    errV_bin(:,ibin) = std(resp_vis(:,peakR_all_sort_ind(ind)),[],2)/sqrt(length(ind));
    errA_bin(:,ibin) = std(resp_aud(:,peakR_all_sort_ind(ind)),[],2)/sqrt(length(ind));
end

colorsA = brewermap(nBins+5,'Greys');
colorsV = brewermap(nBins+5,'YlGn');

tcMsBinFig = figure;
for ibin = 1:nBins
    subplot(5,2,ibin)
    errbar(ttMs,resp_aud_bin(:,ibin),errA_bin(:,ibin),'k')
    hold on
    errbar(ttMs,resp_vis_bin(:,ibin),errV_bin(:,ibin),'g')
    hold on
    binsA(ibin) = plot(ttMs,resp_aud_bin(:,ibin),'k');
    hold on
    binsV(ibin) = plot(ttMs,resp_vis_bin(:,ibin),'g');xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
    xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
%     ylim([-0.01 0.05]);
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    title(['bin ' num2str(ibin)])
    xlabel('ms')
    ylabel('dF/F')
end


print([fnout 'press_align_TCsMsBinned' datasetStr '.pdf'], '-dpdf')

colorsMI = brewermap(nBins,'Reds');

respMod_VsubA = resp_vis_bin-resp_aud_bin;
tcModBinFig = figure;
set(gca,'ColorOrder',colorsMI,'NextPlot','replacechildren');
col = get(gca,'ColorOrder');
hold all
% colormap(brewermap([],'Reds'))
F = plot(ttMs,respMod_VsubA);
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
%     ylim([-0.01 0.05]);
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title('vis-aud, binned')
xlabel('ms')
ylabel('dF/F')
axis square
legend(cellfun(@num2str,num2cell(1:nBins),'unif',false),'Location','northeastoutside')
print([fnout 'press_align_MIsMsBinned' datasetStr '.pdf'], '-dpdf')

%% mean modulation bw aud and vis trials vs. resp to target
% resp_mod = resp_vis-resp_aud;
% r_mod = mean(resp_mod(pre_event_frames:end,:),1);
% r_mod_end = mean(resp_mod(end-(cycTime*2):end,:),1);
% 
% d = [];
% for imouse = 1:size(mouse,2)
%     for iexp = 1:size(mouse(imouse).expt,2)
%         d = [d mouse(imouse).expt(iexp).visTargets];
%         d = unique(d);
%     end
% end
% 
% r_tar_vis = [];
% 
% for imouse = 1:size(mouse,2)
%     for iexp = 1:size(mouse(imouse).expt,2)
%         cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
%         cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);   
%         tempR = mean(mean(mouse(imouse).expt(iexp).align(2).av(1).outcome(1).resp(trans_win,cell_ind,:),3),1)-mean(mean(mouse(imouse).expt(iexp).align(2).av(1).outcome(1).resp(pre_win,cell_ind,:),3),1);
%         r_tar_vis = cat(2, r_tar_vis, tempR);
%     end
% end
% 
% figure;
% for ibin = 1:nBins
%     ind = bin_edges(ibin)+1:bin_edges(ibin+1);
%     subplot(5,2,ibin)
%     scatter(r_mod(peakR_all_sort_ind(ind)),r_tar_vis(peakR_all_sort_ind(ind)),50,'k.')
%     hold on    
% %     plot([-.1:0.01:.1],[-.1:0.01:.1],'k--')
%     title(['bin ' num2str(ibin)])
%     xlabel('antic. MI')
%     ylabel('target resp')
%     xlim([-0.03 0.03])
%     axis square
% end
% 
% figure;
% for ibin = 1:nBins
%     ind = bin_edges(ibin)+1:bin_edges(ibin+1);
%     subplot(5,2,ibin)
%     scatter(r_mod(peakR_all_sort_ind(ind)),r_tar_vis(peakR_all_sort_ind(ind)),50,'k.')
%     hold on    
% %     plot([-.1:0.01:.1],[-.1:0.01:.1],'k--')
%     title(['bin ' num2str(ibin)])
%     xlabel('antic. MI;end')
%     ylabel('target resp')
%     xlim([-0.03 0.03])
%     ylim([0 0.05])
%     axis square
% end
%%
avgCellRespHeatMap
% if cellsInd == 1 | cellsInd == 13 | cellsInd == 12
% nsCorrs
% end
% %% 
% sigRespByCyc
%% plot with time in ms

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

ncells = size(resp_vis,2);
errV = nanstd(resp_vis,[],2)/sqrt(ncells);
errA = nanstd(resp_aud,[],2)/sqrt(ncells);

tcMsFig = figure;
subplot(3,2,1)
errbar(ttMs, nanmean(resp_aud,2),errA, 'k')
hold on
plot(ttMs,nanmean(resp_aud,2),'k')
hold on
errbar(ttMs, nanmean(resp_vis,2),errV, 'g')
hold on
plot(ttMs,nanmean(resp_vis,2),'g')
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
% if cellsInd == 2
%     ylim([-0.01 0.05]);
% else
% ylim([-0.01 0.02]);
% end
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title(['All cells; n = ' num2str(ncells)])
xlabel('ms')
ylabel('dF/F')

subplot(3,2,2)
plot(ttMs,std(resp_vis,[],2), 'g');
hold on
plot(ttMs,std(resp_aud,[],2), 'k');
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
% ylim([-0.01 0.05]);
title(['stdev; n = ' num2str(ncells)])   
xlabel('ms')
ylabel('dF/F')



%*** timecourse of top and bottom responders

%top vis responders tc
percentTop = .1;
nTop = ceil(size(tc_vis_all,2)*percentTop);
topVis_ind = peakR_vis_sort_ind(end-nTop-1:end);
topAud_ind = peakR_aud_sort_ind(end-nTop-1:end);

errV = nanstd(resp_vis(:,topVis_ind),[],2)/sqrt(nTop);
errA = nanstd(resp_aud(:,topVis_ind),[],2)/sqrt(nTop);

figure(tcMsFig);
subplot(3,2,3)
errbar(ttMs, nanmean(resp_aud(:,topVis_ind),2),errA, 'k')
hold on
plot(ttMs,nanmean(resp_aud(:,topVis_ind),2),'k')
hold on
errbar(ttMs, nanmean(resp_vis(:,topVis_ind),2),errV, 'g')
hold on
plot(ttMs,nanmean(resp_vis(:,topVis_ind),2),'g')
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
ylim([-0.01 0.1]);
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title(['top ' num2str(percentTop*100) '% vis responders; n = ' num2str(nTop)])
xlabel('ms')
ylabel('dF/F')

%top aud resp
errV = nanstd(resp_vis(:,topAud_ind),[],2)/sqrt(nTop);
errA = nanstd(resp_aud(:,topAud_ind),[],2)/sqrt(nTop);

figure(tcMsFig);
subplot(3,2,4)
errbar(ttMs, nanmean(resp_aud(:,topAud_ind),2),errA, 'k')
hold on
plot(ttMs,nanmean(resp_aud(:,topAud_ind),2),'k')
hold on
errbar(ttMs, nanmean(resp_vis(:,topAud_ind),2),errV, 'g')
hold on
plot(ttMs,nanmean(resp_vis(:,topAud_ind),2),'g')
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
ylim([-0.01 0.1]);
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title(['top ' num2str(percentTop*100) '% aud responders; n = ' num2str(nTop)])
xlabel('ms')
ylabel('dF/F')

%bottom vis responders tc
botVis_ind = peakR_vis_sort_ind(1:nTop);
botAud_ind = peakR_aud_sort_ind(1:nTop);

errV = nanstd(resp_vis(:,botVis_ind),[],2)/sqrt(nTop);
errA = nanstd(resp_aud(:,botVis_ind),[],2)/sqrt(nTop);

figure(tcMsFig);
subplot(3,2,5)
errbar(ttMs, nanmean(resp_aud(:,botVis_ind),2),errA, 'k')
hold on
plot(ttMs,nanmean(resp_aud(:,botVis_ind),2),'k')
hold on
errbar(ttMs, nanmean(resp_vis(:,botVis_ind),2),errV, 'g')
hold on
plot(ttMs,nanmean(resp_vis(:,botVis_ind),2),'g')
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
ylim([-0.04 0.01]);
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title(['bottom ' num2str(percentTop*100) '% vis responders; n = ' num2str(nTop)])
xlabel('ms')
ylabel('dF/F')

%bot aud resp
errV = nanstd(resp_vis(:,botAud_ind),[],2)/sqrt(nTop);
errA = nanstd(resp_aud(:,botAud_ind),[],2)/sqrt(nTop);

figure(tcMsFig);
subplot(3,2,6)
errbar(ttMs, nanmean(resp_aud(:,botAud_ind),2),errA, 'k')
hold on
plot(ttMs,nanmean(resp_aud(:,botAud_ind),2),'k')
hold on
errbar(ttMs, nanmean(resp_vis(:,botAud_ind),2),errV, 'g')
hold on
plot(ttMs,nanmean(resp_vis(:,botAud_ind),2),'g')
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
ylim([-0.04 0.01]);
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title(['bottom ' num2str(percentTop*100) '% aud responders; n = ' num2str(nTop)])
xlabel('ms')
ylabel('dF/F')

print([fnout 'press_align_TCsMs' datasetStr '.pdf'], '-dpdf')

%% mean across all trials all cells
errAll = std(resp_all,[],2)/sqrt(size(resp_all,2));
allTrialsFig = figure;
suptitle(titleStr)
subplot(3,2,1)
shadedErrorBar(ttMs,mean(resp_all,2),errAll,'k')
hold on
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
ylim([-0.01 0.05]);
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title(['All cells; n = ' num2str(ncells)])
xlabel('ms')
ylabel('dF/F')

%A vs V trials figure
errVis = std(resp_vis,[],2)/sqrt(size(resp_vis,2));
errAud = std(resp_aud,[],2)/sqrt(size(resp_aud,2));
AVoriTCFig = figure;
suptitle(titleStr)
subplot(3,2,1)
shadedErrorBar(ttMs,mean(resp_vis,2),errVis,'g')
hold on
shadedErrorBar(ttMs,mean(resp_aud,2),errAud,'k')
hold on
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
ylim([-0.01 0.05]);
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title(['All cells; n = ' num2str(ncells)])
xlabel('ms')
ylabel('dF/F')

%A vs V standard deviation figure
errStdV = std(std_vis,[],2)/sqrt(size(std_vis,2));
errStdA = std(std_aud,[],2)/sqrt(size(std_aud,2));
AVoriSTDFig = figure;
suptitle([titleStr ', stdev TC']);
subplot(3,2,1)
shadedErrorBar(ttMs,mean(std_vis,2),errStdV,'g')
hold on
shadedErrorBar(ttMs,mean(std_aud,2),errStdA,'k')
hold on
xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
ylim([0 0.09]);
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
title(['All cells; n = ' num2str(ncells)])
xlabel('ms')
ylabel('dF/F')


int_all_ori = cell([1,4]);
resp_all_ori = cell([1,4]);
int_vis_ori = cell([1,4]);
resp_vis_ori = cell([1,4]);
int_aud_ori = cell([1,4]);
resp_aud_ori = cell([1,4]);
std_resp_vis_ori  = cell([1,4]);
std_resp_aud_ori = cell([1,4]);

oriColMat = ['k';'g';'b';'r'];
oriStr = {'0';'45';'90';'135'};

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        for iori = 1:4
            ori_ind = intersect(mouse(imouse).expt(iexp).cells(iori+1).ind,cell_ind);
            
            tempAll = cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,ori_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,ori_ind,:));
            iAll = squeeze(mean(trapz(tempAll(pre_event_frames:end,:,:)),3));
            tcAll = mean(tempAll,3); 
            int_all_ori{iori} = cat(2, int_all_ori{iori}, iAll);
            resp_all_ori{iori} = cat(2, resp_all_ori{iori}, tcAll);

            tempVis = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,ori_ind,:);
            iVis = squeeze(mean(trapz(tempVis(pre_event_frames:end,:,:)),3));
            tcVis = mean(tempVis,3); 
            stdVis = std(tempVis,[],3);
            int_vis_ori{iori} = cat(2, int_vis_ori{iori}, iVis);
            resp_vis_ori{iori} = cat(2, resp_vis_ori{iori}, tcVis);
            std_resp_vis_ori{iori} = cat(2,std_resp_vis_ori{iori},stdVis);

            tempAud = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,ori_ind,:);
            iAud = squeeze(mean(trapz(tempAud(pre_event_frames:end,:,:)),3));
            tcAud = mean(tempAud,3); 
            stdAud = std(tempAud,[],3);
            int_aud_ori{iori} = cat(2, int_aud_ori{iori}, iAud);
            resp_aud_ori{iori} = cat(2, resp_aud_ori{iori}, tcAud);
            std_resp_aud_ori{iori} = cat(2,std_resp_aud_ori{iori},stdAud);
        end
    end
end

%cdf plot of tuned cell respones
figure(allTrialsFig);
subplot(3,2,2)
for iori = 1:4
    if ~isempty(int_all_ori{iori})
    hAll= cdfplot(int_all_ori{iori});
    hAll.Color = oriColMat(iori);
    hAll.LineStyle = '-';
    hold on
    end
end
xlabel('dF/F')
ylabel('cmlv frctn')
title('integral response - ori selective cells, all trials')
legend(oriStr)

figure(allTrialsFig)
errAll_ori = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),resp_all_ori,'unif',false);
plt = 3:6;
for iori = 1:4
    subplot(3,2,plt(iori))
    shadedErrorBar(ttMs,mean(resp_all_ori{iori},2),errAll_ori{iori},'k');
    hold on
    xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
    ylim([-0.01 0.05]);
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    title([oriStr(iori) ' slctv cells; n = ' num2str(size(resp_all_ori{iori},2))])
    xlabel('ms')
    ylabel('dF/F')
end

%cdf for aud and vis trials
figure(AVoriTCFig)
subplot(3,2,2)
for iori = 1:4
    if ~isempty(int_vis_ori{iori})
    hAll= cdfplot(int_vis_ori{iori});
    hAll.Color = oriColMat(iori);
    hAll.LineStyle = '-';
    hold on
    end
end
for iori = 1:4
    if ~isempty(int_aud_ori{iori})
    hAll= cdfplot(int_aud_ori{iori});
    hAll.Color = oriColMat(iori);
    hAll.LineStyle = '--';
    hold on
    end
end
xlabel('dF/F')
ylabel('cmlv frctn')
title('int resp - ori selective cells,vis-soid/aud-dash')
legend(oriStr)

% A vs V trials TC by ori
figure(AVoriTCFig)
errVis_ori = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),resp_vis_ori,'unif',false);
errAud_ori = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),resp_aud_ori,'unif',false);
plt = 3:6;
for iori = 1:4
    subplot(3,2,plt(iori))
    shadedErrorBar(ttMs,mean(resp_vis_ori{iori},2),errVis_ori{iori},'g');
    hold on
    shadedErrorBar(ttMs,mean(resp_aud_ori{iori},2),errAud_ori{iori},'k');
    hold on
    xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
    ylim([-0.01 0.05]);
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    title([oriStr(iori) ' slctv cells; n = ' num2str(size(resp_all_ori{iori},2))])
    xlabel('ms')
    ylabel('dF/F')
end

% stdev TC by ori, A vs V
figure(AVoriSTDFig);
errStdVis_ori = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),std_resp_vis_ori,'unif',false);
errStdAud_ori = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),std_resp_aud_ori,'unif',false);
for iori = 1:4
    subplot(3,2,plt(iori))
    shadedErrorBar(ttMs,mean(std_resp_vis_ori{iori},2),errStdVis_ori{iori},'g');
    hold on
    shadedErrorBar(ttMs,mean(std_resp_aud_ori{iori},2),errStdAud_ori{iori},'k');
    hold on
    xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
    ylim([0 0.09]);
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    title([oriStr(iori) ' slctv cells; n = ' num2str(size(resp_all_ori{iori},2))])
    xlabel('ms')
    ylabel('dF/F')
end

figure(AVoriTCFig)
print([fnout 'press_align_AVTrials_oriSort' datasetStr '.pdf'], '-dpdf')
figure(AVoriSTDFig);
print([fnout 'press_align_AVTrialsSTD_oriSort' datasetStr '.pdf'], '-dpdf')

%% scatter mean std for each cell - A vs V trials
stdAvsVScat = figure;
subplot(1,2,1)
scatter(mean(std_vis(pre_event_frames:(minTrialLengthFrames/2)+pre_event_frames,:),1),mean(std_aud(pre_event_frames:(minTrialLengthFrames/2)+pre_event_frames,:),1),'k');
hold on
scatter(mean(mean(std_vis(pre_event_frames:(minTrialLengthFrames/2)+pre_event_frames,:),1)),mean(mean(std_aud(pre_event_frames:(minTrialLengthFrames/2)+pre_event_frames,:),1)),'ro','markerfacecolor','r')
hold on
plot([0:0.1:0.5],[0:0.1:0.5],'k--')
axis square
xlabel('Vis')
ylabel('Aud')
title('avg std: 1st half of trial')
subplot(1,2,2)
scatter(mean(std_vis((minTrialLengthFrames/2)+pre_event_frames:minTrialLengthFrames+pre_event_frames,:),1),mean(std_aud((minTrialLengthFrames/2)+pre_event_frames:minTrialLengthFrames+pre_event_frames,:),1),'k');
hold on
scatter(mean(mean(std_vis((minTrialLengthFrames/2)+pre_event_frames:minTrialLengthFrames+pre_event_frames,:),1)),mean(mean(std_aud((minTrialLengthFrames/2)+pre_event_frames:minTrialLengthFrames+pre_event_frames,:),1)),'ro','markerfacecolor','r');
plot([0:0.1:0.5],[0:0.1:0.5],'k--')
axis square
xlabel('Vis')
ylabel('Aud')
title('avg std: 2nd half of trial')

stdAvsVScat_CD = figure;
subplot(1,2,1)
hStdV = cdfplot(mean(std_vis(pre_event_frames:(minTrialLengthFrames/2)+pre_event_frames,:),1));
hStdV.Color = 'g'
hold on
hStdA = cdfplot(mean(std_aud(pre_event_frames:(minTrialLengthFrames/2)+pre_event_frames,:),1));
hStdA.Color = 'k'
[h p] = kstest2(mean(std_vis(pre_event_frames:(minTrialLengthFrames/2)+pre_event_frames,:),1),mean(std_aud(pre_event_frames:(minTrialLengthFrames/2)+pre_event_frames,:),1));
title(['avg std: 1st half; p=' num2str(p)])
axis square
subplot(1,2,2)
hStdV = cdfplot(mean(std_vis(pre_event_frames:minTrialLengthFrames+pre_event_frames,:),1));
hStdV.Color = 'g';
hold on
hStdA = cdfplot(mean(std_aud(pre_event_frames:minTrialLengthFrames+pre_event_frames,:),1));
hStdA.Color = 'k';
[h p] = kstest2(mean(std_vis(pre_event_frames:minTrialLengthFrames+pre_event_frames,:),1),mean(std_aud(pre_event_frames:minTrialLengthFrames+pre_event_frames,:),1));
title(['avg std: all; p=' num2str(p)])
axis square
print([fnout 'press_align_stds' datasetStr '.pdf'], '-dpdf')

figure(stdAvsVScat);
print([fnout 'press_align_stds_CD' datasetStr '.pdf'], '-dpdf')

% %% plot avg trace for random subset of cells
% 
% for imouse = 1:size(mouse,2)
%     for iexp = 1:size(mouse(imouse).expt,2)
%         figure;
%         suptitle({mouse(imouse).expt(iexp).date, mouse(imouse).expt(iexp).mouse_name})
%         cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
%         if length(cell_ind) > 9
%             cell_ind = sort(randsample(cell_ind,9));
%         end
%         for icell = 1:length(cell_ind)
%             subplot(3,3,icell)
%             plot(tt,mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind(icell),:),3),'g')
%             hold on
%             plot(tt,mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind(icell),:),3),'k')
%             hold on
%             vline(baseStimFrames,':k')
%             title(['cell# ' num2str(cell_ind(icell))])
%         end
%     end
% end

%% plot scatter of integral for each experiment, across experiments;
antiRespScatters

resp_vis = [];
resp_aud = [];
resp_vis_early = [];
resp_aud_early = [];
resp_vis_late = [];
resp_aud_late = [];
intAvsVFig = figure;
suptitle(titleStr)
i = 1;

for imouse = 1:size(mouse,2)
    resp_vis_mouse{imouse} = [];
    resp_aud_mouse{imouse} = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        figure(intAvsVFig);
        subplot(n,n2,i)
        i = i+1;
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        iV = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_event_frames:end,cell_ind,:)),3));
        iA = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(pre_event_frames:end,cell_ind,:)),3));
        % add STE color axis here
        scatter(iV,iA,50,'k.');
        hold on
        errorbarxy(mean(iV),mean(iA),std(iV)/length(iV),std(iA)/length(iA),{'ro','r','r'});
%         xlim([floor(min([iV iA]))-(floor(min([iV iA]))/2) ceil(max([iV iA]))+(ceil(max([iV iA]))/2)]);
%         ylim([floor(min([iV iA]))-(floor(min([iV iA]))/2) ceil(max([iV iA]))+(ceil(max([iV iA]))/2)]);
            xlim([-5 15])
            ylim([-5 15])
%         
        hold on
        plot([-10:0.1:20],[-10:0.1:20],'k--')
        axis square
        xlabel('Vis Tr Resp')
        ylabel('Aud Tr Resp')
        
        resp_vis = cat(2, resp_vis, iV);
        resp_aud = cat(2, resp_aud, iA);
        
        resp_vis_mouse{imouse} = cat(2,resp_vis_mouse{imouse}, iV);
        resp_aud_mouse{imouse} = cat(2,resp_aud_mouse{imouse}, iA);

        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
    end
end

subplot(n,n2,i)
scatter(resp_vis,resp_aud,50,'k.')
hold on
errorbarxy(mean(resp_vis),mean(resp_aud),std(resp_vis)/length(resp_vis),std(resp_aud)/length(resp_aud),{'ro','r','r'})
hold on
for imouse = 1:size(mouse,2)
   errorbarxy(mean(resp_vis_mouse{imouse}),mean(resp_aud_mouse{imouse}),std(resp_vis_mouse{imouse})/length(resp_vis_mouse{imouse}),std(resp_aud_mouse{imouse})/length(resp_vis_mouse{imouse}), {[av(mouse_ind(imouse)).col_str 'o'],av(mouse_ind(imouse)).col_str,av(mouse_ind(imouse)).col_str});
   hold on
   mouse_mean(imouse) = scatter(mean(resp_vis_mouse{imouse}),mean(resp_aud_mouse{imouse}),av(mouse_ind(imouse)).col_str,'o', 'filled');
   %    mouseLegend{imouse} = str{imouse};
   hold on
end
legend(mouse_mean,str,'Location','southeast')
% xlim([floor(min([resp_vis resp_aud]))-(floor(min([resp_vis resp_aud]))/2) ceil(max([resp_vis resp_aud]))+(ceil(max([resp_vis resp_aud]))/2)]);
% ylim([floor(min([resp_vis resp_aud]))-(floor(min([resp_vis resp_aud]))/2) ceil(max([resp_vis resp_aud]))+(ceil(max([resp_vis resp_aud]))/2)]);
xlim([-5 15])
ylim([-5 15])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')

resp_vis_all_int = resp_vis;
resp_aud_all_int = resp_aud;

title(['All cells; n = ' num2str(size(resp_vis,2))])

figure(intAvsVFig);
print([fnout 'press_align_int' datasetStr '.pdf'], '-dpdf')

%% individual int scatter of all cells with mean
figure;
subplot(3,2,1)
scatter(resp_vis,resp_aud,100,'k.');
hold on
H = errorbarxy(mean(resp_vis),mean(resp_aud),std(resp_vis)/length(resp_vis),std(resp_aud)/sqrt(length(resp_aud)),{'r.','r','r'});
H.hMain.MarkerSize = 20;
hold on
xlim([-5 15])
ylim([-5 15])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
print([fnout 'press_align_int_all' datasetStr '.pdf'], '-dpdf')

%% plot scatter of integral for each experiment, across experiments; ORI SELECT LABELLED 
int_vis_ori = cell([1,4]);
int_aud_ori = cell([1,4]);
resp_vis_ori = cell([1,4]);
resp_aud_ori = cell([1,4]);
std_vis_ori = cell([1,4]);
std_aud_ori = cell([1,4]);
intAvsVFig_ori = figure;
suptitle(titleStr)
i = 1;
oriColMat = ['k';'g';'b';'r'];

for imouse = 1:size(mouse,2)
%     resp_vis_mouse{imouse} = [];
%     resp_aud_mouse{imouse} = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        figure(intAvsVFig_ori);
        subplot(n,n2,i)
        i = i+1;
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        for iori = 1:4
            ori_ind = intersect(mouse(imouse).expt(iexp).cells(iori+1).ind,cell_ind);
            iV = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_event_frames:end,ori_ind,:)),3));
            iA = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(pre_event_frames:end,ori_ind,:)),3));
            tcV = mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,ori_ind,:),3);
            tcA = mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,ori_ind,:),3);
            % add STE color axis here
            scatter(iV,iA,50,[oriColMat(iori) '.']);
            hold on
%             errorbarxy(mean(iV),mean(iA),std(iV)/length(iV),std(iA)/length(iA),{[oriColMat(iori) 'o'],oriColMat(iori),oriColMat(iori)});
            oriMean(iori) = scatter(mean(iV),mean(iA),oriColMat(iori),'o','filled');
            if length([iV iA]) > 0;
%             xlim([floor(min([iV iA]))-(floor(min([iV iA]))/2) ceil(max([iV iA]))+(ceil(max([iV iA]))/2)]);
%             ylim([floor(min([iV iA]))-(floor(min([iV iA]))/2) ceil(max([iV iA]))+(ceil(max([iV iA]))/2)]);
            xlim([-5 15])
            ylim([-5 15])
            end
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('Vis Tr Resp')
            ylabel('Aud Tr Resp')
            
            
            int_vis_ori{iori} = cat(2, int_vis_ori{iori}, iV);
            int_aud_ori{iori} = cat(2, int_aud_ori{iori}, iA);
            resp_vis_ori{iori} = cat(2, resp_vis_ori{iori}, tcV);
            resp_aud_ori{iori} = cat(2, resp_aud_ori{iori}, tcA);
        end

        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
    end
end

std_vis_ori = cellfun(@(x) std(x,[],2),resp_vis_ori,'unif',false);
std_aud_ori = cellfun(@(x) std(x,[],2),resp_aud_ori,'unif',false);

subplot(n,n2,i)
% hold on
% errorbarxy(mean(resp_vis_ori),mean(resp_aud_ori),std(resp_vis_ori)/length(resp_vis_ori),std(resp_aud_ori)/length(resp_aud_ori),{'ro','r','r'})
% hold on
oriLegend = {'0','45','90','135'};
for iori = 1:4
    scatter(int_vis_ori{iori},int_aud_ori{iori},50,oriColMat(iori),'.')
    hold on
    errorbarxy(mean(int_vis_ori{iori}),mean(int_aud_ori{iori}),std(int_vis_ori{iori})/length(int_vis_ori{iori}),std(int_aud_ori{iori})/length(int_aud_ori{iori}), {[oriColMat(iori) 'o'],oriColMat(iori),oriColMat(iori)});
    hold on
    ori_mean(iori) = scatter(mean(int_vis_ori{iori}),mean(int_aud_ori{iori}),oriColMat(iori),'o', 'filled');
    hold on
end
legend(ori_mean,oriLegend,'Location','southeast')
% xlim([floor(min([resp_vis_ori resp_aud_ori]))-(floor(min([resp_vis_ori resp_aud_ori]))/2) ceil(max([resp_vis_ori resp_aud_ori]))+(ceil(max([resp_vis_ori resp_aud_ori]))/2)]);
% ylim([floor(min([resp_vis_ori resp_aud_ori]))-(floor(min([resp_vis_ori resp_aud_ori]))/2) ceil(max([resp_vis_ori resp_aud_ori]))+(ceil(max([resp_vis_ori resp_aud_ori]))/2)]);
xlim([-5 15])
ylim([-5 15])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')

title(['All cells; n = ' num2str(size(int_vis_ori,2))])

%cdf plot of tuned cell respones

respAvsVori = figure;
oriLegend = {'0','45','90','135'};
subplot(3,2,2)
hV = [];
hA = [];
for iori = 1:4
    if ~isempty(int_vis_ori{iori})
    hV= cdfplot(int_vis_ori{iori});
    hV.Color = oriColMat(iori);
    hV.LineStyle = '-';
    hold on
    hA = cdfplot(int_aud_ori{iori});
    hA.Color = oriColMat(iori)
    hA.LineStyle = '--';
    hold on
    end
end
xlabel('dF/F')
ylabel('fraction cells')
title('integral response - ori selective cells')
% legend(oriLegend)
print([fnout 'press_align_intOri_cdf' datasetStr '.pdf'], '-dpdf')

figure;
suptitle({titleStr; 'std timecourse'})
for iori = 1:4
    subplot(3,2,iori)
    plot(ttMs,std_vis_ori{iori},'g')
    hold on
    plot(ttMs,std_aud_ori{iori},'k')
    hold on
    xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
    ylim([-0.01 0.05]);
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    title([oriStr(iori) ' slctv cells; n = ' num2str(size(std_vis_ori{iori},2))])
    xlabel('ms')
    ylabel('std - dF/F')
end

print([fnout 'press_align_stdTCbyOri' datasetStr '.pdf'], '-dpdf')


figure;
for iori = 1:4
    scatter(int_vis_ori{iori},int_aud_ori{iori},50,oriColMat(iori),'.')
    hold on
    errorbarxy(mean(int_vis_ori{iori}),mean(int_aud_ori{iori}),std(int_vis_ori{iori})/length(int_vis_ori{iori}),std(int_aud_ori{iori})/length(int_aud_ori{iori}), {[oriColMat(iori) 'o'],oriColMat(iori),oriColMat(iori)});
    hold on
    ori_mean(iori) = scatter(mean(int_vis_ori{iori}),mean(int_aud_ori{iori}),oriColMat(iori),'o', 'filled');
    hold on
end
legend(ori_mean,oriLegend,'Location','southeast')
% xlim([floor(min([resp_vis_ori resp_aud_ori]))-(floor(min([resp_vis_ori resp_aud_ori]))/2) ceil(max([resp_vis_ori resp_aud_ori]))+(ceil(max([resp_vis_ori resp_aud_ori]))/2)]);
% ylim([floor(min([resp_vis_ori resp_aud_ori]))-(floor(min([resp_vis_ori resp_aud_ori]))/2) ceil(max([resp_vis_ori resp_aud_ori]))+(ceil(max([resp_vis_ori resp_aud_ori]))/2)]);
xlim([-5 10])
ylim([-5 10])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
title(['All cells; n = ' num2str(size(int_vis_ori,1))])
print([fnout 'press_align_intOriall' datasetStr '.pdf'], '-dpdf')

figure(intAvsVFig_ori);
print([fnout 'press_align_intOri' datasetStr '.pdf'], '-dpdf')


figure;
for iori = 1:4
    subplot(3,2,iori)
    if sum(~isnan(int_vis_ori{iori})) > 4
        hV = cdfplot(int_vis_ori{iori});
        hV.Color = 'g';
        hold on
        hA = cdfplot(int_aud_ori{iori});
        hA.Color = 'k';
        [h,p] = kstest2(int_vis_ori{iori},int_aud_ori{iori});
    else
        p = NaN
    end
    title([mouse(imouse).expt(iexp).cells(iori+1).name ' cells; n = ' num2str(sum(cell2mat(cellfun(@length,int_vis_ori,'unif',false)))) '; p = ' num2str(chop(p,2))])
    xlabel('dF/F')
end
resp_vis_all = cell2mat(int_vis_ori);
resp_aud_all = cell2mat(int_aud_ori);
subplot(3,2,5)
if sum(~isnan(resp_vis_all)) > 4
    hV = cdfplot(resp_vis_all);
    hV.Color = 'g';
    hold on
    hA = cdfplot(resp_aud_all);
    hA.Color = 'k';
    [h,p] = kstest2(resp_vis_all,resp_aud_all);
    else
        p = NaN;
end
    title(['all tuned cells; n= ' num2str(length(resp_vis_all)) '; p=' num2str(p)])
    xlabel('dF/F')
    subplot(3,2,6)
if sum(~isnan(resp_vis_all_int)) > 4
    hV = cdfplot(resp_vis_all_int);
    hV.Color = 'g';
    hold on
    hA = cdfplot(resp_aud_all_int);
    hA.Color = 'k';
    [h,p] = kstest2(resp_vis_all_int,resp_aud_all_int);
    else
        p = NaN;
end
    title(['all resp cells; n= ' num2str(length(resp_vis_all_int)) '; p=' num2str(p)])
    xlabel('dF/F')
suptitle([titleStr '- Vis:Green, Aud:Black'])


print([fnout 'press_align_intOri_CD' datasetStr '.pdf'], '-dpdf')


% %% plot histogram of number of excited, inhibited, and non-responsive cells by resp integral
% 
% histEvsIvsNRcells_V = figure;
% histEvsIvsNRcells_A = figure;
% nbins = 100;
% i = 1;
% resp_vis = [];
% resp_aud = [];
% exc_ind_all = [];
% inh_ind_all = [];
% % resp_FA = [];
% % resp_CR = [];
% 
% for imouse = 1:size(mouse,2)
%     for iexp = 1:size(mouse(imouse).expt,2)
%         
%         cell_ind = 1:mouse(imouse).expt(iexp).info.nCells;
%         exc_ind = mouse(imouse).expt(iexp).cells(6).ind;
%         inh_ind = mouse(imouse).expt(iexp).cells(7).ind;
%         
%         iV = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_event_frames:end,cell_ind,:)),3));
%         iA = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(pre_event_frames:end,cell_ind,:)),3));
%         
%         bin_edges = linspace(min([iV iA]),max([iV iA]),nbins);
%         
%         [barVexc barVexcInd] = histc(iV(exc_ind),bin_edges);
%         [barVinh barVinhInd] = histc(iV(inh_ind),bin_edges);
%         [barVnr barVnrInd] = histc(iV(setdiff(cell_ind, [exc_ind inh_ind])),bin_edges);
%         
%         [barAexc barAexcInd] = histc(iA(exc_ind),bin_edges);
%         [barAinh barAinhInd] = histc(iA(inh_ind),bin_edges);
%         [barAnr barAnrInd] = histc(iA(setdiff(cell_ind, [exc_ind inh_ind])),bin_edges);
%         
%         figure(histEvsIvsNRcells_V);
%         subplot(n,n2,i)
%         hE_V = bar(bin_edges,barVexc,'r');
%         hE_V_ch = get(hE_V,'child');
%         set(hE_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
%         hold on
%         hI_V = bar(bin_edges,barVinh,'b');
%         hI_V_ch = get(hI_V,'child');
%         set(hI_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
%         hold on
%         hNR_V = bar(bin_edges,barVnr,'k');
%         hNR_V_ch = get(hNR_V,'child');
%         set(hNR_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
%         hold on
%         vline(0,'k')
%         title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
%         xlabel('resp int')
%         ylabel('# cells')
%         
%         figure(histEvsIvsNRcells_A);
%         subplot(n,n2,i)
%         hE_A = bar(bin_edges,barAexc,'r');
%         hE_A_ch = get(hE_A,'child');
%         set(hE_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
%         hold on
%         hI_A = bar(bin_edges,barVinh,'b');
%         hI_A_ch = get(hI_A,'child');
%         set(hI_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
%         hold on
%         hNR_A = bar(bin_edges,barVnr,'k');
%         hNR_A_ch = get(hNR_A,'child');
%         set(hNR_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
%         hold on
%         vline(0,'k')
%         title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
%         xlabel('resp int')
%         ylabel('# cells')
%         
%         resp_vis = cat(2, resp_vis, iV);
%         resp_aud = cat(2, resp_aud, iA);
%         
%         if i > 1
%             exc_ind_all = cat(2,exc_ind_all,exc_ind+runCells);
%             inh_ind_all = cat(2,inh_ind_all,inh_ind+runCells);
%             runCells = length(cell_ind)+runCells;
%         else
%             exc_ind_all = exc_ind;
%             inh_ind_all = inh_ind;
%             runCells = length(cell_ind);
%         end
%         
%         i = i+1;
%     end
% end
% 
% bin_edges = linspace(min([resp_vis resp_aud]),max([resp_vis resp_aud]),nbins);
%         
% [barVexc barVexcInd] = histc(resp_vis(exc_ind_all),bin_edges);
% [barVinh barVinhInd] = histc(resp_vis(inh_ind_all),bin_edges);
% [barVnr barVnrInd] = histc(resp_vis(setdiff(1:runCells, [exc_ind_all inh_ind_all])),bin_edges);
% 
% [barAexc barAexcInd] = histc(resp_aud(exc_ind_all),bin_edges);
% [barAinh barAinhInd] = histc(resp_aud(inh_ind_all),bin_edges);
% [barAnr barAnrInd] = histc(resp_aud(setdiff(1:runCells, [exc_ind_all inh_ind_all])),bin_edges);
% 
% figure(histEvsIvsNRcells_V);
% subplot(n,n2,i)
% hE_V = bar(bin_edges,barVexc,'r');
% hE_V_ch = get(hE_V,'child');
% set(hE_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
% hold on
% hI_V = bar(bin_edges,barVinh,'b');
% hI_V_ch = get(hI_V,'child');
% set(hI_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
% hold on
% hNR_V = bar(bin_edges,barVnr,'k');
% hNR_V_ch = get(hNR_V,'child');
% set(hNR_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
% %         set(gca,'XTick',bin_edges)
% hold on
% vline(0,'k')
% title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
% xlabel('resp int')
% ylabel('# cells')
% 
% figure(histEvsIvsNRcells_A);
% subplot(n,n2,i)
% hE_A = bar(bin_edges,barAexc,'r');
% hE_A_ch = get(hE_A,'child');
% set(hE_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
% hold on
% hI_A = bar(bin_edges,barVinh,'b');
% hI_A_ch = get(hI_A,'child');
% set(hI_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
% hold on
% hNR_A = bar(bin_edges,barVnr,'k');
% hNR_A_ch = get(hNR_A,'child');
% set(hNR_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
% %         set(gca,'XTick',bin_edges)
% hold on
% vline(0,'k')
% title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
% xlabel('resp int')
% ylabel('# cells')
% 
% figure(histEvsIvsNRcells_V);
% print([fnout 'press_align_Vint_hist' datasetStr '.pdf'], '-dpdf')
% figure(histEvsIvsNRcells_A);
% print([fnout 'press_align_Aint_hist' datasetStr '.pdf'], '-dpdf')
end
