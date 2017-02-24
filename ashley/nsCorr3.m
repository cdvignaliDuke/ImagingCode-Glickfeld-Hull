ialign = 1;
vis = 1;
aud = 2;
hits = 1;
minTrL = 2800;
early_win = [500 1000];
late_win = [2300 2800];

nsCorrSubFig = figure;
nsCorrGrpFig = figure;
expt_colors = jet(length(expt)+5);
expt_colors = expt_colors(3:length(expt)+3,:);
ms_colors = brewermap(size(mouse,2),'Dark2');
i = 1;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)  
    clear neuronCorr
    cell_ind = mouse(imouse).expt(iexp).cells(14).ind;
    expt_name = [mouse(imouse).expt(iexp).mouse_name '-' mouse(imouse).expt(iexp).date];
%how many cycles is a trial at least 3s?
cycTimeMs = mouse(imouse).expt(iexp).info.cyc_time_ms;
cycTime = mouse(imouse).expt(iexp).info.cyc_time;
cyc3s = ceil(minTrL/cycTimeMs);
fr1s = (1000*cycTime)/cycTimeMs;
early_win_frames = floor(early_win*(cycTime/cycTimeMs))+pre_event_frames;
late_win_frames = floor(late_win*(cycTime/cycTimeMs))+pre_event_frames;

if length(mouse(imouse).expt(iexp).align(ialign).av(vis).outcome(hits).cmlvCycResp) >= cyc3s
% get data for 3s trial
RV = mouse(imouse).expt(iexp).align(ialign).av(vis).outcome(hits).cmlvCycResp{cyc3s};
RA = mouse(imouse).expt(iexp).align(ialign).av(aud).outcome(hits).cmlvCycResp{cyc3s};

% chunk out early and late analysis windows
RV_early = RV(early_win_frames(1):early_win_frames(2),:,:);
RA_early = RA(early_win_frames(1):early_win_frames(2),:,:);
RV_late = RV(late_win_frames(1):late_win_frames(2),:,:);
RA_late = RA(late_win_frames(1):late_win_frames(2),:,:);
%line up all trial data
RV_early_alignTrials = reshape(permute(RV_early,[1 3 2]),[size(RV_early,1)*size(RV_early,3) size(RV_early,2)]);
RA_early_alignTrials = reshape(permute(RA_early,[1 3 2]),[size(RA_early,1)*size(RA_early,3) size(RA_early,2)]);
RV_late_alignTrials = reshape(permute(RV_late,[1 3 2]),[size(RV_late,1)*size(RV_late,3) size(RV_late,2)]);
RA_late_alignTrials = reshape(permute(RA_late,[1 3 2]),[size(RA_late,1)*size(RA_late,3) size(RA_late,2)]);
%subtract time-course mean, then correlate
corrV_early = corrcoef(bsxfun(@minus,RV_early_alignTrials,mean(RV_early_alignTrials,2)));
corrA_early = corrcoef(bsxfun(@minus,RA_early_alignTrials,mean(RA_early_alignTrials,2)));
corrV_late = corrcoef(bsxfun(@minus,RV_late_alignTrials,mean(RV_late_alignTrials,2)));
corrA_late = corrcoef(bsxfun(@minus,RA_late_alignTrials,mean(RA_late_alignTrials,2)));

figure;
colormap(brewermap([],'*RdBu'))
suptitle(expt_name)
subplot(2,3,1)
imagesc(corrV_early)
title('vis, early')
subplot(2,3,2)
imagesc(corrA_early)
title('aud, early')
subplot(2,3,3)
imagesc(corrV_early-corrA_early)
title('vis minus aud')
subplot(2,3,4)
imagesc(corrV_late)
title('vis, late')
subplot(2,3,5)
imagesc(corrA_late)
title('aud, late')
subplot(2,3,6)
imagesc(corrV_late-corrA_late)
title('vis, late')
for iplot = 1:6
    subplot(2,3,iplot)
    axis square
    colorbar
    caxis([-1 1])
end

% get mean ns corr for different cell groups;

neuronCorr(1).name = 'all';
neuronCorr(2).name = 'task responsive';
neuronCorr(3).name = 'non-responsive';
neuronCorr(4).name = '0 selective';
neuronCorr(5).name = '45 selective';
neuronCorr(6).name = '90 selective';
neuronCorr(7).name = '135 selective';
neuronCorr(1).ind = cell_ind;
neuronCorr(2).ind = setdiff(1:size(RV,2),cell_ind);
neuronCorr(3).ind = intersect(cell_ind,mouse(imouse).expt(iexp).cells(2).ind);
neuronCorr(4).ind = intersect(cell_ind,mouse(imouse).expt(iexp).cells(3).ind);
neuronCorr(5).ind = intersect(cell_ind,mouse(imouse).expt(iexp).cells(4).ind);
neuronCorr(6).ind = intersect(cell_ind,mouse(imouse).expt(iexp).cells(5).ind);
neuronCorr(7).ind = intersect(cell_ind,mouse(imouse).expt(iexp).cells(6).ind);

for iInd = 1:7
    neuronCorr(iInd).corrV_early = corrV_early(neuronCorr(iInd).ind,neuronCorr(iInd).ind);
    neuronCorr(iInd).corrA_early = corrA_early(neuronCorr(iInd).ind,neuronCorr(iInd).ind);
    neuronCorr(iInd).corrV_late = corrV_late(neuronCorr(iInd).ind,neuronCorr(iInd).ind);
    neuronCorr(iInd).corrA_late = corrA_late(neuronCorr(iInd).ind,neuronCorr(iInd).ind);
end

%plot mean ns corr for each group
figure(nsCorrGrpFig)
for iInd = 1:7;
    
    tempdata = neuronCorr(iInd).corrV_early;
    tempdata = tril(tempdata,-1);
    tempdata = tempdata(:);
    xplt = mean(tempdata(tempdata~=0));
    xste = std(tempdata(tempdata~=0))/sqrt(length(tempdata(tempdata~=0)));
    tempdata = neuronCorr(iInd).corrA_early;
    tempdata = tril(tempdata,-1);
    tempdata = tempdata(:);
    yplt = mean(tempdata(tempdata~=0));
    yste = std(tempdata(tempdata~=0))/sqrt(length(tempdata(tempdata~=0)));

    plotInd = (iInd*2)-1;
    subplot(7,2,plotInd);
    errorbarxy(xplt,yplt,xste,yste,{'ko-', 'k', 'k'});
    hold on
    h = plot(xplt,yplt,'ko');
    h.MarkerEdgeColor = 'non';
    h.MarkerFaceColor = ms_colors(imouse,:);

    tempdata = neuronCorr(iInd).corrV_late;
    tempdata = tril(tempdata,-1);
    tempdata = tempdata(:);
    xplt = mean(tempdata(tempdata~=0));
    xste = std(tempdata(tempdata~=0))/sqrt(length(tempdata(tempdata~=0)));
    tempdata = neuronCorr(iInd).corrA_late;
    tempdata = tril(tempdata,-1);
    tempdata = tempdata(:);
    yplt = mean(tempdata(tempdata~=0));
    yste = std(tempdata(tempdata~=0))/sqrt(length(tempdata(tempdata~=0)));
    
    plotInd = iInd*2;
    subplot(7,2,plotInd);
    errorbarxy(xplt,yplt,xste,yste,{'ko-', 'k', 'k'});
    hold on
    h = plot(xplt,yplt,'ko');
    h.MarkerEdgeColor = 'non';
    h.MarkerFaceColor = ms_colors(imouse,:);

end

% %plot example cell time-course and 1s moving mean time-course 
% iCell = cell_ind(3);
% exCellTC_V = mean(RV(:,iCell,:),3);
% exCellTC_V_smthd = movmean(exCellTC_V,fr1s);
% exCellTC_A = mean(RA(:,iCell,:),3);
% exCellTC_A_smthd = movmean(exCellTC_A,fr1s);
% figure;
% plot(ttMs(1:size(RV,1)),exCellTC_V,'g')
% hold on
% plot(ttMs(1:size(RA,1)),exCellTC_A,'k')
% hold on
% plot(ttMs(1:length(exCellTC_V_smthd)),exCellTC_V_smthd,'g--')
% hold on
% plot(ttMs(1:length(exCellTC_A_smthd)),exCellTC_A_smthd,'k--')
% hold on
% vline(baseStimMs,'k:')
% xlabel('ms')
% ylabel('dF/F')
% title(['Cell ' num2str(iCell)])
% 
% % get sliding window mean for all cells 
% RV_resp = mean(RV,3);
% RV_sw = movmean(RV,fr1s);
% RA_resp = mean(RA,3);
% RA_sw = movmean(RA,fr1s);
% 
% 
% RV_corr = zeros(size(RV,2),size(RV,2),size(RV,1));
% RA_corr = zeros(size(RA,2),size(RA,2),size(RA,1));
% for icorr = 1:size(RA_sw,1)
%     tempdata = squeeze(RV_sw(icorr,:,:));
%     RV_corr(:,:,icorr) = corrcoef(tempdata');
%     tempdata = squeeze(RA_sw(icorr,:,:));
%     RA_corr(:,:,icorr) = corrcoef(tempdata');
% end
% 
% % plot time-course of ns corr for example cell 
% exCellCorr_V = squeeze(mean(RV_corr(iCell,:,:),2));
% exCellCorr_A = squeeze(mean(RA_corr(iCell,:,:),2));
% figure;
% h = plot(ttMs(1:size(RV,1)),exCellCorr_V,'k');
% hold on
% h = plot(ttMs(1:size(RA,1)),exCellCorr_A);
% h.Color = [0.5 0.5 0.5];
% hold on
% vline(baseStimMs,'k:')
% xlabel('ms')
% ylabel('ns corr')
% 
% %plot overlay of all task-responsive cells (visual trials)
% trcCorr_V = squeeze(mean(RV_corr(cell_ind,:,:),2));
% trcCorr_A = squeeze(mean(RA_corr(cell_ind,:,:),2));
% figure;
% title('task responsive cells')
% subplot(1,2,1)
% plot(ttMs(1:size(RV,1)),trcCorr_V);
% hold on
% vline(baseStimMs,'k:')
% xlabel('ms')
% ylabel('ns corr')
% subplot(1,2,2)
% plot(ttMs(1:size(RA,1)),trcCorr_A);
% hold on
% vline(baseStimMs,'k:')
% xlabel('ms')
% ylabel('ns corr')
% 
% % plot mean ns corr of task-responsive cells vs non-responsive cells
% trcCorr_V_mean = mean(trcCorr_V,1);
% nrCorr_V_mean = squeeze(mean(mean(RV_corr(setdiff(1:size(RV,2),cell_ind),:,:),1),2));
% trcCorr_A_mean = mean(trcCorr_A,1);
% nrCorr_A_mean = squeeze(mean(mean(RA_corr(setdiff(1:size(RA,2),cell_ind),:,:),1),2));
% figure;
% suptitle(expt_name)
% subplot(1,2,1)
% plot(ttMs(1:size(RV,1)),trcCorr_V_mean,'k');
% hold on
% plot(ttMs(1:size(RV,1)),nrCorr_V_mean,'k--');
% hold on
% vline(baseStimMs,'k:')
% xlabel('ms')
% ylabel('ns corr')
% legend({'task-repsonsive';'non-responsive'});
% title('visual trials')
% subplot(1,2,2)
% plot(ttMs(1:size(RA,1)),trcCorr_A_mean,'k');
% hold on
% plot(ttMs(1:size(RA,1)),nrCorr_A_mean,'k--');
% hold on
% vline(baseStimMs,'k:')
% xlabel('ms')
% ylabel('ns corr')
% legend({'task-repsonsive';'non-responsive'});
% title('auditory trials')
% 
% % mean ns corr vis vs aud, task responsive cells;
% figure;
% suptitle(expt_name)
% subplot(1,2,1)
% h = plot(ttMs(1:size(RV,1)),trcCorr_V_mean,'k');
% hold on
% h = plot(ttMs(1:size(RA,1)),trcCorr_A_mean,'k');
% h.Color = [0.5 0.5 0.5];
% hold on
% legend({'visual';'auditory'});
% vline(baseStimMs,'k:')
% xlabel('ms')
% ylabel('ns corr')
% title('task-responsive')
% 
% %subtract time-courses, plot overlay of all experiments
% 
% figure(nsCorrSubFig);
% hold on
% H(i) = plot(ttMs(1:size(RV,1)),trcCorr_V_mean-trcCorr_A_mean,'color',expt_colors(i,:));
% hold on
% vline(baseStimMs,'k:')
% xlabel('ms')
% ylabel('ns corr')
% title({'task-responsive cells'; 'vis minus aud'})
% H_leg{i} = expt_name;
% i = i+1;
end
    end
end

% figure(nsCorrSubFig);
% legend(H,H_leg)

figure(nsCorrGrpFig)
for iInd = 1:7
    plotInd = (iInd*2)-1;
    subplot(7,2,plotInd)
    hold on
    plot(-1:0.1:1,-1:0.1:1,'k--');
    title(['early-' neuronCorr(iInd).name])
    xlim([-1 1])
    ylim([-1 1])
    xlabel('visual')
    ylabel('auditory')
    axis square
    
    plotInd = iInd*2;
    subplot(7,2,plotInd)
    hold on
    plot(-1:0.1:1,-1:0.1:1,'k--');
    title(['late'])
    xlim([-1 1])
    ylim([-1 1])
    xlabel('visual')
    ylabel('auditory')
    axis square
end


print([fnout 'nsCorrMeanScatter' datasetStr '.pdf'], '-dpdf','-fillpage')