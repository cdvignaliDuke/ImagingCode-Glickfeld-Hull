earlyCycMin = 3;
lateCycMin = 7;
earlyCycMax = 4;
lateCycMax = 8;
ttCycMs = (-pre_event_frames:((lateCycMax(end)*cycTime))-1)/(cycTime/cycTimeMs);

nAnalysisFrames_early = 0:cycTime*2-1;
nAnalysisFrames_late = 0:(cycTime*2-1)-(trans_win(1)-pre_event_frames-1);
early_win = nAnalysisFrames_early+trans_win(1)+cycTime;
late_win = nAnalysisFrames_late+(trans_win(1)+((lateCycMin-1)*cycTime));

earlyResp_vis = trapz(resp_vis_cmlvCyc{earlyCycMax}(early_win,:));
earlyResp_aud = trapz(resp_aud_cmlvCyc{earlyCycMax}(early_win,:));
[hEarly,pEarly] = ttest(earlyResp_vis,earlyResp_aud,'alpha',0.05/size(earlyResp_vis,2),'tail','both');
[hkEarly,pkEarly] = kstest2(earlyResp_vis,earlyResp_aud);

lateResp_vis = trapz(resp_vis_cmlvCyc{lateCycMax}(late_win,:));
lateResp_aud = trapz(resp_aud_cmlvCyc{lateCycMax}(late_win,:));
[hLate,pLate] = ttest(lateResp_vis,lateResp_aud,'alpha',0.05/size(lateResp_vis,2),'tail','both');
[hkLate,pkLate] = kstest2(lateResp_vis,lateResp_aud);

%% scatter plot aud vs vis responses
antiRespQuant = figure;
suptitle(titleStr)

subplot(4,2,1)
scatter(earlyResp_vis,earlyResp_aud,50,'k.')
hold on
errorbarxy(nanmean(earlyResp_vis,2),nanmean(earlyResp_aud,2), (nanstd(earlyResp_vis,[],2)/sqrt(size(earlyResp_vis,2))) , (nanstd(earlyResp_aud,[],2)/sqrt(size(earlyResp_aud,2))) ,{'ro','r','r'});
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
vline(0,'b')
hline(0,'b')
xlim([-3 8])
ylim([-3 8])
axis square
xlabel('visual')
ylabel('auditory')
title({'early - int of analysis win'; ['pairwise p=' num2str(pEarly)]})

subplot(4,2,2)
scatter(lateResp_vis,lateResp_aud,50,'k.')
hold on
errorbarxy(nanmean(lateResp_vis,2),nanmean(lateResp_aud,2), (nanstd(lateResp_vis,[],2)/sqrt(size(lateResp_vis,2))) , (nanstd(lateResp_aud,[],2)/sqrt(size(lateResp_aud,2))) ,{'ro','r','r'});
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
vline(0,'b')
hline(0,'b')
xlim([-3 8])
ylim([-3 8])
axis square
xlabel('visual')
ylabel('auditory')
title({'late - int of analysis win';['pairwise p=' num2str(pLate)]})

%% cdf plots early and late integral
figure(antiRespQuant)

subplot(4,2,3)
eVis = cdfplot(earlyResp_vis);
eVis.Color = 'g';
hold on
eAud = cdfplot(earlyResp_aud);
eAud.Color = 'k';
xlim([-1 2])
axis square
xlabel('int - dF/F')
ylabel('cmlv frctn')
legend({'visual';'auditory'},'location','southeastoutside')
title(['early - kstest p=' num2str(pkEarly)])

subplot(4,2,4)
lVis = cdfplot(lateResp_vis);
lVis.Color = 'g';
hold on
lAud = cdfplot(lateResp_aud);
lAud.Color = 'k';
xlim([-1 2])
axis square
xlabel('int - dF/F')
ylabel('cmlv frctn')
legend({'visual';'auditory'},'location','southeastoutside')
title(['late - kstest p=' num2str(pkLate)])

%% plot timecourse
tcEV = nanmean(resp_vis_cmlvCyc{earlyCycMax},2);
erEV = nanstd(resp_vis_cmlvCyc{earlyCycMax},[],2)/sqrt(size(resp_vis_cmlvCyc{earlyCycMax},2));
tcEA = nanmean(resp_aud_cmlvCyc{earlyCycMax},2);
erEA = nanstd(resp_aud_cmlvCyc{earlyCycMax},[],2)/sqrt(size(resp_aud_cmlvCyc{earlyCycMax},2));
tcLV = nanmean(resp_vis_cmlvCyc{lateCycMax},2);
erLV = nanstd(resp_vis_cmlvCyc{lateCycMax},[],2)/sqrt(size(resp_vis_cmlvCyc{lateCycMax},2));
tcLA = nanmean(resp_aud_cmlvCyc{lateCycMax},2);
erLA = nanstd(resp_aud_cmlvCyc{lateCycMax},[],2)/sqrt(size(resp_aud_cmlvCyc{lateCycMax},2));


tcEsub = nanmean(resp_vis_cmlvCyc{earlyCycMax}-resp_aud_cmlvCyc{earlyCycMax},2);
erEsub = nanstd(resp_vis_cmlvCyc{earlyCycMax}-resp_aud_cmlvCyc{earlyCycMax},[],2)/sqrt(size(resp_vis_cmlvCyc{earlyCycMax},2));
tcLsub = nanmean(resp_vis_cmlvCyc{lateCycMax}-resp_aud_cmlvCyc{lateCycMax},2);
erLsub = nanstd(resp_vis_cmlvCyc{lateCycMax}-resp_aud_cmlvCyc{lateCycMax},[],2)/sqrt(size(resp_vis_cmlvCyc{lateCycMax},2));

figure(antiRespQuant)

subplot(4,2,5)
shadedErrorBar(ttCycMs(20:length(tcEV)),tcEV(20:end),erEV(20:end),'g');
hold on
shadedErrorBar(ttCycMs(20:length(tcEA)),tcEA(20:end),erEA(20:end),'k');
hold on
xlim([-10 cycTime*earlyCycMax]/(cycTime/cycTimeMs))
if cellsInd == 2 | strcmp(datasetStr,'_audControl') | strcmp(datasetStr,'_V1axonsAL') | strcmp(datasetStr,'_V1axonsPM')
ylim([-0.01 0.06])
else
ylim([-0.01 0.03])
end
vline([ttCycMs(early_win(1)) ttCycMs(early_win(end))], '--r')
vline(baseStimFrames/(cycTime/cycTimeMs),':k');
title([num2str(earlyCycMax) ' cycles, n = ' num2str(size(resp_vis_cmlvCyc{earlyCycMax},2))]);


subplot(4,2,7)
shadedErrorBar(ttCycMs(20:length(tcEsub)),tcEsub(20:end),erEsub(20:end),'k');
hold on
xlim([-10 cycTime*earlyCycMax]/(cycTime/cycTimeMs))
if cellsInd == 2 | strcmp(datasetStr,'_audControl') | strcmp(datasetStr,'_V1axonsAL') | strcmp(datasetStr,'_V1axonsPM')
ylim([-0.01 0.02])
else
ylim([-0.01 0.01])
end
vline([ttCycMs(early_win(1)) ttCycMs(early_win(end))], '--r')
vline(baseStimFrames/(cycTime/cycTimeMs),':k');
title([num2str(earlyCycMax) ' cycles, n = ' num2str(size(resp_vis_cmlvCyc{earlyCycMax},2)) '- V-A sub']);

subplot(4,2,6)
shadedErrorBar(ttCycMs(20:length(tcLV)),tcLV(20:end),erLV(20:end),'g');
hold on
shadedErrorBar(ttCycMs(20:length(tcLA)),tcLA(20:end),erLA(20:end),'k');
hold on
xlim([-10 cycTime*lateCycMax]/(cycTime/cycTimeMs))
if cellsInd == 2 | strcmp(datasetStr,'_audControl') | strcmp(datasetStr,'_V1axonsAL') | strcmp(datasetStr,'_V1axonsPM')
ylim([-0.01 0.06])
else
ylim([-0.01 0.03])
end
vline([ttCycMs(late_win(1)) ttCycMs(late_win(end))], '--r')
vline([0 cycTime:cycTime:cycTime*lateCycMax]/(cycTime/cycTimeMs),':k');
title([num2str(lateCycMax) ' cycles, n = ' num2str(size(resp_vis_cmlvCyc{lateCycMax},2))]);

subplot(4,2,8)
shadedErrorBar(ttCycMs(20:length(tcLsub)),tcLsub(20:end),erLsub(20:end),'k');
hold on
xlim([-10 cycTime*lateCycMax]/(cycTime/cycTimeMs))
if cellsInd == 2 | strcmp(datasetStr,'_audControl') | strcmp(datasetStr,'_V1axonsAL')
ylim([-0.01 0.02])
else
ylim([-0.01 0.01])
end
vline([ttCycMs(late_win(1)) ttCycMs(late_win(end))], '--r')
vline([0 cycTime:cycTime:cycTime*lateCycMax]/(cycTime/cycTimeMs),':k');
title([num2str(lateCycMax) ' cycles, n = ' num2str(size(resp_vis_cmlvCyc{lateCycMax},2)) '- V-A sub']);
%%
% modIndex_vis_aud
%%
figure(antiRespQuant)
print([fnout 'press_align_antiEarlyLateQuant' datasetStr '.pdf'], '-dpdf','-fillpage')
