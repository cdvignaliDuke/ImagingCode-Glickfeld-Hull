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

subplot(3,2,1)
scatter(earlyResp_vis,earlyResp_aud,50,'k.')
hold on
errorbarxy(mean(earlyResp_vis,2),mean(earlyResp_aud,2), (std(earlyResp_vis,[],2)/sqrt(size(earlyResp_vis,2))) , (std(earlyResp_aud,[],2)/sqrt(size(earlyResp_aud,2))) ,{'ro','r','r'});
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

subplot(3,2,2)
scatter(lateResp_vis,lateResp_aud,50,'k.')
hold on
errorbarxy(mean(lateResp_vis,2),mean(lateResp_aud,2), (std(lateResp_vis,[],2)/sqrt(size(lateResp_vis,2))) , (std(lateResp_aud,[],2)/sqrt(size(lateResp_aud,2))) ,{'ro','r','r'});
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

subplot(3,2,3)
eVis = cdfplot(earlyResp_vis);
eVis.Color = 'g';
hold on
eAud = cdfplot(earlyResp_aud);
eAud.Color = 'k';
xlim([-1 2])
axis square
xlabel('int - dF/F')
ylabel('cmlv frctn')
legend({'visual';'auditory'})
title(['early - kstest p=' num2str(pkEarly)])

subplot(3,2,4)
lVis = cdfplot(lateResp_vis);
lVis.Color = 'g';
hold on
lAud = cdfplot(lateResp_aud);
lAud.Color = 'k';
xlim([-1 2])
axis square
xlabel('int - dF/F')
ylabel('cmlv frctn')
legend({'visual';'auditory'})
title(['late - kstest p=' num2str(pkLate)])

%% plot timecourse
tcEV = mean(resp_vis_cmlvCyc{earlyCycMax},2);
erEV = std(resp_vis_cmlvCyc{earlyCycMax},[],2)/sqrt(size(resp_vis_cmlvCyc{earlyCycMax},2));
tcEA = mean(resp_aud_cmlvCyc{earlyCycMax},2);
erEA = std(resp_aud_cmlvCyc{earlyCycMax},[],2)/sqrt(size(resp_aud_cmlvCyc{earlyCycMax},2));
tcLV = mean(resp_vis_cmlvCyc{lateCycMax},2);
erLV = std(resp_vis_cmlvCyc{lateCycMax},[],2)/sqrt(size(resp_vis_cmlvCyc{lateCycMax},2));
tcLA = mean(resp_aud_cmlvCyc{lateCycMax},2);
erLA = std(resp_aud_cmlvCyc{lateCycMax},[],2)/sqrt(size(resp_aud_cmlvCyc{lateCycMax},2));

figure(antiRespQuant)

subplot(3,2,5)
shadedErrorBar(ttCycMs(20:length(tcEV)),tcEV(20:end),erEV(20:end),'g');
hold on
shadedErrorBar(ttCycMs(20:length(tcEA)),tcEA(20:end),erEA(20:end),'k');
hold on
xlim([-10 cycTime*earlyCycMax]/(cycTime/cycTimeMs))
if cellsInd == 2
ylim([-0.01 0.05])
else
ylim([-0.01 0.03])
end
vline([ttCycMs(early_win(1)) ttCycMs(early_win(end))], '--r')
vline(baseStimFrames/(cycTime/cycTimeMs),':k');
title([num2str(earlyCycMax) ' cycles, n = ' num2str(size(resp_vis_cmlvCyc{earlyCycMax},2))]);

subplot(3,2,6)
shadedErrorBar(ttCycMs(20:length(tcLV)),tcLV(20:end),erLV(20:end),'g');
hold on
shadedErrorBar(ttCycMs(20:length(tcLA)),tcLA(20:end),erLA(20:end),'k');
hold on
xlim([-10 cycTime*lateCycMax]/(cycTime/cycTimeMs))
if cellsInd == 2
ylim([-0.01 0.05])
else
ylim([-0.01 0.03])
end
vline([ttCycMs(late_win(1)) ttCycMs(late_win(end))], '--r')
vline([0 cycTime:cycTime:cycTime*lateCycMax]/(cycTime/cycTimeMs),':k');
title([num2str(lateCycMax) ' cycles, n = ' num2str(size(resp_vis_cmlvCyc{lateCycMax},2))]);

%%
figure(antiRespQuant)
print([fnout 'press_align_antiEarlyLateQuant' datasetStr '.pdf'], '-dpdf')
