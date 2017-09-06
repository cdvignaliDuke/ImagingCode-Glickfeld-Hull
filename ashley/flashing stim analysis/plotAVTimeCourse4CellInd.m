function plotAVTimeCourse4CellInd(tcVis,tcAud,cellIndex,nBaselineFr,...
    frameRateHz,responseAxisLimit)

timestampS = ((1:size(tcVis,1))-nBaselineFr)./frameRateHz;
timeTick = chop([-timestampS(end-1), 0 timestampS(end-1)],2);
timeLim = [-0.5 timestampS(end)];

nVisTrials = size(tcVis,3);
nAudTrials = size(tcAud,3);

meanTCVis = mean(mean(tcVis(:,cellIndex,:),3),2);
errTCVis = ste(mean(tcVis(:,cellIndex,:),3),2);
meanTCAud = mean(mean(tcAud(:,cellIndex,:),3),2);
errTCAud = ste(mean(tcAud(:,cellIndex,:),3),2);

legendText_2cyc = cellfun(@(x,y) sprintf([x ' (%s)'],y),{'vis';'aud'},...
    {num2str(nVisTrials);num2str(nAudTrials)},'unif',0);
h = shadedErrorBar(timestampS,meanTCVis,errTCVis,'k');
legendPlots(1) = h.mainLine;
hold on
h = shadedErrorBar(timestampS,meanTCAud,errTCAud,'c');
legendPlots(2) = h.mainLine;
plot(timestampS,meanTCVis,'k-');
figXAxis([],'time(s)',timeLim,timeTick,timeTick)
figYAxis([],'dF/F',responseAxisLimit)
legend(legendPlots,legendText_2cyc,'location','northwest')

end