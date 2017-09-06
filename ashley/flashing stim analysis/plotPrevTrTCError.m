function plotPrevTrTCError(prevVisTCAllTrials,prevAudTCAllTrials,trialType,...
    responseLim,nBaselineFrames,timeLim,frameRateHz)

timestamp = (-nBaselineFrames+1:size(prevVisTCAllTrials,1)-nBaselineFrames)/frameRateHz;
timeAxisLabel = [-0.25:0.25:timeLim(end)];

nTrialsEaCondition = [{size(prevVisTCAllTrials,3)};{size(prevAudTCAllTrials,3)}];
legendText = cellfun(@(x,y) sprintf([x ' (%s)'],num2str(y)),...
    {'prev vis';'prev aud'},nTrialsEaCondition,'unif',0);
if strcmp(trialType,'visual')
    prevVisColor = [0 0 0];
    prevAudColor = [0.75 0.75 0.75];
else
    prevVisColor = [0 0.75 0.75];
    prevAudColor = [0 1 1];    
end

tc = squeeze(mean(mean(prevVisTCAllTrials,2),3));
err = squeeze(ste(mean(prevVisTCAllTrials,2),3));
h = shadedErrorBar(timestamp,tc,err,'k');
h.mainLine.Color = prevVisColor;
lines4legend(1) = h.mainLine;
hold on
tc = squeeze(mean(mean(prevAudTCAllTrials,2),3));
err = squeeze(ste(mean(prevAudTCAllTrials,2),3));
h = shadedErrorBar(timestamp,tc,err,'k');
h.mainLine.Color = prevAudColor;
lines4legend(2) = h.mainLine;
figXAxis([],'time (s)',timeLim,timeAxisLabel,timeAxisLabel)
figYAxis([],'dF/F',responseLim)
vline(0,'k--')
fig = gca;
fig.Box = 'off';
fig.TickDir = 'out';
leg = legend(lines4legend,legendText,'location','northwest');
title(leg,'trials (n)')
end