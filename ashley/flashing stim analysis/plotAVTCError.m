function plotAVTCError(visTCAllTrials,audTCAllTrials,responseLim,...
    nBaselineFrames,timeLim,frameRateHz)

timestamp = (-nBaselineFrames+1:size(visTCAllTrials,1)-nBaselineFrames)/frameRateHz;
timeAxisLabel = [-0.25:0.25:timeLim(end)];

nTrialsEaCondition = [{size(visTCAllTrials,3)};{size(audTCAllTrials,3)}];
legendText = cellfun(@(x,y) sprintf([x ' (%s)'],num2str(y)),{'vis';'aud'},...
    nTrialsEaCondition,'unif',0);

tc = squeeze(mean(mean(visTCAllTrials,2),3));
err = squeeze(ste(mean(visTCAllTrials,2),3));
h = shadedErrorBar(timestamp,tc,err,'k');
lines4legend(1) = h.mainLine;
hold on
tc = squeeze(mean(mean(audTCAllTrials,2),3));
err = squeeze(ste(mean(audTCAllTrials,2),3));
h = shadedErrorBar(timestamp,tc,err,'c');
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