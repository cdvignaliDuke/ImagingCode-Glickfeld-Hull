function longTimeCourseEaDelayFig = plotLongTimeCourseStartAlign(neuronTimeCourse,...
    timeCourseStartAlignMean,timeCourseStartAlignErr,timestamp,delayColors)
nDelay = length(timeCourseStartAlignMean);
longTimeCourseEaDelayFig = figure;
subplot 211
longTimeCourseLegend = [];
for idelay = 1:nDelay
    y = timeCourseStartAlignMean{idelay};
    yerr = timeCourseStartAlignErr{idelay};
    h = shadedErrorBar(timestamp,y,yerr,[],1);
    h.mainLine.Color = delayColors(idelay,:);
    h.mainLine.LineWidth = 2;  
    hold on
    longTimeCourseLegend(idelay) = h.mainLine;
end
figXAxis([],'time(s)',[timestamp(1) timestamp(end)])
figYAxis([],'dF/F',[])
ax = gca;
ax.TickDir = 'out';
ax.Box = 'off';
legend(longTimeCourseLegend,{neuronTimeCourse(1:nDelay).trType},...
    'location','northeastoutside')
subplot 212
longTimeCourseLegend = [];
for idelay = 1:nDelay
    y = timeCourseStartAlignMean{idelay};
    h = plot(timestamp,y);
    h.Color = delayColors(idelay,:);
    h.LineWidth = 2;  
    hold on
    longTimeCourseLegend(idelay) = h;
end
figXAxis([],'time(s)',[timestamp(1) timestamp(end)])
figYAxis([],'dF/F',[])
ax = gca;
ax.TickDir = 'out';
ax.Box = 'off';
legend(longTimeCourseLegend,{neuronTimeCourse(1:nDelay).trType},...
    'location','northeastoutside')
end