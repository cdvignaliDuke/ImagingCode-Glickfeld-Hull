function drawCalciumResponses(tt,ratio,delta,deconv,stimulus,thresh)
%DRAWCALCIUMRESPONSES
%DRAWCALCIUMRESPONSES(TT,RATIO,DELTA)
%DRAWCALCIUMRESPONSES(TT,RATIO,DELTA,DECONV)
%DRAWCALCIUMRESPONSES(TT,RATIO,DELTA,DECONV,STIMULUS)
%DRAWCALCIUMRESPONSES(TT,RATIO,DELTA,DECONV,STIMULUS,THRESHOLD)

if nargin < 6
    thresh = 2.0;
end


hold off;
bg = 0.9;

[nSamples,nTrials]=size(ratio);

yl = [-delta/2 (nTrials)*delta];

% draw stimulus 
if nargin > 4 & ~isempty(stimulus)
    h1 = area(tt,stimulus*yl(2));
    hold on;
    h2 = area(tt,stimulus*yl(1));
    set(h1,'facecolor',bg*[1 1 1],'edgecolor',bg*[1 1 1]);
    set(h2,'facecolor',bg*[1 1 1],'edgecolor',bg*[1 1 1]);
    set(gca,'ycolor','k');
end

% plot firing events as red ticks
if nargin > 3  & ~isempty(deconv)
    events = double(deconv>thresh);
    events(find(~events))=nan;

    xx = [tt;tt];xx = xx(:);
    yy = permute(events,[3 1 2]);
    yy = [-delta/6*yy;-delta*yy/3];
    yy = reshape(yy,[length(tt)*2,nTrials]);
    tcOffsetPlot(xx,yy,delta,'r','linewidth',1);
    hold on;
end

% plot time courses
if ~isempty(ratio);
    tcOffsetPlot(tt,ratio,delta,'k');
    hold on;

    % vertical scale bar
    plot([0 0],[0 25],'k','linewidth',5);
end

ylim(yl);

xlim([min(tt) max(tt)]);

xlabel('Time (s)');
ylabel('\DeltaF/F (%)');    

hold off;

return;
