function drawCalciumHistograms(tt,ratio,delta,deconv,stimulus,scale,binsiz,os)
%DRAWCALCIUMHISTOGRAMS
%DRAWCALCIUMHISTOGRAMS(TT,RATIO,DELTA)
%DRAWCALCIUMHISTOGRAMS(TT,RATIO,DELTA,DECONV)
%DRAWCALCIUMHISTOGRAMS(TT,RATIO,DELTA,DECONV,STIMULUS)
%DRAWCALCIUMHISTOGRAMS(TT,RATIO,DELTA,DECONV,STIMULUS,SCALE)
%DRAWCALCIUMHISTOGRAMS(TT,RATIO,DELTA,DECONV,STIMULUS,SCALE,BINSIZ)
%DRAWCALCIUMHISTOGRAMS(TT,RATIO,DELTA,DECONV,STIMULUS,SCALE,BINSIZ,OS)

if nargin < 6;    scale = 1;end
if nargin < 7;  binsiz = 2;end;
if nargin < 8; os = -delta/2 ;end;

hold off;
bg = 0.9;

[nSamples,nTrials]=size(ratio);

yl = [-delta/2 (nTrials)*delta];

% draw stimulus 
if nargin > 4
    h1 = area(tt,stimulus*yl(2));
    hold on;
    h2 = area(tt,stimulus*yl(1));
    set(h1,'facecolor',bg*[1 1 1],'edgecolor',bg*[1 1 1]);
    set(h2,'facecolor',bg*[1 1 1],'edgecolor',bg*[1 1 1]);
    set(gca,'ycolor','k');
end

% draw spike responses
if nargin > 3  & ~isempty(deconv)
	hh = offsetbar(tt(1:binsiz:end),scale*rasterbin(deconv,binsiz),delta,os,'edgecolor','r','facecolor','r');
    hold on;
end

% plot time courses
if ~isempty(ratio);
    tcOffsetPlot(tt,ratio,delta,'k');
    hold on;

    % vertical scale bar
    plot([0 0],[0 25],'k','linewidth',5);
end

% plot firing events as red ticks
% if nargin > 3  & ~isempty(deconv)
%     events = double(deconv>0);
%     events(find(~events))=nan;
% 
%     xx = [tt;tt];xx = xx(:);
%     yy = permute(events,[3 1 2]);
%     yy = [-delta/6*yy;-delta*yy/3];
%     yy = reshape(yy,[length(tt)*2,nTrials]);
%     tcOffsetPlot(xx,yy,delta,'r','linewidth',1);
%     nshift = nshift + 1;
% end

ylim(yl);

xlim([min(tt) max(tt)]);

xlabel('Time (s)');
ylabel('\DeltaF/F (%)');    

hold off;

return;
