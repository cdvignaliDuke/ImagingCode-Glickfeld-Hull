%% horizontal violin plots
% requires data in the form of:
% x (all cells numbered by area 1-4)
% yPS (all cells prefSize)
% ySI (all cells SI)
% alternative can use with data_cell (4 cells of pS/SI for each area)

% convert yPS to data_cell
for i=1:4
    data_cell{i} = [];
    data_cell{i} = yPS(x==i);
end

y_mean = [mean(data_cell{1}) mean(data_cell{2}) mean(data_cell{3}) mean(data_cell{4})];
y_std = [std(data_cell{1}) std(data_cell{2}) std(data_cell{3}) std(data_cell{4})];

figure(3);clf;set(gcf,'Color','w');
% 1. average
subplot(2,2,1)
errorbar(1:4,y_mean,y_std,'ok');
set(gca,'box','off','TickDir','out')
set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
ylabel('Pref Size (deg)')
xlim([0.5 4.5])
ylim([0 60])
%5 horizontal violin
for i = 1:length(areas)
    [fi xi] = ksdensity(data_cell{i});
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(2,2,2);
colors = get(ax,'ColorOrder');
for i=1:length(areas)
    hold on
    h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
    p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
    h5(i).FaceColor = colors(i,:);
end
axis([0 100 0.5 length(areas)+0.5]);
legend off
ax.YTick = [1:4];
ax.YTickLabel = fliplr(areas);
set(gca,'box','off','TickDir','out')
ylabel('Area')
xlabel('Pref size (deg)')

% for SI
for i=1:4
    data_cell{i} = [];
    data_cell{i} = ySI(x==i);
end

y_mean = [mean(data_cell{1}) mean(data_cell{2}) mean(data_cell{3}) mean(data_cell{4})];
y_std = [std(data_cell{1}) std(data_cell{2}) std(data_cell{3}) std(data_cell{4})];

% 1. average
subplot(2,2,3)
errorbar(1:4,y_mean,y_std,'ok');
set(gca,'box','off','TickDir','out')
set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
ylabel('SI')
xlim([0.5 4.5])
ylim([0 1.1])
%5 sideways violin
for i = 1:length(areas)
    [fi xi] = ksdensity(data_cell{i},'BoundaryCorrection','reflection','Support',[0-eps 1+eps]); %,'BoundaryCorrection','reflection'
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(2,2,4);
colors = get(ax,'ColorOrder');
for i=1:length(areas)
    hold on
    h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
    p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
    h5(i).FaceColor = colors(i,:);
end
axis([0 1 0.5 length(areas)+0.5]);
legend off
ax.YTick = [1:4];
ax.YTickLabel = fliplr(areas);
set(gca,'box','off','TickDir','out')
ylabel('Area')
xlabel('SI')

filename = 'N:\home\kevin\ppts\_paper figs\plots\PS_SI_testplots.pdf';
set(gcf,'PaperSize',[12 6])
%print(filename,'-dpdf','-fillpage')