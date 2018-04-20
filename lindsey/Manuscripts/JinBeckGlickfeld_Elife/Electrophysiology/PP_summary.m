% [Gdata, Gtext, Graw] = xlsread('Y:\Analysis\electrode\ducumentation\PP.xlsx',3);
% 
% layers = [];
% response = [];
% FR = [];
% 
% for i_exp = 1:size(Gdata,1)
%     
%    load(cell2mat(Gtext(i_exp+1,2)))
%    idx = PP.cellanova; 
%    layers = [layers PP.layers(idx)];
%    response = [response;PP.response.normalize.mean(idx,:)];
%    
%    FR = [FR;PP.response.mean(idx,:)];
%    
%    
% 
%   
% 
%    
% end 


%% plot the overlayed recovery across areas for all layers combined or only layer2/3
% not: in V1, layer 3 means layer 2/3. 4means layer 4, 5 means layer 5
% in HVAs: 1 means layer1, 2 means initial sinks, presumably layer2/3 maybe some layer4?,
% 3 and 4 means anything deep layers....
areas={'V1','LM','PM','AL','SC','LGN'};
colors = [0 0 0; lines(5)];
for i=1:length(areas)
    Gdata=[];
    Gtext=[];
    Graw = [];
    [Gdata, Gtext, Graw] = xlsread('Y:\Analysis\electrode\ducumentation\PP.xlsx',i);
    
    layers{i,1} = [];
    response{i,1} = [];
    FR{i,1} = [];
    pt_ratio{i,1} = [];
    pt_ms{i,1} = [];
    pt_ratio_all{i,1} = [];
    pt_ms_all{i,1} = [];
    waveform_all{i,1} = [];
    waveform{i,1} = [];
    layers_all{i,1} = [];
    traces{i,1} = [];
    for i_exp = 1:size(Gdata,1)
        
        load(cell2mat(Gtext(i_exp+1,2)))
        load(cell2mat(Gtext(i_exp+1,5)))
        
        idx = PP.cellanova; 
        %idx = PP.anovastrict;
        
        layers{i,1} = [layers{i,1} PP.layers(idx)];
        response{i,1} = [response{i,1};PP.response.normalize.mean(idx,:)];       
        
        FR{i,1} = [FR{i,1};PP.response.mean(idx,:)];
       
        traces{i,1} = cat(1,traces{i,1},PP.Traces.raw_mean(idx,:,:));
        pt_ratio{i,1} = [ pt_ratio{i,1};spikes.wave_stats.pt_ratio(idx)];
        pt_ms{i,1} = [pt_ms{i,1}; spikes.wave_stats.pt_ms(idx)];
        waveform{i,1} = [waveform{i,1},spikes.template_big(:,idx)];
        
        pt_ratio_all{i,1} = [ pt_ratio_all{i,1};spikes.wave_stats.pt_ratio];
        pt_ms_all{i,1} = [pt_ms_all{i,1}; spikes.wave_stats.pt_ms];
        waveform_all{i,1} = [waveform_all{i,1},spikes.template_big];
        layers_all{i,1} = [layers{i,1} PP.layers];
        
        
        
    end
    
    
end
%% plot the adaptation ratio over responses 
binwidth = PP.binwidth;
edges = PP.edges;
FR_cat = [];
response_cat = [];

for i = 1:4
    FR_cat = cat(1,FR_cat,FR{i,1});
    response_cat = cat(1,response_cat,response{i,1});
    
end 

% caculate adaptation index as (control-adapted)./(control+adapted)
% first change all the negative values into 0
FR_cat_zero = FR_cat;
FR_cat_zero(FR_cat<0)=0;
Adapt_index = [];
Adapt_index = (FR_cat_zero-repmat(FR_cat_zero(:,6),1,6) )./(repmat(FR_cat_zero(:,6),1,6) + FR_cat_zero);

figure
% subplot(1,2,1)
% scatter(FR_cat(:,6)./binwidth,response_cat(:,1),'MarkerEdgeColor',[0 0 0])
% axis square
subplot(1,2,1)
scatter(FR_cat(:,6)./binwidth,Adapt_index(:,1),'MarkerEdgeColor',[0 0 0])
hold on
hline(0,'k:')
ylim([-1 1])
axis square
 [R,P] = corrcoef(FR_cat(:,6)./binwidth,Adapt_index(:,1)); 
 text(60, 0.8,['r = ' num2str(R(2,1)) ' p = ' num2str(P(2,1))])
 ylabel('Adaptation index')
 xlabel('Firing Rate to the first response')
 set(gca, 'TickDir','out')
% subplot(1,2,2)
% histogram(Adapt_index(:,1),20)
% xlim([-1 1])
% axis square
%% plot the adaptation index across areas
hist_edges = -1:0.1:1; 
oder = [1 4 3 2 5];
figure
for i = 1:5
    temp = [];
    temp = FR{i,1};
    temp_zero = temp;
    temp_zero(temp<0)=0;
    adpt{i,1} = [];
    adpt{i,1} = (temp_zero-repmat(temp_zero(:,6),1,6) )./(repmat(temp_zero(:,6),1,6) + temp_zero);
    
    
%     subplot(1,2,1)
%     scatter(FR{i,1}(:,6)./binwidth, adpt{i,1}(:,1),'MarkerEdgeColor',colors(i,:))
%     hold on
%     hline(0,'k:')
%     ylim([-1 1])
%     axis square
    % [R,P] = corrcoef(FR_cat(:,6)./binwidth,Adapt_index(:,1));
    
    subplot(5,1,oder(i))
    histogram(adpt{i,1}(:,5),hist_edges,'FaceColor',colors(i,:),'EdgeColor','none')
    hold on
    scatter(mean(adpt{i,1}(:,5)),30,'MarkerEdgeColor',colors(i,:),'Marker','v')
    
    xlim([-1.1 1.1])
    ylim([0 35])
    text(0.8,30,areas(i))
    set(gca,'TickDir','out')
    
end


%% plot V1 different layers
i=1;
figure
subplot(3,1,1)
idx = layers{i,1}==3;
histogram(adpt{i,1}(idx,1),hist_edges,'FaceColor',colors(i,:),'EdgeColor','none')
hold on
scatter(mean(adpt{i,1}(idx,1)),18,'MarkerEdgeColor',colors(i,:),'Marker','v')

xlim([-1.1 1.1])
ylim([0 20])
text(0.8,18,'L2/3')
set(gca,'TickDir','out')

subplot(3,1,2)
idx = layers{i,1}==4;
histogram(adpt{i,1}(idx,1),hist_edges,'FaceColor',[0 0.6 0.6],'EdgeColor','none')
hold on
scatter(mean(adpt{i,1}(idx,1)),18,'MarkerEdgeColor',[0 0.6 0.6],'Marker','v')

xlim([-1.1 1.1])
ylim([0 20])
text(0.8,18,'L4')
set(gca,'TickDir','out')

subplot(3,1,3)
idx = layers{i,1}>5;
histogram(adpt{i,1}(idx,1),hist_edges,'FaceColor',[1 0.2 0.4],'EdgeColor','none')
hold on
scatter(mean(adpt{i,1}(idx,1)),18,'MarkerEdgeColor',[1 0.2 0.4],'Marker','v')

xlim([-1.1 1.1])
ylim([0 20])
text(0.8,18,'Deep')
set(gca,'TickDir','out')


%% Plot the in the initial sink vs below the sink
for i = 2:4
figure
subplot(2,1,1)
idx = layers{i,1}==2;
histogram(adpt{i,1}(idx,1),hist_edges,'FaceColor',colors(i,:),'EdgeColor','none')
hold on
scatter(mean(adpt{i,1}(idx,1)),18,'MarkerEdgeColor',colors(i,:),'Marker','v')

xlim([-1.1 1.1])
ylim([0 20])
text(0.8,18,'super')
set(gca,'TickDir','out')

subplot(2,1,2)
idx = layers{i,1}>3;
histogram(adpt{i,1}(idx,1),hist_edges,'FaceColor',colors(i,:),'EdgeColor','none')
hold on
scatter(mean(adpt{i,1}(idx,1)),18,'MarkerEdgeColor',colors(i,:),'Marker','v')

xlim([-1.1 1.1])
ylim([0 20])
text(0.8,18,'deep')
set(gca,'TickDir','out')


supertitle(areas{i})
    
end 






%% plot the average response traces for the first pulse
binwidth = PP.binwidth;
edges = PP.edges;
idx = edges<=0;
pulse = 6;
figure

i=5;

temp = [];
temp = traces{i,1}-repmat(squeeze(mean(traces{i,1}(:,pulse,idx),3)),1,6,28);
y_values=squeeze(mean(temp(:,pulse,:),1))./binwidth; 
errorbars = squeeze(std(temp(:,pulse,:),[],1))./(binwidth*sqrt(length(layers{i,1}))); 

shadedErrorBar_ch(edges,y_values,errorbars,{'color',[1 0.8 0.4],'linewidth',1},0);

 hold on

i=1;
temp = [];
temp = traces{i,1}-repmat(squeeze(mean(traces{i,1}(:,pulse,idx),3)),1,6,28);
y_values=squeeze(mean(temp(layers{i,1}==4,pulse,:),1))./binwidth; 
errorbars = squeeze(std(temp(layers{i,1}==4,pulse,:),[],1))./(binwidth*sqrt(sum(layers{i,1}==4))); 

shadedErrorBar_ch(edges,y_values,errorbars,{'color',[0 0.6 0.6],'linewidth',1},0);
hold on

y_values=squeeze(mean(temp(layers{i,1}~=4,pulse,:),1))./binwidth; 
errorbars = squeeze(std(temp(layers{i,1}~=4,pulse,:),[],1))./(binwidth*sqrt(sum(layers{i,1}~=4))); 

shadedErrorBar_ch(edges,y_values,errorbars,{'color',[0 0 0],'linewidth',1},0);
xlim([-0.1 0.3])

hold on

% 
% i=6;
% 
% temp = [];
% temp = traces{i,1}-repmat(squeeze(mean(traces{i,1}(:,pulse,idx),3)),1,6,28);
% y_values=squeeze(mean(temp(:,pulse,:),1))./binwidth; 
% errorbars = squeeze(std(temp(:,pulse,:),[],1))./(binwidth*sqrt(length(layers{i,1}))); 
% 
% shadedErrorBar_ch(edges,y_values,errorbars,{'color',[0.2 0.6 0],'linewidth',1},0);

% cat all the HVAs 
temp = [];
for i =2:4
temp = cat(1, temp, traces{i,1}-repmat(squeeze(mean(traces{i,1}(:,pulse,idx),3)),1,6,28));
end 
y_values=squeeze(mean(temp(:,pulse,:),1))./binwidth; 
errorbars = squeeze(std(temp(:,pulse,:),[],1))./(binwidth*sqrt(size(temp,1))); 

shadedErrorBar_ch(edges,y_values,errorbars,{'color',[1 0.2 0.4],'linewidth',1},0);
text(0.2,16,'V1-L4','Color',[0 0.6 0.6])
text(0.2,15,'V1-L2/3/5','Color',[0 0 0])
text(0.2,14,'HVAs','Color',[1 0.2 0.4])
text(0.2,13,'SC','Color',[1 0.8 0.4])
ylim([-3 18])
set(gca,'XTick',-0.1:0.1:0.3,'YTick',-3:3:18,'TickDir','out')
axis square
xlabel('Time(s)')
ylabel('Firing rate(Hz)')
%% plot the pp recovery
figure
temp = [];
i=1;
temp = response{i,1}(layers{i,1}==4,:);

% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[0 0.6 0.6]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[0 0.6 0.6]);
text(3000,0.6,[areas{i},'-L4',': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',[0 0.6 0.6])

temp = [];
i=1;
temp = response{i,1}(layers{i,1}~=4,:);
% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[0 0 0]);
text(3000,0.5,[areas{i},'-L2/3/5/6',': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',[0 0 0])

% collapsed all HVAs



temp = [];
for i=2:4
temp = cat(1,temp,response{i,1});
end 
% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[1 0.2 0.4]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[1 0.2 0.4]);
text(3000,0.4,['HVAs',': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',[1 0.2 0.4])



temp = [];
i=5;
temp = response{i,1};
% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[1 0.8 0.4]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[1 0.8 0.4]);
text(3000,0.3,[areas{i},': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',[1 0.8 0.4])
ylim([0 1.2])
set(gca,'Xtick',[0 2000 4000 8000])
set(gca,'Xticklabel',{'0';'2';'4';'baseline'},'TickDir','out')
ylabel('Normalized FR')
xlabel('Inter-stimulus interval (s)')
axis square

%% plot the pp of superficial vs deep
figure
temp = [];
i=1;
temp = response{i,1}(layers{i,1}==4,:);

% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[0 0.6 0.6]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[0 0.6 0.6]);
text(3000,0.6,[areas{i},'L4',': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000)],'color',[0 0.6 0.6])

temp = [];

temp = response{i,1}(layers{i,1}==3,:);
% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[0 0 0]);
text(3000,0.5,[areas{i},'L2/3',': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ],'color',[0 0 0])


temp = [];
temp = response{i,1}(layers{i,1}>4,:);
% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[1 0.2 0.4]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[1 0.2 0.4]);
text(3000,0.4,[areas{i},'deep layers',': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ],'color',[1 0.2 0.4])



ylim([0 1.2])
set(gca,'Xtick',[0 2000 4000 8000])
set(gca,'Xticklabel',{'0';'2';'4';'baseline'},'TickDir','out')
ylabel('Normalized FR')
xlabel('Inter-stimulus interval (s)')
axis square


%% plot the pp of superficial vs deep for HVAs
response_HVA = [];
layers_HVA  = [];
for i = 2:4
    
    response_HVA = cat(1,response_HVA,response{i,1});
    layers_HVA = cat(2,layers_HVA,layers{i,1});
end

figure
temp = [];
temp = response_HVA(layers_HVA<=2,:);

% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[0 0.6 0.6]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[0 0.6 0.6]);
text(3000,0.6,['HVA superficial',': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',[0 0.6 0.6])

temp = [];
temp = response_HVA(layers_HVA>2,:);
% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[0 0 0]);
text(3000,0.5,['HVA deep',': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',[0 0 0])




ylim([0 1.2])
set(gca,'Xtick',[0 2000 4000 8000])
set(gca,'Xticklabel',{'0';'2';'4';'baseline'},'TickDir','out')
ylabel('Normalized FR')
xlabel('Inter-stimulus interval (s)')
axis square

%% plot the HVA difference
for i = 1:4
    temp = [];
    temp = response{i,1};
% get the fit
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

fitsingle_all(i,:) =  [tao,A,sse,R_square];
h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',colors(i,:));
h.LineStyle = 'none';
h.Marker = 'o';

hold on

DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,colors(i,:));
text(3000,0.5-(i-1)*0.1,[areas{i},': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',colors(i,:))

end 
ylim([0 1.1])
set(gca,'Xtick',[0 2000 4000 8000])
set(gca,'Xticklabel',{'0';'2';'4';'baseline'},'TickDir','out')
ylabel('Normalized FR')
xlabel('Inter-stimulus interval (s)')
axis square
%% plot the scatter for all sorted neurons
for i =4%:length(areas)
    figure
    subplot(2,2,1)
    scatter(pt_ratio{i,1}, pt_ms{i,1})
    hold on
    data = [pt_ratio{i,1}, pt_ms{i,1}];
    options = statset('Display','final');
    
    obj = fitgmdist(data,2);
    h = ezcontour(@(x,y)pdf(obj,[x y]),[0.1 0.5],[0 1.2]);
    xlabel('peak:trough height ratio')
    ylabel('trough to peak (ms)')
    
    idx = cluster(obj,data);
    cluster1 = data(idx ==1,:);
    cluster2 = data(idx ==2,:);
    
    
    subplot(2,2,2)
    h1 = scatter(cluster1(:,1),cluster1(:,2),10,'ro');
    hold on
    h2 = scatter(cluster2(:,1),cluster2(:,2),10,'bo');
     legend([h1 h2],'Broad spiking','Narrow spiking','Location','NW')
     xlim([0.1 0.5])
     ylim([0 1.2])
    xlabel('peak:trough height ratio')
    ylabel('trough to peak (ms)')
    
    subplot(2,2,3)
    wavedege = 0:1./30000:60./30000;
    plot( wavedege*1000,waveform{i,1}(:,idx ==1)./(-repmat(min(waveform{i,1}(:,idx ==1),[],1),61,1)),'r')
    hold on
    plot( wavedege*1000,waveform{i,1}(:,idx ==2)./(-repmat(min(waveform{i,1}(:,idx ==2),[],1),61,1)),'b')
    xlabel('ms')
    ylabel('normalized amplitude')
    % plot the adaptation broad spiking vs narrow spiking neurons for all
    % layaers
    subplot(2,2,4)
    scatter([PP.Offs,8000], [mean(response{i,1}(idx ==1,:),1) 1],'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1])
    hold on
    h=errorbar([PP.Offs,8000],[mean(response{i,1}(idx ==1,:),1) 1],[std(response{i,1}(idx ==1,:),[],1) 0]./sqrt(sum(idx ==1)),'Color',[1 0 0]);
    h.LineStyle = 'none';
    scatter([PP.Offs,8000], [mean(response{i,1}(idx ==2,:),1) 1],'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[1 1 1])
    hold on
    h=errorbar([PP.Offs,8000],[mean(response{i,1}(idx ==2,:),1) 1],[std(response{i,1}(idx ==2,:),[],1) 0]./sqrt(sum(idx ==2)),'Color',[0 0 1]);
    h.LineStyle = 'none';
    ylim([0 1.2])
    set(gca,'Xtick',[0 2000 4000 8000])
    set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
    ylabel('Normalized FR')
    xlabel('Inter-stimulus interval (s)')
    
end 
%% plot the all cortical neurons broad spiking vs narrow spiking
% concate all cortical areas
pt_ratio_cat = [];
pt_ms_cat = [];
waveform_cat = [];
response_cat = [];

for i = 1:4
    pt_ratio_cat = cat(1,pt_ratio_cat,pt_ratio{i,1});
    pt_ms_cat = cat(1,pt_ms_cat,pt_ms{i,1});
    waveform_cat = cat(2,waveform_cat,waveform{i,1});
    response_cat = cat(1,response_cat,response{i,1});
end 

figure
subplot(2,2,1)
scatter(pt_ratio_cat, pt_ms_cat)
hold on
data = [pt_ratio_cat, pt_ms_cat];
options = statset('Display','final');

obj = fitgmdist(data,2);
h = ezcontour(@(x,y)pdf(obj,[x y]),[0 0.6],[0 1.2]);
xlabel('peak:trough height ratio')
ylabel('trough to peak (ms)')

idx = cluster(obj,data);
cluster1 = data(idx ==1,:);
cluster2 = data(idx ==2,:);
set(gca,'TickDir','out')
axis square


subplot(2,2,2)
h1 = scatter(cluster1(:,1),cluster1(:,2),20,'ko');
hold on
h2 = scatter(cluster2(:,1),cluster2(:,2),20,'ro');
legend([h1 h2],'Broad spiking','Narrow spiking','Location','NW')
xlim([0 0.6])
ylim([0 1.2])
xlabel('peak:trough height ratio')
ylabel('trough to peak (ms)')
set(gca,'TickDir','out')
axis square

subplot(2,2,3)
wavedege = 0:1./30000:60./30000;
plot( wavedege*1000,waveform_cat(:,idx ==1)./(-repmat(min(waveform_cat(:,idx ==1),[],1),61,1)),'k')
hold on
plot( wavedege*1000,waveform_cat(:,idx ==2)./(-repmat(min(waveform_cat(:,idx ==2),[],1),61,1)),'r')
xlabel('ms')
ylabel('normalized amplitude')
% plot the adaptation broad spiking vs narrow spiking neurons for all
% layaers
set(gca,'TickDir','out')
axis square
subplot(2,2,4)
% scatter([PP.Offs,8000], [mean(response_cat(idx ==1,:),1) 1],'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1])

h=errorbar([PP.Offs,8000],[mean(response_cat(idx ==1,:),1) 1],[std(response_cat(idx ==1,:),[],1) 0]./sqrt(sum(idx ==1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.Marker = 'o';
temp = [];
temp = response_cat(idx ==1,:);
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);
hold on
DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[0 0 0]);
text(3000,0.6,[': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' R^2:' sprintf('%.3f',R_square)],'color',[0 0 0])



hold on
h=errorbar([PP.Offs,8000],[mean(response_cat(idx ==2,:),1) 1],[std(response_cat(idx ==2,:),[],1) 0]./sqrt(sum(idx ==2)),'Color',[1 0 0]);
h.LineStyle = 'none';
h.Marker = 'o';

temp = [];
temp = response_cat(idx ==2,:);
[tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);
hold on
DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,[1 0 0]);
text(3000,0.5,[': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' R^2:' sprintf('%.3f',R_square)],'color',[1 0 0])

ylim([0 1])
set(gca,'Xtick',[0 2000 4000 8000])
set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
ylabel('Normalized FR')
xlabel('Inter-stimulus interval (s)')
set(gca,'TickDir','out')
axis square
%% try to plot the putative excitatory neuron vs inhibitory neuron
for i = 2%:length(areas)
    figure
    subplot(2,2,1)
    scatter(pt_ratio{i,1}, pt_ms{i,1})
    hold on
    data = [pt_ratio{i,1}, pt_ms{i,1}];
    options = statset('Display','final');
    
    obj = fitgmdist(data,2);
    h = ezcontour(@(x,y)pdf(obj,[x y]),[0.1 0.5],[0 1.2]);
    xlabel('peak:trough height ratio')
    ylabel('trough to peak (ms)')
    
    idx = cluster(obj,data);
    cluster1 = data(idx ==1,:);
    cluster2 = data(idx ==2,:);
    
    
    subplot(2,2,2)
    h1 = scatter(cluster1(:,1),cluster1(:,2),10,'ro');
    hold on
    h2 = scatter(cluster2(:,1),cluster2(:,2),10,'bo');
     legend([h1 h2],'Broad spiking','Narrow spiking','Location','NW')
     xlim([0.1 0.5])
     ylim([0 1.2])
    xlabel('peak:trough height ratio')
    ylabel('trough to peak (ms)')
    
    subplot(2,2,3)
    wavedege = 0:1./30000:60./30000;
    plot( wavedege*1000,waveform{i,1}(:,idx ==1)./(-repmat(min(waveform{i,1}(:,idx ==1),[],1),61,1)),'r')
    hold on
    plot( wavedege*1000,waveform{i,1}(:,idx ==2)./(-repmat(min(waveform{i,1}(:,idx ==2),[],1),61,1)),'b')
    xlabel('ms')
    ylabel('normalized amplitude')
    % plot the adaptation broad spiking vs narrow spiking neurons for all
    % layaers
    subplot(2,2,4)
    scatter([PP.Offs,8000], [mean(response{i,1}(idx ==1,:),1) 1],'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1])
    hold on
    h=errorbar([PP.Offs,8000],[mean(response{i,1}(idx ==1,:),1) 1],[std(response{i,1}(idx ==1,:),[],1) 0]./sqrt(sum(idx ==1)),'Color',[1 0 0]);
    h.LineStyle = 'none';
    scatter([PP.Offs,8000], [mean(response{i,1}(idx ==2,:),1) 1],'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[1 1 1])
    hold on
    h=errorbar([PP.Offs,8000],[mean(response{i,1}(idx ==2,:),1) 1],[std(response{i,1}(idx ==2,:),[],1) 0]./sqrt(sum(idx ==2)),'Color',[0 0 1]);
    h.LineStyle = 'none';
    ylim([0 1.2])
    set(gca,'Xtick',[0 2000 4000 8000])
    set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
    ylabel('Normalized FR')
    xlabel('Inter-stimulus interval (s)')
    
end 
%% plot one area superficial vs deep
figure
for i =1%1: length(areas)
    
    
    
    h=errorbar([PP.Offs,8000],[mean(response{i,1}(layers{i,1}<4,:),1) 1],[std(response{i,1}(layers{i,1}<4,:),[],1) 0]./sqrt(size(response{i,1},1)),'Color',[0 0 0]);
    h.LineStyle = 'none';
    h.Marker = 'o';
    
    
    hold on
     
    hold on
    h=errorbar([PP.Offs,8000],[mean(response{i,1}(layers{i,1}==4,:),1) 1],[std(response{i,1}(layers{i,1}==4,:),[],1) 0]./sqrt(size(response{i,1},1)),'Color',[0 0 0]);
    h.LineStyle = 'none';
    
    
    ylim([0 1.2])
    set(gca,'Xtick',[0 2000 4000 8000])
    set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
    ylabel('Normalized FR')
    xlabel('Inter-stimulus interval (s)')
    
    
    
    
end




%% plot the single exponential fit for superficial layers
figure
%  subplot(2,2,1)
for i = 1%:length(areas)
    
    if i==1
         % V1 3 mean layer2/3, SC 3 means SO
        temp = response{i,1}(layers{i,1}==3,:);
    elseif i==5
        temp = response{i,1}(layers{i,1}~=4,:); % for SC superficial layers 
        
    else
        temp = response{i,1}(layers{i,1}==2,:);
        
    end
    
    
    % get the fit 
    [tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1])

   
    
    fitsingle_super(i,:) =  [tao,A,sse,R_square];
    h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',colors(i,:));
    h.LineStyle = 'none';
    
    hold on
    scatter([PP.Offs,8000], [mean(temp,1) 1],'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1])
    
    DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,colors(i,:));
    
  
    
    
    ylim([0 1.4])
    set(gca,'Xtick',[0 2000 4000 8000])
    set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
    ylabel('Normalized FR')
    xlabel('Inter-stimulus interval (s)')
    text(3000,0.6-i*0.1,[areas{i},': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',colors(i,:))
    
    title('Superficial layers only; Fit:single exponential')
    
end

% one way anova test 
% Temp = [repmat(1,30,1),temp];
% [p,tbl,stats] = anova1(Temp)
% c= multcompare(stats);

%% plot difference across deep layers

subplot(2,2,2)
for i = 1:length(areas)
    
    if i==1
         % V1 3 mean layer2/3, SC 3 means SO
        temp = response{i,1}(layers{i,1}>3,:);
    elseif i==5
        temp = response{i,1}(layers{i,1}==4,:); % for SC superficial layers 
        
    else
        temp = response{i,1}(layers{i,1}>2,:);
        
    end
    
    
    % get the fit 
    [tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1])

   
    
    fitsingle_deep(i,:) =  [tao,A,sse,R_square];
    h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',colors(i,:));
    h.LineStyle = 'none';
    
    hold on
    scatter([PP.Offs,8000], [mean(temp,1) 1],'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:))
    
    DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,colors(i,:));
    
  
    
    
    ylim([0 1.2])
    set(gca,'Xtick',[0 2000 4000 8000])
    set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
    ylabel('Normalized FR')
    xlabel('Inter-stimulus interval (s)')
    text(3000,0.6-i*0.1,[areas{i},': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',colors(i,:))
    
    title('Deep layers only; Fit:single exponential')
    
end

%% plot difference combine all layers

 figure
for i = 1:5%length(areas)
    
   
        temp = response{i,1};
        
  
    
    
    % get the fit 
    [tao,A,sse,R_square] = SingleExponFit([PP.Offs,10000],[mean(temp,1),1]);

   
    
    fitsingle_all(i,:) =  [tao,A,sse,R_square];
    h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',colors(i,:));
    h.LineStyle = 'none';
    
    hold on
    scatter([PP.Offs,8000], [mean(temp,1) 1],'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:))
    
    DoubleExponFitPlot([0:10:8000],tao,1,1,A,0,colors(i,:));
    
  
    
    
    ylim([0 1.1])
    set(gca,'Xtick',[0 2000 4000 8000])
    set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
    ylabel('Normalized FR')
    xlabel('Inter-stimulus interval (s)')
    text(3000,0.6-i*0.1,[areas{i},': n = ' num2str(size(temp,1)) ': Tao:' sprintf('%.2f',tao./1000) ' A:' sprintf('%.2f',A)  ' R^2:' sprintf('%.3f',R_square)],'color',colors(i,:))
    
    title('All layers combined; Fit:single exponential')
    
end


%% plot define layers 
% shankchannel=[2,9,3,8,4,7,14,6,5,13,10,1,11,12,15,16 ; 31 19 29 20 28 26 27 24 25 23 30 32 22 21 18 17];
% e_distance(1,:) = [25:50:775]; % shank1
% e_distance(2,:) = [0:50:750]; % shank2
% baseline=0.2;
% edges = [-baseline:1/10000:(0.1+0.4)];
% i=1;
% figure
% for i_exp = 1:length(Gdata)
%     
% load(cell2mat(Gtext(i_exp+1,4)))
% subplot(2,2,i)
% i=i+1;
%  if i_exp ==2 || i_exp ==5
%     shank= 1;
%  else
%      shank=2;
%  end
%     d=-diff(squeeze(LFP.mean(shank,:,:)),2,1);% approximate second spatial derivative
%     d_i = interp2(d); % doubles the samples in every direction
%     d_f = d_i(:,1:2:size(d_i,2)); % downsamples the times again
%     test = linspace(min(e_distance(shank,:)),max(e_distance(shank,:)),size(d_i,1));
%       
% 
%     pcolor(edges,(-test),d_f);
%   
%     xlim([-0.1 0.3])
%     shading interp;
%     colormap(jet)
%     caxis manual
%     
%     % caxis([-200 200]) % for V1
%     caxis([-150 150]) % for HVAs
%     colorbar
%     hold on
%      hline(-LFP.layers,'k:')
%    
%     title(['shank:', num2str(shank)])
% %     hold on
% %     scatter(0.05,-400,'MarkerEdgeColor',[1 1 1])
% 
%   
% 
%     
% end

%% scatter plot 
minvalue=[-2]
binwidth = 0.020;
low=-5;
high=80;
figure

i=1;

for i_layers= [2 3 4]%[3 4 5]
    
    subplot(3,2,i)
    temp=[];
    temp = FR(layers==i_layers,:)./binwidth;
    h= scatter(temp(:,1),temp(:,6),'MarkerEdgeColor',[0 0 0]);
    h.SizeData = 15;
    xlim([low high])
    ylim([low high])
    hold on
    plot([-5 80],[-5 80],'k:')
    xlabel('FR-250ms')
    ylabel('FR-baseline')
    if i_layers ==2 %3
        title(['layer' '2/3'])
    else
        title(['layer' num2str(i_layers+1) '?'])
        %         title(['layer' num2str(i_layers)])
    end
    i=i+1;
    subplot(3,2,i)
    
    h=scatter(temp(:,3),temp(:,6),'MarkerEdgeColor',[0 0 0]);
    h.SizeData = 15;
    
    xlim([low high])
    ylim([low high])
    hold on
    plot([-5 80],[-5 80],'k:')
    xlabel('FR-1s')
    ylabel('FR-baseline')
    if i_layers ==2 %3
        title(['layer' '2/3'])
    else
        title(['layer' num2str(i_layers+1) '?'])
        %         title(['layer' num2str(i_layers)])
    end
    i=i+1;
    
end

%% plot the normalized FR over intervals 

 figure
hold on

h=errorbar([PP.Offs,8000],[mean(response,1) 1],[std(response,[],1) 0]./sqrt(size(response,1)),'Color',[1 0 0]);
h.LineStyle = 'none';

hold on
scatter([PP.Offs,8000], [mean(response,1) 1],'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1])


ylim([0 1.2])
set(gca,'Xtick',[0 2000 4000 8000])
set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
ylabel('Normalized FR')
xlabel('Inter-stimulus interval (s)')
text(100,1.2,['n = ' num2str(size(response,1))],'color',[1 0 0])
title('Layer 2/3 and layer 4?')
%% plot the difference between layers
Colorchoice = [0 0 0; 0 0 0;0 0 0];
figure

for i_layers = [2 3 4]%[3 4 5]
    subplot(2,2,i_layers-1)
    temp = response(layers==i_layers,:);
    h=errorbar([PP.Offs,8000],[mean(temp,1) 1],[std(temp,[],1) 0]./sqrt(size(temp,1)),'Color',[0 0 0]);
    h.LineStyle = 'none';
    h.Color = Colorchoice(i_layers-1,:);
    
    hold on
    h=scatter([PP.Offs,8000], [mean(temp,1) 1],'MarkerEdgeColor',Colorchoice(i_layers-1,:),'MarkerFaceColor',[1 1 1])
    
    
   ylim([0 1.4])
    set(gca,'Xtick',[0 2000 4000 8000])
    set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
    ylabel('Normalized FR')
    xlabel('Inter-stimulus interval (s)')
    text(5000,0.4,['n = ' num2str(size(temp,1))],'Color',Colorchoice(i_layers-1,:))
    if i_layers ==2 %3
        title(['layer' '2/3'])
    else
        title(['layer' num2str(i_layers+1) '?'])
        %         title(['layer' num2str(i_layers)])
    end
end

subplot(2,2,4)
h=errorbar([PP.Offs,8000],[mean(response,1) 1],[std(response,[],1) 0]./sqrt(size(response,1)),'Color',[0 0 0]);
h.LineStyle = 'none';

hold on
scatter([PP.Offs,8000], [mean(response,1) 1],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1])


ylim([0 1.2])
set(gca,'Xtick',[0 2000 4000 8000])
set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
ylabel('Normalized FR')
xlabel('Inter-stimulus interval (s)')
text(5000,0.4,['n = ' num2str(size(response,1))])
title('All layers')