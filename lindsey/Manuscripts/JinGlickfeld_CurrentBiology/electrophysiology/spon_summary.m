%need to separate out the power
[Gdata, Gtext, Graw] = xlsread('Z:\All_Staff\home\miaomiao\Analysis\electrode\ducumentation\optogenetic effects.xlsx',6);
Dir = Graw(2:end,11);
LED_power = cell2mat(Graw(2:end,10));
Background = Graw(2:end,8);
opsin = Graw(2:end,7);
LED_center =  cell2mat(Graw(2:end,9));
Center_area = Graw(2:end,12); 
power = unique(LED_power);
shankchannel=[7 14 4 6 8 5 3 13;11 1 12 10 15 2 16 9;32 22 30 21 31 18 19 17;27 26 24 28 25 20 23 29];
shanks = [1:4];
electrodes=8;
% pull out the cells that are dirven or inhibit by led stimultion

excite.recovery = [];
inhibit.recovery = [];
excite.all = [];
inhibit.all = [];
excite.template = [];
inhibit.template = [];
excite.basebin = [];
inhibit.basebin = [];
excite.ALLbasebin = [];
inhibit.ALLbasebin = [];
excite.sebasebin = [];
inhibit.sebasebin = [];

excite.basebinTrl = {};
inhibit.basebinTrl = {};

excite.template_min = [];
inhibit.template_min = [];
excite.wavestats.pt_ms = [];
inhibit.wavestats.pt_ms = [];
excite.wavestats.pt_ratio = [];
inhibit.wavestats.pt_ratio = [];
LED_spreadcell = cell(8,7);
LED_foldcell =cell(8,4);
excite.test = cell(8,7);
excite.test_fold = cell(8,4);
inhibit.test = cell(8,7);
inhibit.test_fold = cell(8,4);

LED_spread = NaN(8,7,length(LED_center)); % treat each experiments as a value
LED_fold = NaN(8,4,length(LED_center)); % fold the led for symetricity
adjust = [3,2,1,0];
areas = {'V1','LM','AL','PM'}; 
C_pchg = []; 
C_area = {}; 
all.inhibitPSTH = [];
all.inhibitbase = [];
for i = 1:size( Dir,1)
    load(Dir{i,1})
   
    % collapsed p change in the center shanks
    idx =  ismember(spon.spikes.chID,shankchannel(LED_center(i),:))'&(~spon.stattest.excite);
    C_pchg = [C_pchg spon.stattest.pchange(idx)]; 
    C_area = [C_area repmat(Center_area(i),1,sum(idx))]; 
    % get the response traces
    all.inhibitPSTH = [all.inhibitPSTH;spon.psth.mean(idx,:) ];
    all.inhibitbase = [all.inhibitbase; mean(spon.psth.basebin(idx,:),2)];
    
    % criterion pairttest tail significant and in the led center shank
    idx =ismember(spon.spikes.chID,shankchannel(LED_center(i),:))'& spon.stattest.excite&(spon.stattest.pchange>0); % 
    excite.recovery =  [excite.recovery; spon.psth.off.mean(idx,:)];
    excite.all = [excite.all;spon.psth.mean(idx,:) ];
    excite.basebin = [excite.basebin; mean(spon.psth.basebin(idx,:),2)];
    excite.ALLbasebin = [excite.ALLbasebin; mean(spon.psth.all.basebin(idx,:),2)];
    %excite.sebasebin = [excite.sebasebin; mean([spon.psth.basebin(idx,:) spon.psth.postled.basebin(idx,:)],2)];
    excite.sebasebin = [excite.sebasebin; mean([mean(spon.psth.basebin(idx,:),2) mean(spon.psth.postled.basebin(idx,:),2)],2)];
    excite.basebinTrl{i} = spon.psth.basebin(idx,:);
  
    excite.template = [excite.template; spon.spikes.template_big(:,idx)'];
    excite.template_min = [excite.template_min; min(spon.spikes.template_big(:,idx)',[],2)];
    excite.wavestats.pt_ms = [excite.wavestats.pt_ms; spon.spikes.wave_stats.pt_ms(idx)];
    excite.wavestats.pt_ratio = [ excite.wavestats.pt_ratio; spon.spikes.wave_stats.pt_ratio(idx)];
    
     
    
    idx = ismember(spon.spikes.chID,shankchannel(LED_center(i),:))'& spon.stattest.inhibit&(spon.stattest.pchange<0); % (spon.stattest.pchange<0)
    
    inhibit.recovery =  [inhibit.recovery; spon.psth.off.mean(idx,:)];
    inhibit.all = [inhibit.all;spon.psth.mean(idx,:) ];
    inhibit.basebin = [inhibit.basebin; mean(spon.psth.basebin(idx,:),2)];
    inhibit.ALLbasebin = [inhibit.ALLbasebin; mean(spon.psth.all.basebin(idx,:),2)];
    %inhibit.sebasebin = [inhibit.sebasebin; mean([spon.psth.basebin(idx,:) spon.psth.postled.basebin(idx,:)],2)];
    inhibit.sebasebin = [inhibit.sebasebin; mean([mean(spon.psth.basebin(idx,:),2) mean(spon.psth.postled.basebin(idx,:),2)],2)];
    inhibit.basebinTrl{i} = spon.psth.basebin(idx,:);
    inhibit.template = [inhibit.template; spon.spikes.template_big(:,idx)'];
    inhibit.template_min = [inhibit.template_min; min(spon.spikes.template_big(:,idx)',[],2)];
    inhibit.wavestats.pt_ms = [inhibit.wavestats.pt_ms; spon.spikes.wave_stats.pt_ms(idx)];
    inhibit.wavestats.pt_ratio = [inhibit.wavestats.pt_ratio; spon.spikes.wave_stats.pt_ratio(idx)];
    
    
%     %check for waveform shapes
%     if sum(idx)>0
%         figure
%        
%         plot(spon.spikes.template_big(:,idx))
%          title(['exp:' num2str(i)])
%     end 
%     % check for waveform shapes
    
    shank_plus = adjust(LED_center(i));
    for shank=1:4%shanks       
        for j=1:electrodes
            elec=shankchannel(shank,j);
            idx_temp = (spon.spikes.chID==elec )& (~spon.stattest.excite)';         
            if sum(idx_temp)>0
               LED_spread(j,shank+shank_plus,i) = mean(spon.stattest.pchange(idx_temp));
               LED_spreadcell{j,shank+shank_plus} = [LED_spreadcell{j,shank+shank_plus}  spon.stattest.pchange(idx_temp)];
               
              

            end
            idx_temp = spon.spikes.chID==elec ;
            excite.test{j,shank+shank_plus} = [excite.test{j,shank+shank_plus} spon.stattest.excite(idx_temp)];
            inhibit.test{j,shank+shank_plus} = [inhibit.test{j,shank+shank_plus} spon.stattest.inhibit(idx_temp)];
            idx_temp = [];
        end
        
    end
    if LED_center(i)==4
        LED_fold(1:8,1,i) = squeeze(LED_spread(1:8,4,i));
        LED_fold(1:8,2,i) = squeeze(LED_spread(1:8,3,i));
        LED_fold(1:8,3,i) = squeeze(LED_spread(1:8,2,i));
        LED_fold(1:8,4,i) = squeeze(LED_spread(1:8,1,i));
        
        
    elseif LED_center(i)==1
        LED_fold(1:8,1:4,i) = LED_spread(1:8,4:7,i);
    elseif LED_center(i)==2
        LED_fold(1:8,1,i) = squeeze(LED_spread(1:8,4,i));
        LED_fold(1:8,2,i) = nanmean(LED_spread(1:8,[3 5],i),2);
        LED_fold(1:8,3,i) = squeeze(LED_spread(1:8,6,i));
    elseif LED_center(i)==3   
        LED_fold(1:8,1,i) = squeeze(LED_spread(1:8,4,i));
        LED_fold(1:8,2,i) = nanmean(LED_spread(1:8,[3 5],i),2);
        LED_fold(1:8,3,i) = squeeze(LED_spread(1:8,2,i));
    end 
 
    
    
    
end


LED_foldcell(1:8,1) = LED_spreadcell(:,4);
LED_foldcell(1:8,2) = cellfun(@(x,y) cat(2,x,y),LED_spreadcell(:,3),LED_spreadcell(:,5),'un',0);
LED_foldcell(1:8,3) = cellfun(@(x,y) cat(2,x,y),LED_spreadcell(:,2),LED_spreadcell(:,6),'un',0);
LED_foldcell(1:8,4) = cellfun(@(x,y) cat(2,x,y),LED_spreadcell(:,1),LED_spreadcell(:,7),'un',0);


inhibit.test_fold(1:8,1) = inhibit.test(:,4);
inhibit.test_fold(1:8,2) = cellfun(@(x,y) cat(2,x,y),inhibit.test(:,3),inhibit.test(:,5),'un',0);
inhibit.test_fold(1:8,3) = cellfun(@(x,y) cat(2,x,y),inhibit.test(:,2),inhibit.test(:,6),'un',0);
inhibit.test_fold(1:8,4) = cellfun(@(x,y) cat(2,x,y),inhibit.test(:,1),inhibit.test(:,7),'un',0);

excite.test_fold(1:8,1) = excite.test(:,4);
excite.test_fold(1:8,2) = cellfun(@(x,y) cat(2,x,y),excite.test(:,3),excite.test(:,5),'un',0);
excite.test_fold(1:8,3) = cellfun(@(x,y) cat(2,x,y),excite.test(:,2),excite.test(:,6),'un',0);
excite.test_fold(1:8,4) = cellfun(@(x,y) cat(2,x,y),excite.test(:,1),excite.test(:,7),'un',0);

edges = spon.edges;
offedges = spon.offedges;
%% get the fit for each experiments 
spread_exp = squeeze(nanmean(LED_fold,1));
tau_exp = NaN(1,size(Dir,1));
colors = gray(16);
figure
% get the fit for each experiments
for i= 1:size(Dir,1)
    temp = [];
    temp = spread_exp(:,i)';
    whichnan = isnan(temp);
    % nan values cannot be in the center shank
    if (whichnan(1)~=1) || sum(whichnan)==0 
     xvalues = 0:0.4:1.2;
     xvalues(whichnan)=[];
     temp(whichnan) = [];
     
    f=fit(xvalues',-temp','exp1');
    DoubleExponFitPlot(0:0.1:1.4,-1./f.b,1,0,-f.a,0,colors(i,:));
    hold on
    tau_exp(i)=-1./f.b;
    
    end

end 


%%
figure
subplot(1,2,1)
for i_power = 1:length(power)
    
    idx = [];
    idx = LED_power==power(i_power);
    scatter(repmat(power(i_power),1,sum(idx)),tau_exp(idx))
    hold on
end 
ylabel('tau(mm)')
xlabel('LED power(mW)')
subplot(1,2,2)
for i_power = 1:length(power)
    
    idx = [];
    idx = LED_power==power(i_power);
    scatter(repmat(power(i_power),1,sum(idx)),-spread_exp(1,idx))
    hold on
end 
ylabel('(Base-LED)/Base')
xlabel('LED power(mW)')

%% plot the spread by each units, BASE-LED/BASE
% plot the fraction of cells shows significant inhibited...
  
LED_spread_mean = cell2mat(cellfun(@mean,LED_foldcell,'un',0));
n_exp = cell2mat(cellfun(@length,LED_foldcell,'un',0));
LED_spread_sem = cell2mat(cellfun(@(x) std(x)./sqrt(length(x)),LED_foldcell,'un',0));

for i = 1:4
    tanspread{i} = cell2mat(LED_foldcell(:, i)');
end


tan_spread_mean = cell2mat(cellfun(@mean,tanspread,'un',0));
tan_n = cell2mat(cellfun(@length,tanspread,'un',0));
tan_spread_sem = cell2mat(cellfun(@(x) std(x)./sqrt(length(x)),tanspread,'un',0));

figure
subplot(2,5,1:3)


imagesc(LED_spread_mean,[-1 1])

colormap(flipud(brewermap([],'RdBu')))
colorbar
set(gca,'TickDir','out','XTick',1:1:4,'XTickLabel',{'0','0.4','0.8','1.2'},'YTick',1:1:8,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'})
ylabel('Distance: superficial to deep (mm)')
xlabel('Distance to LED center (mm)')

subplot(2,5,6:7)
errorbar(0:0.1:0.7,-LED_spread_mean(:,1),LED_spread_sem(:,1),'k')
hold on
scatter(0:0.1:0.7,-LED_spread_mean(:,1),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
ylim([-0.2 1])
xlim([-0.05 0.75])
set(gca,'TickDir','out','XTick',0:0.1:0.7)
ylabel('Ratio of inhibition')
xlabel('Distance: superficial to deep (mm)')
subplot(2,5,9:10)
h = errorbar(0:0.4:1.2,-tan_spread_mean,tan_spread_sem,'k');
h.LineStyle = 'none';
% get the single exponential fit
hold on
f=fit((0:0.4:1.2)',-tan_spread_mean','exp1');
DoubleExponFitPlot(0:0.1:1.4,-1./f.b,1,0,-f.a,0,[0 0 0]);
text(0.1,0.9,['Tau = ' num2str(-1./f.b)])
hold on
scatter(0:0.4:1.2,-tan_spread_mean,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])


ylim([0 1])
xlim([-0.05 1.3])
set(gca,'TickDir','out','XTick',0:0.4:1.2)
ylabel('Ratio of inhibition')
xlabel('Distance to LED center (mm)')
%% get the fit with confidence intervals 
tan_spread_mean = cell2mat(cellfun(@mean,tanspread,'un',0));
tan_n = cell2mat(cellfun(@length,tanspread,'un',0));
rng(0)
n_boot = 1000; 
spread_boot = nan(n_boot,4); 
T_boot = nan(n_boot,1); 
for i_boot= 1:n_boot
    for i_shank = 1:4
       spread_boot(i_boot,i_shank) =  mean(randsample(tanspread{i_shank},tan_n(i_shank),1));        
    end 
    F{i_boot}=fit((0:0.4:1.2)',-spread_boot(i_boot,:)','exp1');
    T_boot(i_boot) = -1./F{i_boot}.b; 
    
end 

% confidence intervals
T_sort = [];
T_sort = sort(T_boot);
T_5 = T_sort(round(0.025*n_boot));
T_95 = T_sort(round(0.975*n_boot));
   


%% stats
num_cell = n_exp(:,1);
group = [];
 
data = [];
for i=1:8

    data = [data;-LED_foldcell{i,1}'];
   
    group = [group;ones(num_cell(i),1)+i ]; 
end 

[p1,tbl,stats] = kruskalwallis(data,group);
   b= multcompare(stats);
   %% stats
num_cell = tan_n;
group = [];
 
data = [];
for i=1:4

    data = [data;-tanspread{1,i}'];
   
    group = [group;ones(num_cell(i),1)+i ]; 
end 

[p1,tbl,stats] = kruskalwallis(data,group);
   b= multcompare(stats);
%% plot the spread of perecent inhibited
% calculate num of inhibit and excite and other
inhibit_num = cell2mat(cellfun(@sum,inhibit.test_fold,'un',0));
excite_num =  cell2mat(cellfun(@sum,excite.test_fold,'un',0));
total_num = cell2mat(cellfun(@length,excite.test_fold,'un',0));
% percent suppression exclude the excited cells 
tot_noexcite = total_num - excite_num; 
[phat,pci] = binofit(inhibit_num(:),tot_noexcite(:));
p_inhibit = reshape(phat,[8,4]);
pci_inhibit = reshape(pci,[8,4,2]);
% spread of percent inhibit

[tan_phat,tan_pci] = binofit(sum(inhibit_num,1),sum(tot_noexcite,1));

figure
subplot(2,5,1:3)


imagesc(p_inhibit,[-1 1])

colormap(brewermap([],'RdBu'))
colorbar
set(gca,'TickDir','out','XTick',1:1:4,'XTickLabel',{'0','0.4','0.8','1.2'},'YTick',1:1:8,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'})
ylabel('Distance: superficial to deep (mm)')
xlabel('Distance to LED center (mm)')

subplot(2,5,6:7)
for i = 1:8
    plot([0.1*i-0.1 0.1*i-0.1],pci(i,:),'k')
    hold on
end 

hold on
scatter(0:0.1:0.7,p_inhibit(:,1),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
ylim([-0.2 1.2])
xlim([-0.05 0.75])
set(gca,'TickDir','out','XTick',0:0.1:0.7)
ylabel('Percent of significant inhibit')
xlabel('Distance: superficial to deep (mm)')
subplot(2,5,9:10)
for i = 1:4
    plot([0.4*i-0.4 0.4*i-0.4],tan_pci(i,:),'k')
    hold on
end 
% get the single exponential fit

f=fit((0:0.4:1.2)',tan_phat','exp1');
DoubleExponFitPlot(0:0.1:1.4,-1./f.b,1,0,-f.a,0,[0 0 0]);
text(0.1,0.9,['Tau = ' num2str(-1./f.b)])
hold on
scatter(0:0.4:1.2,tan_phat,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])


ylim([0 1])
xlim([-0.05 1.3])
set(gca,'TickDir','out','XTick',0:0.4:1.2)
ylabel('Percent of significant inhibit')
xlabel('Distance to LED center (mm)')
%% plot the spread by experiments
% plot the fraction of cells shows significant inhibited...
LED_spread_mean = nanmean(LED_fold,3);
n_exp = sum(~isnan(LED_fold),3);
LED_spread_sem = nanstd(LED_fold,[],3)./sqrt(n_exp);

tan_spread = squeeze(nanmean(LED_fold,1));
tan_spread_mean = nanmean(tan_spread,2);
tan_n = sum(~isnan(tan_spread),2);
tan_spread_sem = nanstd(tan_spread,[],2)./sqrt(tan_n);

figure
subplot(1,2,1)
errorbar(0:0.1:0.7,-LED_spread_mean(:,1),LED_spread_sem(:,1),'k')
hold on
scatter(0:0.1:0.7,-LED_spread_mean(:,1),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
ylim([-0.2 1])
xlim([-0.05 0.75])
set(gca,'TickDir','out','XTick',0:0.1:0.7)
ylabel('Ratio of inhibition')
xlabel('Distance: superficial to deep (um)')
subplot(1,2,2)
h = errorbar(0:0.4:1.2,-tan_spread_mean,tan_spread_sem,'k');
h.LineStyle = 'none';
% get the single exponential fit
hold on
f=fit((0:0.4:1.2)',-tan_spread_mean,'exp1');
DoubleExponFitPlot(0:0.1:1.4,-1./f.b,1,0,-f.a,0,[0 0 0]);

hold on
scatter(0:0.4:1.2,-tan_spread_mean,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])


ylim([0 1])
xlim([-0.05 1.3])
set(gca,'TickDir','out','XTick',0:0.4:1.2,'YTick',0:0.2:1)
ylabel('Ratio of inhibition')
xlabel('Distance to LED center (um)')

%% plot all the cells excited/inhibit by led stimulation
binwidth = spon.binwidth;
LED_duration = spon.LED_duration;
figure
subplot(2,2,1)
shadedErrorBar_ch(edges, mean(excite.all,1)./binwidth, (std(excite.all,[],1)./binwidth)./sqrt(size(excite.all,1)), {'Color',[1 0 0.4],'markerfacecolor',[1 0 0.4]})
hold on
vline(0,'k:')
vline(LED_duration,'k:')
title(['Cells excited' '  n= ' num2str(size(excite.all,1))])
xlim([-1 5])
xlabel('Time (S)')
ylabel('Spikes/S')
subplot(2,2,2)
shadedErrorBar_ch(edges, mean(inhibit.all,1)./binwidth, (std(inhibit.all,[],1)./binwidth)./sqrt(size(inhibit.all,1)), {'Color',[0 0.4 1],'markerfacecolor',[0 0.4 1]})
hold on
vline(0,'k:')
vline(LED_duration,'k:')
xlim([-1 5])
xlabel('Time (S)')
ylabel('Spikes/S')
title(['Cells inhibited' '  n= ' num2str(size(inhibit.all,1))])
% plot the normalized values

%%
binwidth = spon.binwidth;
LED_duration = spon.LED_duration;
figure
subplot(1,2,1)
shadedErrorBar_ch(edges,  mean(excite.all./repmat(excite.basebin,1,size(excite.all,2)),1), std(excite.all./repmat(excite.basebin,1,size(excite.all,2)),[],1)./sqrt(size(excite.all,1)), {'Color',[1 0 0.4],'markerfacecolor',[1 0 0.4]})
hold on

title(['Cells excited' '  n= ' num2str(size(excite.all,1))])
ylim([0 30])
xlim([-1.1 5])
hline(1,'k:')
vline(0,'k:')
vline(LED_duration,'k:')

xlabel('Time (S)')
ylabel('Norm. to baseline')
set(gca,'TickDir','out')
axis square
subplot(1,2,2)
shadedErrorBar_ch(edges,  mean(inhibit.all./repmat(inhibit.basebin,1,size(inhibit.all,2)),1), std(inhibit.all./repmat(inhibit.basebin,1,size(inhibit.all,2)),[],1)./sqrt(size(inhibit.all,1)), {'Color',[0 0.4 1],'markerfacecolor',[0 0.4 1]})
hold on

xlim([-1.1 5])
ylim([0 2])
hline(1,'k:')
vline(0,'k:')
vline(LED_duration,'k:')
set(gca,'TickDir','out')
axis square
xlabel('Time (S)')
ylabel('Norm. to baseline')
title(['Cells inhibited' '  n= ' num2str(size(inhibit.all,1))])
%% stat tests, 3.6 for LED
Norm_temp = []; 
%Norm_temp = excite.all./repmat(excite.basebin,1,size(excite.all,2)); 
Norm_temp = inhibit.all./repmat(inhibit.basebin,1,size(inhibit.all,2)); 
[p1,tbl,stats] = kruskalwallis(Norm_temp(:,edges>=0&edges<=3.3)./repmat(Norm_temp(:,find(edges==0)),1,sum(edges>=0&edges<=3.3)));
   b= multcompare(stats);
%% overlay the waveforms
template_edges = 0:1./30:60./30;
figure
subplot(2,2,1)
plot(template_edges, (excite.template./repmat(-excite.template_min,1,size(excite.template,2)))', 'Color', [1 0 0.4])
ylabel('uV')
xlabel('ms')
hold on
vline(1.3,'k:')
vline(1.6,'k:')
ylim([-1 1.2])
axis square
subplot(2,2,2)
plot(template_edges, (inhibit.template./repmat(-inhibit.template_min,1,size(inhibit.template,2)))', 'Color', [0 0.6 1])
hold on
vline(1.3,'k:')
vline(1.6,'k:')
ylabel('uV')
xlabel('ms')
ylim([-1 1.2])
axis square
subplot(2,2,3)
% plot the average waveform
temp = (excite.template./repmat(-excite.template_min,1,size(excite.template,2)))';
temp_mean = mean(temp,2);
temp_sem = std(temp,[],2)./sqrt(size(temp,2));
plot(template_edges,temp_mean,'Color',[1 0 0.4])
hold on
plot(template_edges,temp_mean+temp_sem,'Color',[1 0 0.4],'Linestyle',':')
plot(template_edges,temp_mean-temp_sem,'Color',[1 0 0.4],'Linestyle',':')

temp = (inhibit.template./repmat(-inhibit.template_min,1,size(inhibit.template,2)))';
temp_mean = mean(temp,2);
temp_sem = std(temp,[],2)./sqrt(size(temp,2));
plot(template_edges,temp_mean,'Color',[0 0.6 1])
hold on
plot(template_edges,temp_mean+temp_sem,'Color',[0 0.6 1],'Linestyle',':')
plot(template_edges,temp_mean-temp_sem,'Color',[0 0.6 1],'Linestyle',':')
ylabel('uV')
xlabel('ms')

subplot(2,2,4)

scatter(inhibit.wavestats.pt_ratio, inhibit.wavestats.pt_ms,'MarkerEdgeColor',[0 0.6 1],'MarkerFaceColor',[0 0.6 1],'SizeData',40)
hold on
scatter(excite.wavestats.pt_ratio, excite.wavestats.pt_ms,'MarkerEdgeColor',[1 0 0.4],'MarkerFaceColor',[1 0 0.4],'SizeData',40)
xlabel('peak:trough height ratio')
ylabel('trough to peak (ms)')

%% zoom in the plot the decay

figure
subplot(2,1,1)
shadedErrorBar_ch(offedges,  mean(excite.recovery./repmat(excite.basebin,1,size(excite.recovery,2)),1), std(excite.recovery./repmat(excite.basebin,1,size(excite.recovery,2)),[],1)./sqrt(size(excite.recovery,1)), {'Color',[1 0 0.4],'markerfacecolor',[1 0 0.4]})
hold on
vline(0,'k:')

title(['Cells excited' '  n= ' num2str(size(excite.recovery,1))])
xlim([-0.1 3.9])
ylim([0 2])
hline(1,'k:')
xlabel('Time after LED off (S)')
ylabel('Norm. to baseline')
subplot(2,1,2)
shadedErrorBar_ch(offedges,  mean(inhibit.recovery./repmat(inhibit.basebin,1,size(inhibit.recovery,2)),1), std(inhibit.recovery./repmat(inhibit.basebin,1,size(inhibit.recovery,2)),[],1)./sqrt(size(inhibit.recovery,1)), {'Color',[0 0.4 1],'markerfacecolor',[0 0.4 1]})


xlim([-0.1 3.9])
ylim([0 2])
hold on
vline(0,'k:')
hline(1,'k:')
xlabel('Time after LED off(S)')
ylabel('Norm. to baseline')
title(['Cells inhibited' '  n= ' num2str(size(inhibit.recovery,1))])

%% plot the decay without normalize



figure
subplot(2,1,1)
shadedErrorBar_ch(offedges,  mean(excite.recovery,1)./binwidth, (std(excite.recovery,[],1)./sqrt(size(excite.recovery,1)))./binwidth, {'Color',[1 0 0.4],'markerfacecolor',[1 0 0.4]})
hold on
vline(0,'k:')

title(['Cells excited' '  n= ' num2str(size(excite.recovery,1))])
 xlim([-0.1 3.9])
ylim([0 10])
hline(1,'k:')
xlabel('Time after LED off (S)')
ylabel('Norm. to baseline')
subplot(2,1,2)
shadedErrorBar_ch(offedges,  mean(inhibit.recovery,1)./binwidth, (std(inhibit.recovery,[],1)./sqrt(size(inhibit.recovery,1)))./binwidth, {'Color',[0 0.4 1],'markerfacecolor',[0 0.4 1]})


 xlim([-0.1 3.9])
ylim([0 10])
hold on
vline(0,'k:')
hline(1,'k:')
xlabel('Time after LED off(S)')
ylabel('Norm. to baseline')
title(['Cells inhibited' '  n= ' num2str(size(inhibit.recovery,1))])
%% test rebound excitation
rebound_temp = [];
rebound_temp = mean(inhibit.recovery(:,offedges>=0&offedges<=0.2),2); 
[p,h]=signrank(inhibit.basebin,rebound_temp);

%% recover with new normalization
figure
subplot(2,1,1)
shadedErrorBar_ch(offedges,  mean(excite.recovery./repmat(excite.sebasebin,1,size(excite.recovery,2)),1), std(excite.recovery./repmat(excite.sebasebin,1,size(excite.recovery,2)),[],1)./sqrt(size(excite.recovery,1)), {'Color',[1 0 0.4],'markerfacecolor',[1 0 0.4]})
hold on
vline(0,'k:')

title(['Cells excited' '  n= ' num2str(size(excite.recovery,1))])
%  xlim([-0.1 3.5])
ylim([0 2])
hline(1,'k:')
xlabel('Time after LED off (S)')
ylabel('Norm. to baseline')
subplot(2,1,2)
shadedErrorBar_ch(offedges,  mean(inhibit.recovery./repmat(inhibit.sebasebin,1,size(inhibit.recovery,2)),1), std(inhibit.recovery./repmat(inhibit.sebasebin,1,size(inhibit.recovery,2)),[],1)./sqrt(size(inhibit.recovery,1)), {'Color',[0 0.4 1],'markerfacecolor',[0 0.4 1]})


%xlim([-0.1 3.5])
ylim([0 2])
hold on
vline(0,'k:')
hline(1,'k:')
xlabel('Time after LED off(S)')
ylabel('Norm. to baseline')
title(['Cells inhibited' '  n= ' num2str(size(inhibit.recovery,1))])



%% stat comparison between two methods 
PV=load('Z:\All_Staff\home\miaomiao\labmeeting\HVAs\LED spread\PV_stats.mat'); 
VGAT = load('Z:\All_Staff\home\miaomiao\labmeeting\HVAs\LED spread\vgat_stats.mat'); 
% using 2 way anova
num_cell = PV.tan_n;
group1 = [];

data1 = [];
for i=1:4

    data1 = [data1;-PV.tanspread{1,i}'];
   
    group1= [group1;ones(num_cell(i),1)+i ]; 
end 

num_cell = VGAT.tan_n;
group2 = [];

data2 = [];
for i=1:4

    data2 = [data2;-VGAT.tanspread{1,i}'];
   
    group2= [group2;ones(num_cell(i),1)+i ]; 
end 

Data = [data1; data2]; 
g1 = [repmat('P', sum(PV.tan_n),1); repmat('V', sum(VGAT.tan_n),1)]; 
g2 = [group1; group2]; 
[p1,tbl,stats] = anovan(Data,{g1 g2}, 'model', 'interaction','varnames',{'method','distance'});

%% 
[p,h]=ranksum(-PV.tanspread{1,1},-VGAT.tanspread{1,1});

%% plot the center shank aross areas


for i_area = 1:4
    idx = [];
    idx = strcmp(C_area, areas{i_area}); 
    C_change{i_area} = C_pchg(idx); 
    
end 
%% 
Colors(1,:) = [0 0 0];
Colors(2:4,:) = lines(3); 
n_cell = cellfun(@length, C_change); 
figure
plotSpread(C_change, 'distributionColors', Colors);
hold on

errorbar(1:4, cellfun(@mean,C_change ), cellfun(@std,C_change)./sqrt(n_cell), 'LineStyle','none','Marker','o','Color',[0 0 0])
ylabel('Norm. suppression')
xlim([0.5 4.5])
set(gca,'XTicklabel',areas, 'TickDir','out')
title('VGAT-ChR2')
axis square

%% load and concate
Viral = load('Z:\All_Staff\home\miaomiao\labmeeting\HVAs\LED spread\191212-viral'); 
Vgat = load('Z:\All_Staff\home\miaomiao\labmeeting\HVAs\LED spread\191212-Vgat'); 

for i_area=1:4
    All{i_area} = [-Viral.C_change{i_area} -Vgat.C_change{i_area}]; 
    
end 
n_cell = cellfun(@length, All); 
figure
plotSpread(All, 'distributionColors', Viral.Colors);
hold on

errorbar(1:4, cellfun(@mean,All ), cellfun(@std,All)./sqrt(n_cell), 'LineStyle','none','Marker','o','Color',[0 0 0])
ylabel('Norm. suppression')
xlim([0.5 4.5])
set(gca,'XTicklabel',Vgat.areas, 'TickDir','out')
title('Combined')
axis square
%% 
num_cell = n_cell; 

data = [];
for i=1:4
    
    %data = [data; Adpt_idx_new{1,i}];
         data = [data; All{i}'];
    %data = [data;All.FT_adpt_idx{i,1}- All.FA_adpt_idx{i,1}];
end

group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3 ];
%group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3; ones(num_cell(5),1)+4; ones(num_cell(6),1)+5 ];

[p1,tbl,stats] = kruskalwallis(data,group);
b= multcompare(stats,'CType','dunn-sidak');
% [p,h] = ranksum(All{1},All{4} ); 


%% plot normalize suppression across areas
edges=Viral.edges; 

Viral = load('Z:\All_Staff\home\miaomiao\labmeeting\HVAs\LED spread\new_analysis\viral.mat');
Vgat = load('Z:\All_Staff\home\miaomiao\labmeeting\HVAs\LED spread\new_analysis\vgat.mat');
areas = {'V1','LM','AL','PM'};

for i_area = 1:4
    idx = [];
    idx = strcmp(Viral.C_area, areas{i_area});
    Viral.C_change{i_area} = Viral.C_pchg(idx);
    Viral.C_psth{i_area} = Viral.all.inhibitPSTH(idx,:);
    Viral.C_base{i_area} = Viral.all.inhibitbase(idx,:); 
    
    idx = [];
    idx = strcmp(Vgat.C_area, areas{i_area});
    Vgat.C_change{i_area} = Vgat.C_pchg(idx);
    Vgat.C_psth{i_area} = Vgat.all.inhibitPSTH(idx,:);
    Vgat.C_base{i_area} = Vgat.all.inhibitbase(idx,:); 
    
    
end
%% concatenate two approaches
for i_area = 1:4
    C_change{i_area} = [Viral.C_change{i_area} Vgat.C_change{i_area}]; 
    C_psth{i_area} = [Viral.C_psth{i_area}; Vgat.C_psth{i_area}];
    C_base{i_area} = [Viral.C_base{i_area}; Vgat.C_base{i_area}]; 
end 

Colors(1,:) = [0 0 0];
Colors(2:4,:) = lines(3); 

%%

LED_duration = Viral.spon.LED_duration;
figure
for i_area=1:4

subplot(2,2,i_area)
shadedErrorBar_ch(edges,  mean(C_psth{i_area}./repmat(C_base{i_area},1,size(C_psth{i_area},2)),1), std(C_psth{i_area}./repmat(C_base{i_area},1,size(C_psth{i_area},2)),[],1)./sqrt(size(C_psth{i_area},1)), {'Color',Colors(i_area,:),'markerfacecolor',Colors(i_area,:)})
hold on
xlim([-1.1 5])
ylim([0 2])
vline(0,'k:')
if i_area==4
    vline(Vgat.spon.LED_duration,'k:')
else
vline(LED_duration,'k:')
end 
set(gca,'TickDir','out')
% axis square
xlabel('Time (S)')
ylabel('Norm. to baseline')
title([areas{i_area} '  n= ' num2str(size(C_psth{i_area},1))])

end 