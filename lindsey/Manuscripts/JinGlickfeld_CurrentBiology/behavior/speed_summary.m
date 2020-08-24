homeDir = 'Z:\All_Staff\home\miaomiao\Analysis\Behavior\select and analysis';
xlsfile = fullfile(homeDir, 'speed.xlsx');
%[Gdata, Gtext, Graw] = xlsread(xlsfile,10);% contra
[Gdata, Gtext, Graw] = xlsread(xlsfile,14);%ipsi
IDs =Gdata(:,1);
Dir = Gtext(2:end,5);
region = Gtext(2:end,2);
power = Gdata(:,3);
background = Gtext(2:end,6);
Background = Gtext(2:end,8);
virus = Gtext(2:end,7);
Thresh = NaN(length(IDs),2); %  first collumn is the control, second is led
FA =  NaN(length(IDs),2);
FA_CI = NaN(length(IDs),2,2); % confidence interval for each mice
slope = NaN(length(IDs),2);
Lapse = NaN(length(IDs),2); % store lapse rate

dprime = NaN(length(IDs),5,2);
criterion = NaN(length(IDs),5,2);
speeds = NaN(length(IDs),6,2);
FA_RT = NaN(length(IDs),2);
FA_RT_sem = NaN(length(IDs),2);
RT = NaN(length(IDs),5,2); % store the reaction time
RT_sem = NaN(length(IDs),5,2); % store the sem 
LEDpre_thresh = NaN(length(IDs),4); % all combined: LED/Ctrl; LED/LED; Ctrl/Ctrl; Ctrl/LED
LEDpre_FA = NaN(length(IDs),4); % previouse led trial, followed by control/or followed by led

Pre_thresh = NaN(length(IDs),4,3);% mice * ctrl/ctrl combinations * pre-success, pre-early, pre-miss
Pre_FA = NaN(length(IDs),4,3);

Pre_fit = {};
TrialID = {};
FA_bintrial = NaN(length(IDs),4,2); % only has four bins
Thresh_bintrial = NaN(length(IDs),4,2);
bootStats = {}; % for each individual mice


fitSall = {};
Hit_2 = NaN(length(IDs),2);
d_2 = NaN(length(IDs),2);
c_2 = NaN(length(IDs),2);

for i_exp = 1:length(IDs)
    load(Dir{i_exp})
     % get the threshold for previouse LED, current trial is control (LED/Ctrl)
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.LED_pre.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.LED_pre.HT_num + Output{1}.Outcome{1}.target.LED_pre.Miss_num;
    LEDpre_fit{i_exp,1} = weibullFitLG(Output{1}.Outcome{1}.target.LED_pre.Speed, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,1) = LEDpre_fit{i_exp,1}.thresh;
    LEDpre_FA(i_exp,1) = Output{1}.Outcome{1}.FA.LED_pre.all.FA;
     
    
    % get the threshold for previouse LED, current trial is LED(LED/LED)
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.Same_pre.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.Same_pre.HT_num + Output{1}.Outcome{2}.target.Same_pre.Miss_num;
    LEDpre_fit{i_exp,2} = weibullFitLG(Output{1}.Outcome{2}.target.Same_pre.Speed, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,2) = LEDpre_fit{i_exp,2}.thresh;
    LEDpre_FA(i_exp,2) = Output{1}.Outcome{2}.FA.Same_pre.all.FA; 
    
    % get the threshold for previouse control, current trial is control (Ctrl/Ctrl)
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.Same_pre.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.Same_pre.HT_num + Output{1}.Outcome{1}.target.Same_pre.Miss_num;
    LEDpre_fit{i_exp,3} = weibullFitLG(Output{1}.Outcome{1}.target.Same_pre.Speed, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,3) = LEDpre_fit{i_exp,3}.thresh;
    LEDpre_FA(i_exp,3) = Output{1}.Outcome{1}.FA.Same_pre.all.FA;
     
    
    % get the threshold for previouse Control, current trial is LED(Ctrl/LED)
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.LED_pre.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.LED_pre.HT_num + Output{1}.Outcome{2}.target.LED_pre.Miss_num;
    LEDpre_fit{i_exp,4} = weibullFitLG(Output{1}.Outcome{2}.target.LED_pre.Speed, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,4) = LEDpre_fit{i_exp,4}.thresh;
    LEDpre_FA(i_exp,4) = Output{1}.Outcome{2}.FA.LED_pre.all.FA;
    
    
    
    
    for i = 1:2
        Hits = [];
        trialAll = [];
        Hits = Output{i}.target.c_hit;
        trialAll = Output{i}.target.HT_num + Output{i}.target.Miss_num;
        fitSall{i_exp,i} = weibullFitLG(Output{i}.Infor.speed, Hits',1, 1, {'nTrials', trialAll'});
        Thresh(i_exp,i) = fitSall{i_exp,i}.thresh;
        slope(i_exp,i) = fitSall{i_exp,i}.slope;
        % add bootstrap for each mouse
        bootStats{i_exp,i}=BootstrapWeibullFit(trialAll',Hits',1000,Output{i}.Infor.speed,1, 1);
        
        FA(i_exp,i) = Output{i}.FA.all.FA;
        FA_CI(i_exp,i,:) = Output{i}.FA.all.FA_confi; 
        
        dprime(i_exp,1:length(Output{i}.Infor.speed),i) = Output{i}.sdt.dprime;
        criterion(i_exp,1:length(Output{i}.Infor.speed),i) = Output{i}.sdt.criterion;
        speeds(i_exp,1:length(Output{i}.Infor.speed),i) = Output{i}.Infor.speed;
        RT(i_exp,1:length(Output{i}.Infor.speed),i) = Output{i}.target.RTonHit_mean(1,:);
        RT_sem(i_exp,1:length(Output{i}.Infor.speed),i) = Output{i}.target.RTonHit_mean(2,:);
        Lapse(i_exp,i) = 1-Hits(end);
        FA_RT(i_exp,i) = mean(Output{i}.FA.FA_RT); 
        FA_RT_sem(i_exp,i) = std(Output{i}.FA.FA_RT)./sqrt(length(Output{i}.FA.FA_RT)); 
        % interpolate the 2 deg/s for d' and c
        Hit_2(i_exp,i) = fitSall{i_exp,i}.modelFun(fitSall{i_exp,i}.coefEsts, 2); 
        [d_2(i_exp,i), c_2(i_exp,i)] = dprime_simple(Hit_2(i_exp,i),FA(i_exp,i));

        FA_bintrial(i_exp,1:length(Output{i}.bintrial.FA),i)  = Output{i}.bintrial.FA;
        Thresh_bintrial(i_exp,1:length( Output{i}.bintrial.thresh),i)  =  Output{i}.bintrial.thresh;
    end
end
cbin = Output{1}.bintrial.cbin.Hit; 
%% get the RT for easiest target or near threshold 
RT_easy = [];
RT_thresh = []; 
for i = 1:size(speeds,1)
% for easy trial, find the nearest value to the easiest control orien
temp1 = [];
temp1 = squeeze(speeds(i,:,1)); 
temp2 = [];
temp2 = squeeze(speeds(i,:,2)); 
[easy_con(i), idx1] = max(temp1); 
[c, idx2] = min(abs(temp2-easy_con(i))); 
RT_easy (i,1) = RT(i,idx1,1); 
RT_easy (i,2) = RT(i,idx2,2); 

[c, idx3(i)] = min(abs(temp1-Thresh(i,1))); 
thresh_con(i) = temp1(idx3(i)); 
[c, idx4(i)] = min(abs(temp2-thresh_con(i))); 
RT_thresh(i,1) = RT(i,idx3(i),1);
RT_thresh(i,2) = RT(i,idx4(i),2);

end 

%% only select mouse that has at least 20 trials for each bin/ori 

Mice_select = []; % for trial bin analysis
for i_exp = 1:length(IDs)
    load(Dir{i_exp})
    trial_c=min( Output{1}.bintrial.Trial_num); 
    trial_l=min( Output{2}.bintrial.Trial_num);
    if length(trial_c)>5
       Mice_select(i_exp) = trial_c(2)>=20; 
    else
       Mice_select(i_exp) = trial_c(1)>=20; 
    end 
   
end 
%% save individual plot 
ledcolor = {[0 0 0] [0.2 0.6 1] [1 0 0.6] [0 0.6 0.2]};
for i = 1:length(IDs)
    figure
    suptitle(sprintf('%d-%s-%3.2fmW-%s',IDs(i),region{i},power(i),Background{i}))  
     
    subplot(2,2,1)
    for i_led = 1:2
        
        Orien = speeds(i,:,i_led);
        Orien = Orien(~isnan(Orien)); 
        
        maxI = max(Orien);
        minI = min(Orien);
        xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
        h=line(xgrid, fitSall{i,i_led}.modelFun(fitSall{i,i_led}.coefEsts, xgrid), 'Color',ledcolor{i_led});
        hold on;
        plot(fitSall{i,i_led}.intensityX,fitSall{i,i_led}.fractCorrY, 'o','Color',ledcolor{i_led});
        %thresh = coefEsts(1)*[1 1];
        plot(fitSall{i,i_led}.thresh*[1 1], [0 fitSall{i,i_led}.threshY], '--','Color',ledcolor{i_led});
        plot(bootStats{i,i_led}.ci95, fitSall{i,i_led}.threshY*[1 1], 'Color',ledcolor{i_led});
       
        % set limits correctly
        xLim = [min(xgrid) max(xgrid)].* [0.75 1.25];
        xLim = 10.^ceil(log10(xLim) - [1 0]);    
        
    end
    axis([0.2 32 0 1]) 
    set(gca,'xscale','log','TickDir','out')
    xlabel('Speed Increment')
    ylabel('Hit Rate')
    axis square
    
    subplot(2,2,2)
    for i_led = 1:2
        
        scatter(0, FA(i,i_led),'MarkerEdgeColor',ledcolor{i_led})
        hold on
        plot([0 0], squeeze(FA_CI(i,i_led,:)), 'Color',ledcolor{i_led})
        
    end
    if max(FA(i,:))<0.2
        ylim([0 0.2])
    else
        ylim([0 ceil(max(FA(i,:))*10)./10])
    end
    set(gca,'XTick', [0], 'XTickLabel', {'Collapsed FA'},'TickDir','out')
    ylabel('FA rate')
    axis square
    
    subplot(2,2,3)
    
    for i_led = 1:2
        Orien = speeds(i,:,i_led);
        Orien = Orien(~isnan(Orien));
        errorbar(Orien, squeeze(RT(i,:,i_led)),squeeze(RT_sem(i,:,i_led)),'Color',ledcolor{i_led},'Marker','o')
        hold on
        errorbar(0,FA_RT(i,i_led),FA_RT_sem(i,i_led),'Color',ledcolor{i_led},'Marker','o')
    end
    
    
    box off
    xlim([-5 34])
    ylim([300 600])
    set(gca,'TickDir','out','XTick',0:10:30)
    xlabel('Speed Increment (deg/s)')
    ylabel('mean RT (ms)')
    axis square
    
    print(sprintf('%d-%s-%3.2fmW-%s.pdf',IDs(i),region{i},power(i),Background{i}),'-dpdf','-fillpage')
    close all
end
%% 




Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
all_region = {'V1','LM','AL','PM'};
figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Thresh(idx,:);
[h_t(i_region),p_t(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')

hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
 text(i_region,2,['n=' num2str(exp)])
ylim([-2 2])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
axis square

subplot(2,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA(idx,:);
[h_f(i_region),p_f(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-0.2 0.2])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
axis square

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = d_2(idx,:);
[h_d(i_region),p_d(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
%scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.7,['n=' num2str(exp)])
ylim([-1 1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta d prime-2 (LED - Ctrl)')
axis square

subplot(2,2,4)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = c_2(idx,:);
[h_c(i_region),p_c(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
%scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
%text(i_region,0.7,['n=' num2str(exp)])
ylim([-1 1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta c-2 (LED - Ctrl)')
axis square
%% plot only two region
h_t = [];
p_t = []; 
Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
all_region = {'AL','PM'};
figure
subplot(2,2,1)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = Thresh(idx,:);
[h_t(i_region),p_t(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region+2,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region+2,:),'Marker','o')

hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
 text(i_region,2,['n=' num2str(exp)])
ylim([-2 2])

end 

xlim([0 3])
hline(0,'k:')
set(gca,'XTick',1:1:2,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
axis square

subplot(2,2,2)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = FA(idx,:);
[h_f(i_region),p_f(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region+2,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region+2,:),'Marker','o')
hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-0.2 0.2])

end 

xlim([0 3])
hline(0,'k:')
set(gca,'XTick',1:1:2,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
axis square

subplot(2,2,3)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = d_2(idx,:);
[h_d(i_region),p_d(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region+2,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region+2,:),'Marker','o')
hold on
%scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.7,['n=' num2str(exp)])
ylim([-1 1])

end 

xlim([0 3])
hline(0,'k:')
set(gca,'XTick',1:1:2,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta d prime-2 (LED - Ctrl)')
axis square

subplot(2,2,4)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = c_2(idx,:);
[h_c(i_region),p_c(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region+2,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region+2,:),'Marker','o')
hold on
%scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
%text(i_region,0.7,['n=' num2str(exp)])
ylim([-1 1])

end 

xlim([0 3])
hline(0,'k:')
set(gca,'XTick',1:1:2,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta c-2 (LED - Ctrl)')
axis square


%% two way anova test
Ctrl = RT(idx,:,1);
LED = RT(idx,:,2);
[p,tbl,stats] = anova2([Ctrl;LED],3)
c= multcompare(stats,'Estimate','row');

%% summary of lapse rate, slope,
all_region = {'AL','PM'};
figure
subplot(1,2,1)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = Lapse(idx,:);
exp = sum(idx);
[h_l(i_region),p_l(i_region)]= ttest(temp(:,1),temp(:,2));
plot([i_region-0.3 i_region+0.3],temp, 'color',Colors(i_region+2,:))

hold on
errorbar([i_region-0.3 i_region+0.3],[mean(temp(:,1)) mean(temp(:,2))],[std(temp(:,1))./sqrt(exp) std(temp(:,2))./sqrt(exp)],'Color',Colors(i_region+2,:),'Marker','o','LineStyle','none')

end 

 ylim([0 0.5])
xlim([0.2 2.8])
set(gca,'TickDir','out','XTick',1:1:2, 'XTickLabel',all_region)%{'V1', 'LM', 'AL','PM'}
ylabel('lapse rate')
axis square

subplot(1,2,2)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = slope(idx,:);
exp = sum(idx);
[h_s(i_region),p_s(i_region)]= ttest(temp(:,1),temp(:,2));
plot([i_region-0.3 i_region+0.3],temp, 'color',Colors(i_region+2,:))

hold on
errorbar([i_region-0.3 i_region+0.3],[mean(temp(:,1)) mean(temp(:,2))],[std(temp(:,1))./sqrt(exp) std(temp(:,2))./sqrt(exp)],'Color',Colors(i_region+2,:),'Marker','o','LineStyle','none')

end 

 ylim([0 2])
xlim([0.2 2.8])
set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',all_region)%{'V1', 'LM', 'AL','PM'}
ylabel('Slope of fitted function')
axis square

%%  summary of reaction time for 90deg, reaction time near threshold orien for hits, reaction time for FAs 
new_c = [0.8 0.8 0.8; brighten(ledcolor{2},0.8)]; 

figure
% reaction time for FA rate
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = FA_RT(idx,:);
exp = sum(idx);
[h_frt(i_region),p_frt(i_region)]= ttest(temp(:,1),temp(:,2));
plot([i_region-0.3 i_region+0.3],temp, 'color',Colors(i_region+2,:))

hold on
errorbar([i_region-0.3 i_region+0.3],[mean(temp(:,1)) mean(temp(:,2))],[std(temp(:,1))./sqrt(exp) std(temp(:,2))./sqrt(exp)],'Color',Colors(i_region+2,:),'Marker','o','LineStyle','none')

end 

 ylim([300 600])
xlim([0.2 2.8])
set(gca,'TickDir','out','XTick',1:1:2, 'XTickLabel',all_region)
ylabel('FA-RT (ms)')
axis square

figure
subplot(2,2,1)
% Plot all RTs for AL
i_region = 1;
for i = 1:2
    idx = strcmp(region,all_region(i_region));
    s_temp = squeeze(speeds(idx,:,i));
    r_temp = squeeze(RT(idx,:,i));
    s_mean = mean(s_temp,1);
    s_sem = std(s_temp,[],1)./sqrt(sum(idx));
    r_mean = mean(r_temp,1);
    r_sem = std(r_temp,[],1)./sqrt(sum(idx));
    hold on
    for i_mouse = 1:sum(idx)
    plot(s_temp(i_mouse,:),r_temp(i_mouse,:),'color',new_c(i,:))
    hold on
    end 
    if i==1
        errorbarxy(s_mean, r_mean, s_sem, r_sem, {'ko-','k','k'})
    else
        hold on
        errorbarxy(s_mean, r_mean, s_sem, r_sem, {'bo-','b','b'})
    end
   
    
    
end
xlim([-5 34])
ylim([300 600])
set(gca,'TickDir','out','XTick',0:10:30)
xlabel('Speed Increment (deg/s)')
ylabel('AL mean RT (ms)')
axis square

subplot(2,2,2)
% normalized by the easiest target in control condition
i_region = 1;
for i = 1:2
    idx = strcmp(region,all_region(i_region));
    s_temp = squeeze(speeds(idx,:,i));
    r_temp = squeeze(RT(idx,:,i));
    r_temp = r_temp./squeeze(RT(idx,5,1)); 
    %r_temp = r_temp./(r_temp(:,5))
    s_mean = mean(s_temp,1);
    s_sem = std(s_temp,[],1)./sqrt(sum(idx));
    r_mean = mean(r_temp,1);
    r_sem = std(r_temp,[],1)./sqrt(sum(idx));
    hold on
    for i_mouse = 1:sum(idx)
    plot(s_temp(i_mouse,:),r_temp(i_mouse,:),'color',new_c(i,:))
    hold on
    end 
    if i==1
        errorbarxy(s_mean, r_mean, s_sem, r_sem, {'ko-','k','k'})
    else
        hold on
        errorbarxy(s_mean, r_mean, s_sem, r_sem, {'bo-','b','b'})
    end
   
    
    
end
xlim([-5 34])
ylim([0.9 1.5])
set(gca,'TickDir','out','XTick',0:10:30)
xlabel('Speed Increment (deg/s)')
ylabel('AL Norm. RT')
axis square

 subplot(2,2,3)
% Plot all RTs for AL
i_region = 2;
for i = 1:2
    idx = strcmp(region,all_region(i_region));
    s_temp = squeeze(speeds(idx,:,i));
    r_temp = squeeze(RT(idx,:,i));
    s_mean = mean(s_temp,1);
    s_sem = std(s_temp,[],1)./sqrt(sum(idx));
    r_mean = mean(r_temp,1);
    r_sem = std(r_temp,[],1)./sqrt(sum(idx));
    hold on
    for i_mouse = 1:sum(idx)
    plot(s_temp(i_mouse,:),r_temp(i_mouse,:),'color',new_c(i,:))
    hold on
    end 
    if i==1
        errorbarxy(s_mean, r_mean, s_sem, r_sem, {'ko-','k','k'})
    else
        hold on
        errorbarxy(s_mean, r_mean, s_sem, r_sem, {'bo-','b','b'})
    end
   
    
    
end
xlim([-5 34])
ylim([300 600])
set(gca,'TickDir','out','XTick',0:10:30)
xlabel('Speed Increment (deg/s)')
ylabel('PM mean RT (ms)')
axis square

subplot(2,2,4)
% normalized by the easiest target in control condition
i_region = 2;
for i = 1:2
    idx = strcmp(region,all_region(i_region));
    s_temp = squeeze(speeds(idx,:,i));
    r_temp = squeeze(RT(idx,:,i));
    r_temp = r_temp./squeeze(RT(idx,5,1));
    %r_temp = r_temp./(r_temp(:,5))
    s_mean = mean(s_temp,1);
    s_sem = std(s_temp,[],1)./sqrt(sum(idx));
    r_mean = mean(r_temp,1);
    r_sem = std(r_temp,[],1)./sqrt(sum(idx));
    hold on
    for i_mouse = 1:sum(idx)
    plot(s_temp(i_mouse,:),r_temp(i_mouse,:),'color',new_c(i,:))
    hold on
    end 
    if i==1
        errorbarxy(s_mean, r_mean, s_sem, r_sem, {'ko-','k','k'})
    else
        hold on
        errorbarxy(s_mean, r_mean, s_sem, r_sem, {'bo-','b','b'})
    end
   
    
    
end
xlim([-5 34])
ylim([0.9 1.5])
set(gca,'TickDir','out','XTick',0:10:30)
xlabel('Speed Increment (deg/s)')
ylabel('PM Norm. RT')
axis square
 

%% stats
i_region = 2;
% idx = strcmp(region,all_region(i_region));
% temp = RT(idx,:,:);
% data = [temp(:,1:2);temp(:,3:4)]; 
idx = strcmp(region,all_region(i_region))& Mice_select';
temp = squeeze(FA_bintrial(idx,:,:)); 
exp = sum(idx);
%data = [squeeze(temp(:,:,1))./squeeze(temp(:,1,1)); squeeze(temp(:,:,2))./squeeze(temp(:,1,2))]; 
data = [squeeze(temp(:,:,1)); squeeze(temp(:,:,2))]; 
[p,tbl,stats] = anova2(data,exp);

%% Only for mice that has enough trials 

figure
for i_region = 1:2
    
    subplot(1,2,i_region)
    idx = strcmp(region,all_region(i_region))& Mice_select';
    temp = squeeze(Thresh_bintrial(idx,:,1));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(Thresh_bintrial(idx,:,2));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(1.1,3.8,[all_region(i_region) num2str(sum(idx))],'color',Colors(i_region,:))
    
    ylim([0 4])
    xlim([0.7 4.2])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Thresh')
    axis square
end

figure
for i_region = 1:2
    
    subplot(1,2,i_region)
    idx = strcmp(region,all_region(i_region))& Mice_select';
    temp = squeeze(FA_bintrial(idx,:,1));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(FA_bintrial(idx,:,2));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(1.1,0.19,[all_region(i_region) num2str(sum(idx))],'color',Colors(i_region,:))
    ylim([0 0.3])
    xlim([0.7 4.2])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('FA rate')
    axis square
end

 %% for selected mice norm to the first 
figure
for i_region = 1:2
    
    subplot(1,2,i_region)
    idx = strcmp(region,all_region(i_region))& Mice_select';
    temp = squeeze(Thresh_bintrial(idx,:,1));
    temp = temp./repmat(temp(:,1),1,4);
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(Thresh_bintrial(idx,:,2));
    temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(1,1.5,[all_region(i_region) num2str(sum(idx))],'color',Colors(i_region,:))
    ylim([0.6 2])
    xlim([0.7 4.2])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Norm. Thresh')
    axis square
end

figure
for i_region = 1:2
    
    subplot(1,2,i_region)
    idx = strcmp(region,all_region(i_region))& Mice_select';
    temp = squeeze(FA_bintrial(idx,:,1));
    temp = temp./repmat(temp(:,1),1,4);
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(FA_bintrial(idx,:,2));
    temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(1,0.6,[all_region(i_region) num2str(sum(idx))],'color',Colors(i_region,:))
    ylim([0.4 1.2])
    xlim([0.7 4.2])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Norm. FA rate')
    axis square
end

%% stats
i_region = 2;
% idx = strcmp(region,all_region(i_region));
% temp = Pretrl_FA(idx,:);
% data = [temp(:,1:2);temp(:,3:4)]; 
idx = strcmp(region,all_region(i_region))& Mice_select';
temp = squeeze(Thresh_bintrial(idx,:,:)); 
exp = sum(idx);
%data = [squeeze(temp(:,:,1))./squeeze(temp(:,1,1)); squeeze(temp(:,:,2))./squeeze(temp(:,1,2))]; 
data = [squeeze(temp(:,:,1)); squeeze(temp(:,:,2))]; 
data(:,2:3) =[];

[p,tbl,stats] = anova2(data,exp);
%c = multcompare(stats,'Estimate','Column');

%% plot threshold and FA rate seperated by preceding trials 
% reorganize it so that Ctrl/Ctrl LED/Ctrl Ctrl/LED LED/LED
Pretrl_FA(:,1) = LEDpre_FA(:,3);
Pretrl_FA(:,2) = LEDpre_FA(:,1);
Pretrl_FA(:,3) = LEDpre_FA(:,4);
Pretrl_FA(:,4) = LEDpre_FA(:,2);

Pretrl_thresh(:,1) = LEDpre_thresh(:,3);
Pretrl_thresh(:,2) = LEDpre_thresh(:,1);
Pretrl_thresh(:,3) = LEDpre_thresh(:,4);
Pretrl_thresh(:,4) = LEDpre_thresh(:,2);

%% thresh
figure
for i_region = 1:2
    subplot(1,2,i_region)
    
    idx = strcmp(region,all_region(i_region));
    temp = Pretrl_thresh(idx,:);
    exp = sum(idx);
    for i=1:exp
    plot(1:4,temp(i,:), 'color',Colors(i_region,:))
    
    hold on
    end
    errorbar(1:4,mean(temp,1),std(temp,[],1)./sqrt(exp),'Color',Colors(i_region,:),'Marker','o','LineStyle','none')
    text(0.7,3.8,all_region(i_region),'color',Colors(i_region,:))
    ylim([0 4])
    xlim([0.5 4.5])
    set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',{'Ctrl|Ctrl','LED/Ctrl','Ctrl/LED','LED/LED'})
    ylabel('Threshold (deg)')
    axis square
    xtickangle(30)
end

%% FA
figure
for i_region = 1:2
    subplot(1,2,i_region)
    
    idx = strcmp(region,all_region(i_region));
    temp = Pretrl_FA(idx,:);
    exp = sum(idx);
    for i=1:exp
    plot(1:4,temp(i,:), 'color',Colors(i_region,:))
    
    hold on
    end 
    errorbar(1:4,mean(temp,1),std(temp,[],1)./sqrt(exp),'Color',Colors(i_region,:),'Marker','o','LineStyle','none')
    text(0.7,0.28,all_region(i_region),'color',Colors(i_region,:))
     ylim([0 0.3])
    xlim([0.5 4.5])
    set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',{'Ctrl|Ctrl','LED/Ctrl','Ctrl/LED','LED/LED'})
    ylabel('FA rate')
    axis square
    xtickangle(30)
end

%% get the mean threshold and FA rate from unique mice

uniqueID = unique(IDs(strcmp(region,'PM')|strcmp(region,'AL'))); 
for i=1:length(uniqueID)
    idx=[];
    idx = IDs==uniqueID(i); 
    ThreshU(i) = mean(Thresh(idx,1)); 
    FAU(i) = mean(FA(idx,1));
    lapseU(i) = mean(Lapse(idx,1));
    d22U(i)=mean(d_2(idx,1)); 
end 
mean(ThreshU)
std(ThreshU)./sqrt(length(ThreshU))
mean(FAU)
std(FAU)./sqrt(length(FAU))

mean(d22U)
std(d22U)./sqrt(length(d22U))


%% stats chose the unique mice
mouse = IDs(logical(Mice_select)); 
data = [];
data = squeeze(FA_bintrial(logical(Mice_select),:,1)); 
Uni_select = logical([1 0 1 1 0 0 1 1 0]); 
%data =  data(Uni_select,:);
data = data(Uni_select,:)./(data(Uni_select,1)); 
[p,tbl,stats] = anova1(data); 


%% plot within each area, LED effects seperate by suppression method
back_color = [0 0 0; 0.5 0.5 0.5]; %[1 0.6 0.6;0 0 0; 0.5 0.5 0.5; 0.6 0.8 0.8];  % black PV-ChR2, gray PV-Chronos, red:GAD-ChR2; green:vgat-ChR2
B_unique = unique(Background); 

figure
for i_region = 1:2
    subplot(1,2,i_region)
    idx = strcmp(region,all_region(i_region));
    for i_b =1:size(B_unique,1)
        idx2 = strcmp(Background,B_unique(i_b));
        idx3 = idx&idx2;
        temp = Thresh(idx3,:);
        if ~isempty(temp)
            scatter(temp(:,1),temp(:,2),'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)
            hold on
        end
    end
    text(0.2,2.8,all_region(i_region),'color',Colors(i_region,:))
    ylim([0 3])
    xlim([0 3])
    hold on
    plot([0 3],[0 3],'k:')
    set(gca,'TickDir','out','XTick',0:1:3, 'YTick',0:1:3)
    ylabel('Thresh-LED')
    xlabel('Thresh-Ctrl')
    axis square
    
    
end

figure

for i_region = 1:2
    subplot(1,2,i_region)
    idx = strcmp(region,all_region(i_region));
    for i_b =1:size(B_unique,1)
        idx2 = strcmp(Background,B_unique(i_b));
        idx3 = idx&idx2;
        temp = FA(idx3,:);
        if ~isempty(temp)
            scatter(temp(:,1),temp(:,2),'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)
            hold on
        end
    end
    text(0.02,0.28,all_region(i_region),'color',Colors(i_region,:))
    ylim([0 0.3])
    xlim([0 0.3])
    hold on
    plot([0 0.3],[0 0.3],'k:')
    set(gca,'TickDir','out','XTick',0:0.1:0.3, 'YTick',0:0.1:0.3)
    ylabel('FA-LED')
    xlabel('FA-Ctrl')
    axis square
    
    
end

%% plot delta change over light power

figure
subplot(2,2,1)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',Colors(i_region+2,:),'MarkerFaceColor',Colors(i_region+2,:),'SizeData',40)

hold on


end 
ylim([-2 2])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
xlabel('light power (mW)')
axis square

subplot(2,2,2)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',Colors(i_region+2,:),'MarkerFaceColor',Colors(i_region+2,:),'SizeData',40)

hold on


end 
ylim([-0.2 0.2])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
xlabel('light power (mW)')
axis square
% lapse rate difference
subplot(2,2,3)
for i_region = 1:2
idx = strcmp(region,all_region(i_region));
temp = Lapse(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',Colors(i_region+2,:),'MarkerFaceColor',Colors(i_region+2,:),'SizeData',40)

hold on


end 
ylim([-0.1 0.1])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta Lapse (LED - Ctrl)')
xlabel('light power (mW)')
axis square
%% plot delta change over light power seperate by methods 
figure
subplot(2,2,1)
for i_b = 1:2
idx = strcmp(Background,B_unique(i_b)) & (strcmp(region,'AL')|strcmp(region,'PM'));
% not include other than AL and PM

temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)

hold on


end 
ylim([-2 2])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
xlabel('light power (mW)')
axis square
subplot(2,2,2)
for i_b = 1:2
idx = strcmp(Background,B_unique(i_b)) & (strcmp(region,'AL')|strcmp(region,'PM'));
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)

hold on


end 
ylim([-0.2 0.2])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
xlabel('light power (mW)')
axis square
% lapse rate difference
subplot(2,2,3)
for i_b = 1:2
idx = strcmp(Background,B_unique(i_b)) & (strcmp(region,'AL')|strcmp(region,'PM'));
temp = Lapse(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)

hold on


end 
ylim([-0.1 0.1])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta Lapse (LED - Ctrl)')
xlabel('light power (mW)')
axis square