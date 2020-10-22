input = inputcombined;
Success = strcmp(input.trialOutcomeCell,'success');
trialnum = length(Success);
ReactTime = double(cell2mat(input.reactTimesMs));
Early = strcmp(input.trialOutcomeCell,'failure');
Miss = strcmp(input.trialOutcomeCell,'ignore');
HoldTimeMs = double (cell2mat(input.holdTimesMs));
reqHoldMs =  double (cell2mat(input.tTotalReqHoldTimeMs));
raw_RT = ReactTime; 
tspeed = double(cell2mat(input.tDotSpeedDPS)) - double(cell2mat(input.tBaseDotSpeedDPS));
basespeed = unique(double(cell2mat(input.tBaseDotSpeedDPS))); 
if length(unique(tspeed))>7
    % manually adjust for each mouse
    tspeed(tspeed>0.4 & tspeed <1) = 0.9;
    tspeed(tspeed>1 & tspeed <2) = 1.8;
    
    tspeed(tspeed>2 & tspeed <4) = 3.5;
    
    tspeed(tspeed>4 & tspeed <8) = 7.2;
    
      
end

speed = unique(tspeed);

idx_HM = find(Early==0); % index of Hit and Miss trials
% just use ReactTime in HDC condtions as new RT


RT_HM = zeros(1,length(ReactTime));
RT_HM(idx_HM) = ReactTime(idx_HM);

% adjust the Success, Miss trial conditions
% throw too fast response into early trials
SuccessN = zeros(1,length(Success));
MissN = zeros(1,length(Miss));
EarlyN = Early;
% reaction window is 200-700 ms
for i=1:length(idx_HM)
    idx_temp = [];
    idx_temp = idx_HM(i);% trial of hit or miss trial
    if RT_HM(idx_temp)>=200 && RT_HM(idx_temp)<=700
        SuccessN(idx_temp) = 1;   % new success trials
    end
    if RT_HM(idx_temp)<200
        EarlyN(idx_temp) = 1;   % add too fast time into early trials
    end
    if RT_HM(idx_temp)>700
        MissN(idx_temp) = 1;   % new miss trials
    end
end
% get the reaction time distribution after the stimulus onset regardless of
% too fast or misses
for i_speed = 1: length(speed)
    idx = [];
    idx = tspeed==speed(i_speed);
    Ouput.rawRT{i_speed,1} = raw_RT(idx); 
end 
% get the reaction time distribution on window of 200-700 ms  
[Output.target] = Contrast_HR(speed,SuccessN,MissN,tspeed,RT_HM);
% get the FA rate
[Output.FA] = spd_FA(EarlyN,trialnum,HoldTimeMs,600); % 600 is the boundaries where the earliest change detect would occur
% get the response for short vs long trials 
% artificially define short: <=2s; long: >=3s
short_idx = reqHoldMs<=2000; 
long_idx = reqHoldMs>=3000; 
[Output.Tshort] = Contrast_HR(speed,SuccessN(short_idx),MissN(short_idx),tspeed(short_idx),RT_HM(short_idx));
[Output.Tlong] = Contrast_HR(speed,SuccessN(long_idx),MissN(long_idx),tspeed(long_idx),RT_HM(long_idx));


% get the boostrap 

fitSall = {};
bootStatsall={};
NBootstrapReps = 1000;

Hits = Output.target.c_hit;
trialAll = Output.target.HT_num + Output.target.Miss_num;
fitSall = weibullFitLG(speed, Hits',1, 1, {'nTrials', trialAll'});

[bootStatsall]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,speed,1, 1);
Output.target.thresh = fitSall.thresh;
Output.Infor.speed = speed;
Output.Infor.ID = input.ID; 

% fit_short = {};
% bootStats_short={};
% 
% Hits = Output.Tshort.c_hit;
% trialAll = Output.Tshort.HT_num + Output.Tshort.Miss_num;
% fit_short = weibullFitLG(speed, Hits',1, 1, {'nTrials', trialAll'});
% 
% [bootStats_short]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,speed,1, 1);
% fit_long = {};
% bootStats_long={};
% 
% Hits = Output.Tlong.c_hit;
% trialAll = Output.Tlong.HT_num + Output.Tlong.Miss_num;
% fit_long = weibullFitLG(speed, Hits',1, 1, {'nTrials', trialAll'});
% 
% [bootStats_long]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,speed,1, 1);

%get the fit of the model for RT distributions
Color = flipud(copper(5));
RT_bins = 200:50:701;
N_short = [];
N_long = [];
N = [];
for i_speed = 1:length(speed)
     N_short(i_speed,:) = (hist(Output.Tshort.RTonHit{i_speed},RT_bins)./length(Output.Tshort.RTonHit{i_speed}))./50; % transfer to density
     N_long(i_speed,:) = (hist(Output.Tlong.RTonHit{i_speed},RT_bins)./length(Output.Tlong.RTonHit{i_speed}))./50;
     N(i_speed,:) = (hist(Output.target.RTonHit{i_speed},RT_bins)./length(Output.target.RTonHit{i_speed}))./50; % transfer to density
end 
% [beta_short,R_short,MSE_short]=eLATER(N_short,Color);
% 
% [beta_long,R_long,MSE_long]=eLATER(N_long,Color);
[beta,R,MSE]=eLATER(N,Color);
% % fit the hyperbolar with reaction time 
% x=speed;
% y=Output.target.RTonHit_mean(1,:)./1000;
% hyprb = @(b,x) b(1) + b(2)./x;                 % Generalised Hyperbola
% NRCF = @(b) norm(y - hyprb(b,x));              % Residual Norm Cost Function
% B0 = [1; 1];
% B = fminsearch(NRCF, B0);                       % Estimate Parameters

save(['i' num2str(input.ID) '_' 'speed.mat'])
%% plot the difference between long and short trials
figure
subplot(2,2,1)
maxI = max(speed);
minI = min(speed);
xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
h=line(xgrid, fit_short.modelFun(fit_short.coefEsts, xgrid), 'Color',[0.5 0.5 0.5]);
hold on;
plot(fit_short.intensityX,fit_short.fractCorrY, 'o','Color',[0.5 0.5 0.5]);
%thresh = coefEsts(1)*[1 1];
plot(fit_short.thresh*[1 1], [0 fit_short.threshY], '--','Color',[0.5 0.5 0.5]);
plot(bootStats_short.ci95, fit_short.threshY*[1 1], 'Color',[0.5 0.5 0.5]);

h=line(xgrid, fit_long.modelFun(fit_long.coefEsts, xgrid), 'Color',[0 0 0]);
hold on;
plot(fit_long.intensityX,fit_long.fractCorrY, 'o','Color',[0 0 0]);
%thresh = coefEsts(1)*[1 1];
plot(fit_long.thresh*[1 1], [0 fit_long.threshY], '--','Color',[0 0 0]);
plot(bootStats_long.ci95, fit_long.threshY*[1 1], 'Color',[0 0 0]);
% set limits correctly
ylim([0 1])
xlim([0 40])

set(gca,'xscale','log','TickDir','out')

xlabel('speed change (DPS)')
ylabel('Hit Rate')
axis square
title('Black: long; Gray: short')

subplot(2,2,2)
scatter(0,Output.FA.S_FA,'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
plot([0 0],Output.FA.S_FA_confi, 'Color',[0.5 0.5 0.5])
scatter(0,Output.FA.L_FA,'MarkerEdgeColor',[0 0 0])

plot([0 0],Output.FA.L_FA_confi, 'Color',[0 0 0])
xlim([-0.5 0.5])
ylim([0 0.2])
set(gca,'TickDir','out')
xlabel('speed change (DPS)')
ylabel('FA rate')
axis square

subplot(2,2,3)
h=errorbar(speed,Output.Tshort.RTonHit_mean(1,:),Output.Tshort.RTonHit_mean(2,:));
h.Marker = 'o';
h.Color = [0.5 0.5 0.5];
hold on
h=errorbar(speed,Output.Tlong.RTonHit_mean(1,:),Output.Tlong.RTonHit_mean(2,:));
h.Marker = 'o';
h.Color = [0 0 0];

set(gca,'xscale','log','TickDir','out')
xlim([xgrid(1) 40])
ylim([300 600])
xlabel('speed change (DPS)')
ylabel('mean RT (ms)')
axis square

subplot(2,2,4)
scatter(speed,beta_short(1:5)./beta_short(6:10),'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
scatter(speed,beta_long(1:5)./beta_long(6:10),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
box off;
xlim([0 40])
ylim([0 15])
set(gca,'TickDir','out')
ylabel('Drift rate: mean/sd')
xlabel('speed change (DPS)')
axis square

text(10,2,['Short:' sprintf('%.2f',beta_short(end))],'color',[0.5 0.5 0.5])
text(10,1,['Long:' sprintf('%.2f',beta_long(end))],'color',[0 0 0])



%% plot the speed change detection
figure

ou = Output.target;

supertitle(['Base speed:' num2str(basespeed) '  Cohrence: 0%'])
subplot(2,2,1)
maxI = max(speed);
minI = min(speed);
xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
h=line(xgrid, fitSall.modelFun(fitSall.coefEsts, xgrid), 'Color',[0 0 0]);
hold on;
plot(fitSall.intensityX,fitSall.fractCorrY, 'o','Color',[0 0 0]);
%thresh = coefEsts(1)*[1 1];
plot(fitSall.thresh*[1 1], [0 fitSall.threshY], '--','Color',[0 0 0]);
plot(bootStatsall.ci95, fitSall.threshY*[1 1], 'Color',[0 0 0]);
title(['Thresh: ' num2str(fitSall.thresh)])
% set limits correctly

axis([min(xgrid) max(xgrid) 0 1])
set(gca,'xscale','log','XTick',[0.1 1 10],'TickDir','out')

xlabel('speed change (DPS)')
ylabel('Hit Rate')
axis square

subplot(2,2,2)
for i_led = 1:2
    
   
    scatter(0, Output.FA.all.FA,'MarkerEdgeColor','k')
    hold on
    plot([0 0], Output.FA.all.FA_confi, 'k')
    
end
ylim([0 0.2])

set(gca,'XTick', 0,'TickDir','out')
axis square
ylabel('FA')

subplot(2,2,3)
colors =[0.8 0.8 0.8; 0.6 0.6 0.6; 0.4 0.4 0.4; 0.2 0.2 0.2;0 0 0]; 
for i_speed =1:length(speed)
    
    h(i_speed) =cdfplot(Ouput.rawRT{i_speed,1}(Ouput.rawRT{i_speed,1}>=0&Ouput.rawRT{i_speed,1}<=700));
    %h(i_speed) =cdfplot(ou.RTonHit{i_speed});
    h(i_speed).Color = colors(i_speed,:); 
    hold on
    grid off
end 
ylim([0 1])
xlim([0 800])
vline([200 700],'k:')

set(gca,'TickDir','out')
xlabel('Reaction time (ms)')
axis square

subplot(2,2,4)
errorbar(speed,ou.RTonHit_mean(1,:),ou.RTonHit_mean(2,:),'ko');
hold on
% xplot = 0:1:40;
% plot(xplot,hyprb(B,xplot)*1000,'k')
xlim([0 40])
ylim([100 700])
set(gca,'TickDir','out')
xlabel('speed change (DPS)')
ylabel('mean RT (ms)')
axis square


print(['W:\All_Staff\home\miaomiao\Analysis\Behavior\select and analysis\speed\' num2str(input.ID)],'-dpdf','-fillpage')
%% plot the reaction time first 
Color = flipud(copper(5));
RT_bins = 200:50:701;
N = [];
figure
for i_speed = 1:length(speed)
    N(i_speed,:) = (hist(ou.RTonHit{i_speed},RT_bins)./length(ou.RTonHit{i_speed}))./50; % transfer to density
    x=[];
    x = sort(ou.RTonHit{i_speed});
    n=numel(ou.RTonHit{i_speed});
    y=(1:n)./n;
    p  = [1 2 5 10:10:90 95 98 99]/100; % some arbitrary probabilities
    q=quantile(ou.RTonHit{i_speed},p); 
    
    subplot(1,2,1)
    plot(RT_bins, N(i_speed,:),'color',Color(i_speed,:),'LineWidth',1.5);
    hold on
    xlim([0 800])
    xlabel('Reaction time (ms)')
    ylabel('Probability')
    box off;
    set(gca,'TickDir','out')
    axis square
    
    subplot(1,2,2)
    plot(x,y,'color',Color(i_speed,:),'LineWidth',1.5)
    hold on
    scatter(q,p,'MarkerFaceColor',Color(i_speed,:), 'MarkerEdgeColor',Color(i_speed,:))
    xlim([0 800])
    xlabel('Reaction Time (ms)')
    ylabel('Cumulative probability')
    box off;
    set(gca,'TickDir','out')
    axis square
    
end
%% fit the RT pdf with extended LATER model and compare with normal LATER model

% [beta1,R1,MSE1]=LATER(N,Color);
%% plot the extended later model
figure
subplot(2,2,1)
t1 = 1:10:800;
for i = 1:5
    scatter(RT_bins,N(i,:),'SizeData',20,'MarkerFaceColor',Color(i,:),'MarkerEdgeColor',Color(i,:))
    hold on
line(t1, ((t1 + (beta(i).*((beta(11)./beta(5+i)).^2)))./(((t1.^2)+ ((beta(11)./beta(5+i)).^2)).^1.5)).*(1./(sqrt(2.*pi).*beta(5+i))).*exp(-(1./(2.*(beta(5+i).^2))).*((1./t1 - beta(i)).^2).*((t1.^2)./(t1.^2 +(beta(11)./beta(5+i)).^2))),'color',Color(i,:),'linewidth',1)
end
set(gca,'TickDir','out')
xlim([0 800])
ylim([0 0.009])
xlabel('Reaction Time (ms)')
ylabel('Probability')
title(['MSE= ' num2str(MSE) 'evid= ' num2str(beta(end))])
axis square
subplot(2,2,2)
% mean drift rate over mean evidence
gscatter(speed,beta(1:5),[1 2 3 4 5],Color)
box off;
legend off;
xlim([0 40])
ylim([0 0.003])
set(gca,'TickDir','out')
ylabel('mean rate/mean bounds')
axis square

subplot(2,2,3)
gscatter(speed,beta(6:10),[1 2 3 4 5],Color)
box off;
legend off;
xlim([0 40])
ylim([0 0.0006])
set(gca,'TickDir','out')
ylabel('SD rate/mean bounds')
axis square
subplot(2,2,4)
gscatter(speed,beta(1:5)./beta(6:10),[1 2 3 4 5],Color)
box off;
legend off;
xlim([0 40])
ylim([0 20])
set(gca,'TickDir','out')
ylabel('Drift rate: mean/sd')
axis square
%% plot the Later model
figure
subplot(2,2,1)
t1 = 1:10:800;
for i = 1:5
    scatter(RT_bins,N(i,:),'SizeData',20,'MarkerFaceColor',Color(i,:),'MarkerEdgeColor',Color(i,:))
    hold on
line(t1, (1./(t1.^2)).*(1./(sqrt(2.*pi).*beta1(5+i))).*exp(-(1./(2.*(beta1(5+i).^2))).*(((1./t1) - beta1(i)).^2)),'color',Color(i,:),'linewidth',1)
end
set(gca,'TickDir','out')
xlim([0 800])
ylim([0 0.009])
xlabel('Reaction Time (ms)')
ylabel('Probability')
title(['MSE= ' num2str(MSE1)])
axis square
subplot(2,2,2)
% mean drift rate over mean evidence
gscatter(speed,beta1(1:5),[1 2 3 4 5],Color)
box off;
legend off;
xlim([0 40])
ylim([0 0.003])
set(gca,'TickDir','out')
ylabel('mean rate/mean bounds')
axis square

subplot(2,2,3)
gscatter(speed,beta1(6:10),[1 2 3 4 5],Color)
box off;
legend off;
xlim([0 40])
ylim([0 0.0006])
set(gca,'TickDir','out')
ylabel('SD rate/mean bounds')
axis square
subplot(2,2,4)
gscatter(speed,beta(1:5)./beta(6:10),[1 2 3 4 5],Color)
box off;
legend off;
xlim([0 40])
ylim([0 10])
set(gca,'TickDir','out')
ylabel('Drift rate: mean/sd')
axis square

%% plot the inverse RT and cumulative probability

inv_bins = 0:0.2:8; 
N_inv = [];
figure
for i_speed = 1:length(speed)
    RT_inv = [];
    RT_inv = (1./(ou.RTonHit{i_speed}))*1000;
    N_inv(i_speed,:) = hist(RT_inv,inv_bins)./(length(ou.RTonHit{i_speed}));
    x=[];
    x = sort(RT_inv);
    n=numel(RT_inv);
    y=(1:n)./n;
    p  = [1 2 5 10:10:90 95 98 99]/100; % some arbitrary probabilities
    q=quantile(RT_inv,p); 
  
    
    subplot(1,2,1)
    plot(inv_bins, N_inv(i_speed,:),'color',Color(i_speed,:),'LineWidth',1.5);
    hold on
    xlim([0 8])
    xlabel('Promptness (s^-^1)')
    ylabel('Probability')
    box off;
    set(gca,'TickDir','out')
    axis square
    
   
    subplot(1,2,2)
    plot(x,y,'color',Color(i_speed,:),'LineWidth',1.5)
    hold on
    scatter(q,p,'MarkerFaceColor',Color(i_speed,:), 'MarkerEdgeColor',Color(i_speed,:))
    xlim([0 8])
    xlabel('Promptness (s^-^1)')
    ylabel('Cumulative probability')
    box off;
    set(gca,'TickDir','out')
    axis square
    
    
  
end

%% pribit scale
figure
for i_speed = 1:length(speed)
    rt = ou.RTonHit{i_speed};
    rtinv = [];
    rtinv = (1./(rt));
    x= -1./sort(rt); % multiply by -1 to mirror abscissa
    n = numel(rtinv);
    y= pa_probit((1:n)./n);
    % quantiles
    p    = [1 2 5 10:10:90 95 98 99]/100;
    probit  = pa_probit(p);
    q    = quantile(rt,p);
    q    = -1./q;
    xtick  = sort(-1./(150+[0 pa_oct2bw(50,-1:5)])); % some arbitrary xticks
    
    plot(x,y,'color',Color(i_speed,:),'LineWidth',1.5)
    hold on
    scatter(q,probit,'MarkerFaceColor',Color(i_speed,:), 'MarkerEdgeColor',Color(i_speed,:));
 
    set(gca,'XTick',xtick,'XTickLabel',-1./xtick);
    xlim([min(xtick) max(xtick)]);
    set(gca,'YTick',probit,'YTickLabel',p*100);
    ylim([pa_probit(0.1/100) pa_probit(99.9/100)]);
    axis square;
    box off
    xlabel('Reaction time (ms)');
    ylabel('Cumulative probability');
    title('Reciprobit plot');
    
   if i_speed ==5
    x = q(1:3);
    y = probit(1:3);
    b = regstats(y,x);
    h = pa_regline(b.beta,'r-.');
    x = q(4:15);
    y = probit(4:15);
    b = regstats(y,x);
    h = pa_regline(b.beta,'r-.');
    
   end 
     
    
end 
%% summary the threshold, FA rate and reaction time across mice

Allmice = {'i418','i423','i595','i747'};
FA = zeros(size(Allmice,2),1);
Thresh = zeros(size(Allmice,2),1);
Speed = zeros(size(Allmice,2),5);
RT = zeros(size(Allmice,2),5);


for i_mice = 1: size(Allmice,2)
    name = [Allmice{i_mice} '_' 'speed.mat'];
    load(name)
    FA(i_mice) = Output.FA.all.FA;
    Thresh(i_mice) = Output.target.thresh; 
    Speed(i_mice,:) = Output.Infor.speed;
    RT(i_mice,:) = Output.target.RTonHit_mean(1,:); 
   
end 

%% plot the individual one and summayr
figure
subplot(2,2,1)
for i_mice = 1:size(Allmice,2)
    scatter(1,Thresh(i_mice),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',10)
    hold on
    
end
errorbar(1,mean(Thresh),std(Thresh)./sqrt(size(Allmice,2)),'Color',[0 0 0],'Marker','o')
ylim([0 4])
ylabel('Thresh (deg/s)')
set(gca,'TickDir','out')
axis square

subplot(2,2,2)
for i_mice = 1:size(Allmice,2)
    scatter(1,FA(i_mice),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',10)
    hold on
   
end
 errorbar(1,mean(FA),std(FA)./sqrt(size(Allmice,2)),'Color',[0 0 0],'Marker','o')
ylim([0 0.2])
ylabel('FA rate')
set(gca,'TickDir','out')
axis square

subplot(2,2,3)

for i_mice = 1:size(Allmice,2)
    plot(Speed(i_mice,:),RT(i_mice,:),'color',[0.5 0.5 0.5])
     hold on
    if i_mice ==3
         plot(Speed(i_mice,:),RT(i_mice,:),'color',[1 0 0])
    end 
    
   
   
end
 errorbar(mean(Speed,1),mean(RT,1),std(RT,[],1)./sqrt(size(Allmice,2)),'Color',[0 0 0],'Marker','o')
ylim([300 600])
xlim([-2 34])
ylabel('Mean RT')
set(gca,'TickDir','out','XTick',0:10:30)
axis square

%% overlay different baseline speed 
All{1}=load('W:\All_Staff\home\miaomiao\Analysis\Behavior\select and analysis\speed\i747_speed.mat');
All{2} = load('W:\All_Staff\home\miaomiao\Analysis\Behavior\select and analysis\speed\base10\i747_speed.mat');
%%
colors = [0 0 0; 1 0 0]; 
figure
supertitle('Base speed:0.5 dps (black) vs 10 dps(red)')
for i=1:2
ou = All{i}.Output.target;

subplot(2,2,1)
maxI = max(All{i}.speed);
minI = min(All{i}.speed);
xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
h=line(xgrid, All{i}.fitSall.modelFun(All{i}.fitSall.coefEsts, xgrid), 'Color',colors(i,:));
hold on;
plot(All{i}.fitSall.intensityX,All{i}.fitSall.fractCorrY, 'o','Color',colors(i,:));
%thresh = coefEsts(1)*[1 1];
plot(All{i}.fitSall.thresh*[1 1], [0 All{i}.fitSall.threshY], '--','Color',colors(i,:));
plot(All{i}.bootStatsall.ci95, All{i}.fitSall.threshY*[1 1], 'Color',colors(i,:));

text(0.2,1-0.1*i,['Thresh: ' num2str(All{i}.fitSall.thresh)],'Color',colors(i,:))
% set limits correctly
end
axis([min(xgrid) max(xgrid) 0 1])
set(gca,'xscale','log','XTick',[0.1 1 10],'TickDir','out')

xlabel('speed change (DPS)')
ylabel('Hit Rate')
axis square

subplot(2,2,2)

for i = 1:2
    
   
    scatter(0, All{i}.Output.FA.all.FA,'MarkerEdgeColor',colors(i,:))
    hold on
    plot([0 0],  All{i}.Output.FA.all.FA_confi, 'color',colors(i,:))
    
end
ylim([0 0.2])

set(gca,'XTick', 0,'TickDir','out')
axis square
ylabel('FA')

% subplot(2,2,3)
% colors =[0.8 0.8 0.8; 0.6 0.6 0.6; 0.4 0.4 0.4; 0.2 0.2 0.2;0 0 0]; 
% for i_speed =1:length(speed)
%     
%     h(i_speed) =cdfplot(Ouput.rawRT{i_speed,1}(Ouput.rawRT{i_speed,1}>=0&Ouput.rawRT{i_speed,1}<=700));
%     %h(i_speed) =cdfplot(ou.RTonHit{i_speed});
%     h(i_speed).Color = colors(i_speed,:); 
%     hold on
%     grid off
% end 
% ylim([0 1])
% xlim([0 800])
% vline([200 700],'k:')
% 
% set(gca,'TickDir','out')
% xlabel('Reaction time (ms)')
% axis square

subplot(2,2,3)
for i=1:2
ou = All{i}.Output.target;
errorbar(All{i}.speed,ou.RTonHit_mean(1,:),ou.RTonHit_mean(2,:),'marker','o', 'color',colors(i,:));
hold on
end
xlim([0 40])
ylim([100 700])
set(gca,'TickDir','out')
xlabel('speed change (DPS)')
ylabel('mean RT (ms)')
axis square
print(['W:\All_Staff\home\miaomiao\Analysis\Behavior\select and analysis\speed\base10\747compare'],'-dpdf','-fillpage')
