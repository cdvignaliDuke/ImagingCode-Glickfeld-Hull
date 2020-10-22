%% try to fit the threshold for different off condiiond

fitSall = {};
Threshall = [];
bootStats={};
bootStatsall={};
bootStatspre = {};
NBootstrapReps = 1000;
Threshall_ci = [];

for i = 1:2

    Hits = Output{i}.target.c_hit;
    trialAll = Output{i}.target.HT_num + Output{i}.target.Miss_num;
    fitSall{i,1} = weibullFitLG(Output{i}.Infor.speed, Hits',1, 1, {'nTrials', trialAll'});
    Threshall(i) = fitSall{i,1}.thresh;
    
    [bootStatsall{i,1}]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,Output{i}.Infor.speed,1, 1);
    Threshall_ci(i,1:2) = bootStatsall{i,1}.ci95;
    
end
speed = Output{1}.Infor.speed;
RT_bins = 200:50:701;
Color = flipud(copper(5));
RT_pdf = [];
beta = [];
R = [];
MSE = [];
for i_led = 1:2
for i_speed = 1:length(speed)
     RT_pdf(i_led,i_speed,:) = (hist(Output{i_led}.target.RTonHit{i_speed},RT_bins)./length(Output{i_led}.target.RTonHit{i_speed}))./50; % transfer to density
end 
%[beta(i_led,:),R(i_led,:),MSE(i_led,:)]=eLATER( squeeze(RT_pdf(i_led,:,:)),Color);
end 

%% convert all variables into structure
FitStats = v2struct(fitSall,bootStats,bootStatsall,Threshall_ci);

% plot LED trials vs non led trials

colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
ledcolor = {[0 0 0] [0.2 0.6 1] [1 0 0.6] [0 0.6 0.2]};

%% plot the speed detection
F=figure;

supertitle(['i' num2str(Output{1}.Infor.ID)]);
% supertitle({'  ',['i' num2str(Output{1}.Infor.ID) '-' 'target']});


subplot(2,2,1)
for i_led = 1:2
    
    ou = Output{i_led}.target;
    speed = Output{i_led}.Infor.speed;
    
    
    maxI = max(speed);
    minI = min(speed);
    xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
    h=line(xgrid, fitSall{i_led,1}.modelFun(fitSall{i_led,1}.coefEsts, xgrid), 'Color',ledcolor{i_led});
    hold on;
    plot(fitSall{i_led,1}.intensityX,fitSall{i_led,1}.fractCorrY, 'o','Color',ledcolor{i_led});
    %thresh = coefEsts(1)*[1 1];
    plot(fitSall{i_led,1}.thresh*[1 1], [0 fitSall{i_led,1}.threshY], '--','Color',ledcolor{i_led});
    plot(bootStatsall{i_led,1}.ci95, fitSall{i_led,1}.threshY*[1 1], 'Color',ledcolor{i_led});
    
    % set limits correctly
    xLim = [min(xgrid) max(xgrid)].* [0.75 1.25];
    xLim = 10.^ceil(log10(xLim) - [1 0]);
 
    if i_led ==1
        text(0.2,1-i_led.*0.2,'Control','Color',ledcolor{i_led})
    else
        text(0.2,1-i_led.*0.2,'LM-0.5mW','Color',ledcolor{i_led})
    end
    
    
    
end

axis([min(xgrid) max(xgrid) 0 1])
set(gca,'xscale','log','XTick',[0.1 1 10],'TickDir','out')

xlabel('Speed Increment (dps)')
ylabel('Hit Rate')
axis square

subplot(2,2,2)
for i_led = 1:2
    
   % scatter(0, Output{i_led}.FA_percent.FA,'MarkerEdgeColor',ledcolor{i_led})
    scatter(0, Output{i_led}.FA.all.FA,'MarkerEdgeColor',ledcolor{i_led})
    hold on
    plot([0 0], Output{i_led}.FA.all.FA_confi, 'Color',ledcolor{i_led})
    
end
ylim([0 0.4])

set(gca,'XTick', 0,'TickDir','out')
axis square
ylabel('FA')

% plot the reaction time pdf? 
% subplot(2,2,3)
% temp = flipud(gray(6));
% colors(1,:,:) =temp(2:end,:); 
% temp =brewermap(6,'Blues'); 
% colors(2,:,:) = temp(2:end,:);
% 
% 
% for i_led = 1:2
%     for i_speed = 1:length(speed)
%         h = cdfplot(Output{i_led}.target.RTonHit{i_speed});
%         set(h, 'Color',squeeze(colors(i_led,i_speed,:)));
%         hold on
%     end 
% end 
% xlim([200 700])
% ylim([0 1])
% xlabel('Reaction Time (ms)')
% ylabel('Cumulative Probability')
% box off
% grid off
% set(gca,'TickDir','out')
% axis square


subplot(2,2,4)
for i_led = 1:2
    errorbar(Output{i_led}.Infor.speed, Output{i_led}.target.RTonHit_mean(1,:),Output{i_led}.target.RTonHit_mean(2,:),'Color',ledcolor{i_led},'Marker','o')
    hold on
end 
box off
xlim([-2 34])
ylim([300 600])
set(gca,'TickDir','out','XTick',0:10:30)
xlabel('Speed Increment (sps)')
ylabel('mean RT (ms)')
axis square
%% plot the fitted reaction time
figure
subplot(2,2,1)
temp = flipud(gray(6));
colors(1,:,:) =temp(2:end,:); 
temp =brewermap(6,'Blues'); 
colors(2,:,:) = temp(2:end,:);
t1 = 1:10:800;
for i_led = 1:2
   
for i = 1:length(speed)
    scatter(RT_bins,squeeze(RT_pdf(i_led,i,:)),'o','MarkerFaceColor',squeeze(colors(i_led,i,:)),'MarkerEdgeColor','none')
    hold on
line(t1, ((t1 + (beta(i_led,i).*((beta(i_led,11)./beta(i_led,5+i)).^2)))./(((t1.^2)+ ((beta(i_led,11)./beta(i_led,5+i)).^2)).^1.5)).*(1./(sqrt(2.*pi).*beta(i_led,5+i))).*exp(-(1./(2.*(beta(i_led,5+i).^2))).*((1./t1 - beta(i_led,i)).^2).*((t1.^2)./(t1.^2 +(beta(i_led,11)./beta(i_led,5+i)).^2))),'color',squeeze(colors(i_led,i,:)),'linewidth',1)
end

end 
set(gca,'TickDir','out')
xlim([0 800])
ylim([0 0.009])
xlabel('Reaction Time (ms)')
ylabel('Probability')
axis square

subplot(2,2,2)
% plot drift rate mean/sd
for i_led = 1:2
gscatter(speed,squeeze(beta(i_led,1:5))./squeeze(beta(i_led,6:10)),[1 2 3 4 5],squeeze(colors(i_led,:,:)))
hold on
box off;
legend off;
end 
xlim([0 40])
% ylim([0 20])
set(gca,'TickDir','out')
ylabel('Drift rate: mean/sd')
axis square

subplot(2,2,3)
for i_led = 1:2
scatter(1,1./beta(i_led,end),'MarkerFaceColor',squeeze(colors(i_led,5,:)),'MarkerEdgeColor','none')
hold on
end 
% ylim([0 15])
ylabel('Evidence: mean/sd')
set(gca,'TickDir','out')
axis square
%% two way anova of reaction time difference

Ctrl = RT(idx,:,1);
LED = RT(idx,:,2);
[p,tbl,stats] = anova2([Ctrl;LED],3)
c= multcompare(stats,'Estimate','row');