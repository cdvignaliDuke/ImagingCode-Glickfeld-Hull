%% try to fit the threshold for different off condiiond
fitS = {};
fitSall = {};
Threshall = [];
fitPre = {};
Thresh = [];
Threshpre = [];
bootStats={};
bootStatsall={};
bootStatspre = {};
NBootstrapReps = 1000;
Thresh_ci = [];
Threshpre_ci = [];
for i = 1:2

    Hits = Output{i}.target.c_hit;
    trialAll = Output{i}.target.HT_num + Output{i}.target.Miss_num;
    fitSall{i,1} = weibullFitLG(Output{i}.Infor.Contrast, Hits',1, 1, {'nTrials', trialAll'});
    Threshall(i) = fitSall{i,1}.thresh;
    
    [bootStatsall{i,1}]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,Output{i}.Infor.Contrast,1, 1);
    Threshall_ci(i,1:2) = bootStatsall{i,1}.ci95;
    
end
%% convert all variables into structure
FitStats = v2struct(fitS,fitSall,fitPre,Thresh,Threshpre,bootStats,bootStatspre,Thresh_ci,Threshpre_ci,Hits,bootStatsall,Threshall_ci);

% plot LED trials vs non led trials

colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
ledcolor = {[0 0 0] [0.2 0.6 1] [1 0 0.6] [0 0.6 0.2]};

%% plot the contrast detection
F=figure;
set(F,'position',[100 50 600 600])
supertitle(['i' num2str(Output{1}.Infor.ID)]);
% supertitle({'  ',['i' num2str(Output{1}.Infor.ID) '-' 'target']});


subplot(2,2,1)
for i_led = 1:2
    
    ou = Output{i_led}.target;
    Contrast = Output{i_led}.Infor.Contrast;
    
    
    maxI = max(Contrast);
    minI = min(Contrast);
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
    
    %     for i_Contrast = 1:length(Contrast)
    %
    %         line([Contrast(i_Contrast) Contrast(i_Contrast)],ou.confi(i_Contrast,1:2), 'Color',ledcolor{i_led})
    %         hold on
    %     end
    
    if i_led ==1
        text(0.01,1-i_led.*0.1,'Control','Color',ledcolor{i_led})
    else
        text(0.01,1-i_led.*0.1,'PM-0.8mW','Color',ledcolor{i_led})
    end
    
    
    
end
axis([0 1 0 1])

set(gca,'xscale','log','XTick',0.1:0.1:1, 'XTickLabel', {'0.1';' ';' ';' ';' 0.5';' ';' ';' ';'';'1 '},'TickDir','out')

xlabel('Contrast')
ylabel('Hit Rate')

subplot(2,2,2)
for i_led = 1:2
    
   % scatter(0, Output{i_led}.FA_percent.FA,'MarkerEdgeColor',ledcolor{i_led})
    scatter(0, Output{i_led}.FA.all.FA,'MarkerEdgeColor',ledcolor{i_led})
    hold on
    plot([0 0], Output{i_led}.FA.all.FA_confi, 'Color',ledcolor{i_led})
    
end
ylim([0 0.3])

set(gca,'XTick', [0], 'XTickLabel', {'Collapsed FA'},'TickDir','out')

ylabel('FA /All trials')

subplot(2,2,3)
% plot the trial trial length dependence stuff
for i_led = 1:2
    
    plot(Output{i_led}.Infor.Contrast, Output{i_led}.sdt.dprime,'Color',ledcolor{i_led})
    hold on
    
    
end
axis([0 1 -0.5 4])

set(gca,'xscale','log','XTick',0.1:0.1:1, 'XTickLabel', {'0.1';' ';' ';' ';' 0.5';' ';' ';' ';'';'1 '},'TickDir','out')

xlabel('Contrast')
ylabel('d prime')


subplot(2,2,4)
% plot criterion
for i_led = 1:2
    
    plot(Output{i_led}.Infor.Contrast, Output{i_led}.sdt.criterion,'Color',ledcolor{i_led})
    hold on
    
    
end
axis([0 1 -0.5 4])

set(gca,'xscale','log','XTick',0.1:0.1:1, 'XTickLabel', {'0.1';' ';' ';' ';' 0.5';' ';' ';' ';'';'1 '},'TickDir','out')

xlabel('Contrast')
ylabel('criterion')

%% plot contrast detection for trial history 

Colors = [0 0 0;0.2 0.8 0;1 0 0.6]; % control, early, miss
ID = Output{1}.Infor.ID;
a=0;
figure
% control trial hit rate dependend on previous outcome
supertitle(['Mouse: i' num2str(ID) '  Led region:PM' ]) %' Fixed-ISI'
subplot(3,2,1)
Contrast = Output{1}.Outcome{1}.target.S_pre.Contrast;
Hit = Output{1}.Outcome{1}.target.S_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{1}.target.S_pre.confi(i_contrast,:),'Color',Colors(1,:))
end 
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(1,:),'MarkerFaceColor',[1 1 1])

Contrast = Output{1}.Outcome{1}.target.E_pre.Contrast;
Hit = Output{1}.Outcome{1}.target.E_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{1}.target.E_pre.confi(i_contrast,:),'Color',Colors(2,:))
end 
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',[1 1 1])

Contrast = Output{1}.Outcome{1}.target.M_pre.Contrast;
Hit = Output{1}.Outcome{1}.target.M_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{1}.target.M_pre.confi(i_contrast,:),'Color',Colors(3,:))
end
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',[1 1 1])
% text(10.5,0.95,'Pre-Success','Color',Colors(1,:))
% text(10.5,0.8,'Pre-Early','Color',Colors(2,:))
% text(10.5,0.65,'Pre-Miss','Color',Colors(3,:))
ylim([0 1])
xlim([a 1])
set(gca,'xscale','log','XTick', 0.10:0.1:1,'XTicklabel',{'0.1','0.2','0.3','','0.5','','','','','1'},'TickDir','out')
xlabel('Contrast')
ylabel('Hit rate')
title('Ctrl/Ctrl')

subplot(3,2,2)
Contrast = Output{1}.Outcome{2}.target.CtrlS_pre.Contrast;
Hit = Output{1}.Outcome{2}.target.CtrlS_pre.c_hit;
hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{2}.target.CtrlS_pre.confi(i_contrast,:),'Color',Colors(1,:))
end 
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(1,:),'MarkerFaceColor',[1 1 1])

Contrast = Output{1}.Outcome{2}.target.CtrlE_pre.Contrast;
Hit = Output{1}.Outcome{2}.target.CtrlE_pre.c_hit;
hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{2}.target.CtrlE_pre.confi(i_contrast,:),'Color',Colors(2,:))
end 
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',[1 1 1])

Contrast = Output{1}.Outcome{2}.target.CtrlM_pre.Contrast;
Hit = Output{1}.Outcome{2}.target.CtrlM_pre.c_hit;
hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{2}.target.CtrlM_pre.confi(i_contrast,:),'Color',Colors(3,:))
end
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',[1 1 1])

% text(10.5,0.95,'Pre-Success','Color',Colors(1,:))
% text(10.5,0.85,'Pre-Early','Color',Colors(2,:))
% text(10.5,0.75,'Pre-Miss','Color',Colors(3,:))
ylim([0 1])
xlim([a 1])
set(gca,'xscale','log','XTick', 0.10:0.1:1,'XTicklabel',{'0.1','0.2','0.3','','0.5','','','','','1'},'TickDir','out')
xlabel('Contrast')
ylabel('Hit rate')
title('Ctrl\Led')

subplot(3,2,3)
Contrast = Output{1}.Outcome{1}.target.LEDS_pre.Contrast;
Hit = Output{1}.Outcome{1}.target.LEDS_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{1}.target.LEDS_pre.confi(i_contrast,:),'Color',Colors(1,:))
end 
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(1,:),'MarkerFaceColor',[1 1 1])

Contrast = Output{1}.Outcome{1}.target.LEDE_pre.Contrast;
Hit = Output{1}.Outcome{1}.target.LEDE_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{1}.target.LEDE_pre.confi(i_contrast,:),'Color',Colors(2,:))
end 
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',[1 1 1])

Contrast = Output{1}.Outcome{1}.target.LEDM_pre.Contrast;
Hit = Output{1}.Outcome{1}.target.LEDM_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{1}.target.LEDM_pre.confi(i_contrast,:),'Color',Colors(3,:))
end
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',[1 1 1])
% text(10.5,0.95,'Pre-Success','Color',Colors(1,:))
% text(10.5,0.85,'Pre-Early','Color',Colors(2,:))
% text(10.5,0.75,'Pre-Miss','Color',Colors(3,:))
ylim([0 1])
xlim([a 1])
set(gca,'xscale','log','XTick', 0.10:0.1:1,'XTicklabel',{'0.1','0.2','0.3','','0.5','','','','','1'},'TickDir','out')
xlabel('Contrast')
ylabel('Hit rate')
title('Led\Ctrl')

subplot(3,2,4)
Contrast = Output{1}.Outcome{2}.target.S_pre.Contrast;
Hit = Output{1}.Outcome{2}.target.S_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{2}.target.S_pre.confi(i_contrast,:),'Color',Colors(1,:))
end 
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(1,:),'MarkerFaceColor',[1 1 1])
Contrast = Output{1}.Outcome{2}.target.E_pre.Contrast;
Hit = Output{1}.Outcome{2}.target.E_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{2}.target.E_pre.confi(i_contrast,:),'Color',Colors(2,:))
end 
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',[1 1 1])
Contrast = Output{1}.Outcome{2}.target.M_pre.Contrast;
Hit = Output{1}.Outcome{2}.target.M_pre.c_hit;

hold on
for i_contrast = 1:length(Contrast)
    plot([Contrast(i_contrast) Contrast(i_contrast)],Output{1}.Outcome{2}.target.M_pre.confi(i_contrast,:),'Color',Colors(3,:))
end
scatter(Contrast,Hit,'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',[1 1 1])
% text(10.5,0.95,'Pre-Success','Color',Colors(1,:))
% text(10.5,0.85,'Pre-Early','Color',Colors(2,:))
% text(10.5,0.75,'Pre-Miss','Color',Colors(3,:))
ylim([0 1])
xlim([a 1])
set(gca,'xscale','log','XTick', 0.10:0.1:1,'XTicklabel',{'0.1','0.2','0.3','','0.5','','','','','1'},'TickDir','out')
xlabel('Contrast')
ylabel('Hit rate')
title('LED/LED')
subplot(3,2,5:6)
hold on
plot([1 1],Output{1}.Outcome{1}.FAB.S_pre.all.FA_confi,'Color',Colors(1,:))
scatter(1,Output{1}.Outcome{1}.FAB.S_pre.all.FA,'MarkerEdgeColor',Colors(1,:),'MarkerFaceColor',[1 1 1])
plot([1 1],Output{1}.Outcome{1}.FAB.E_pre.all.FA_confi,'Color',Colors(2,:))
scatter(1,Output{1}.Outcome{1}.FAB.E_pre.all.FA,'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',[1 1 1])
plot([1 1],Output{1}.Outcome{1}.FAB.M_pre.all.FA_confi,'Color',Colors(3,:))
scatter(1,Output{1}.Outcome{1}.FAB.M_pre.all.FA,'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',[1 1 1])

plot([2 2],Output{1}.Outcome{2}.FAB.CtrlS_pre.all.FA_confi,'Color',Colors(1,:))
scatter(2,Output{1}.Outcome{2}.FAB.CtrlS_pre.all.FA,'MarkerEdgeColor',Colors(1,:),'MarkerFaceColor',[1 1 1])
plot([2 2],Output{1}.Outcome{2}.FAB.CtrlE_pre.all.FA_confi,'Color',Colors(2,:))
scatter(2,Output{1}.Outcome{2}.FAB.CtrlE_pre.all.FA,'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',[1 1 1])
plot([2 2],Output{1}.Outcome{2}.FAB.CtrlM_pre.all.FA_confi,'Color',Colors(3,:))
scatter(2,Output{1}.Outcome{2}.FAB.CtrlM_pre.all.FA,'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',[1 1 1])

plot([3 3],Output{1}.Outcome{1}.FAB.LEDS_pre.all.FA_confi,'Color',Colors(1,:))
scatter(3,Output{1}.Outcome{1}.FAB.LEDS_pre.all.FA,'MarkerEdgeColor',Colors(1,:),'MarkerFaceColor',[1 1 1])
plot([3 3],Output{1}.Outcome{1}.FAB.LEDE_pre.all.FA_confi,'Color',Colors(2,:))
scatter(3,Output{1}.Outcome{1}.FAB.LEDE_pre.all.FA,'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',[1 1 1])
plot([3 3],Output{1}.Outcome{1}.FAB.LEDM_pre.all.FA_confi,'Color',Colors(3,:))
scatter(3,Output{1}.Outcome{1}.FAB.LEDM_pre.all.FA,'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',[1 1 1])


plot([4 4],Output{1}.Outcome{2}.FAB.S_pre.all.FA_confi,'Color',Colors(1,:))
scatter(4,Output{1}.Outcome{2}.FAB.S_pre.all.FA,'MarkerEdgeColor',Colors(1,:),'MarkerFaceColor',[1 1 1])
plot([4 4],Output{1}.Outcome{2}.FAB.E_pre.all.FA_confi,'Color',Colors(2,:))
scatter(4,Output{1}.Outcome{2}.FAB.E_pre.all.FA,'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',[1 1 1])
plot([4 4],Output{1}.Outcome{2}.FAB.M_pre.all.FA_confi,'Color',Colors(3,:))
scatter(4,Output{1}.Outcome{2}.FAB.M_pre.all.FA,'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',[1 1 1])
text(4.5,0.16,'Pre-Success','Color',Colors(1,:))
text(4.5,0.12,'Pre-Early','Color',Colors(2,:))
text(4.5,0.08,'Pre-Miss','Color',Colors(3,:))

ylim([0 0.2])
ylabel('FA rate')
xlim([0.5 5.5])
set(gca,'XTick', 1:1:4,'XTicklabel',{'Ctrl\Ctrl','Ctrl\LED','LED\Ctrl','LED\LED'},'TickDir','out','YTick',0:0.1:0.5)