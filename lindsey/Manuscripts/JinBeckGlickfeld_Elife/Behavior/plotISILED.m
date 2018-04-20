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
    for i_off = 1:length(Output{i}.Infor.Off)
        trialVec = (Output{i}.target.HT_num{i_off,1}+Output{i}.target.Miss_num{i_off,1})';
        trialVpre = (Output{i}.pre.HT_num{i_off,1}+Output{i}.pre.Miss_num{i_off,1})';
        fitS{i,1}{i_off,1} = weibullFitLG(Output{i}.Infor.Orien, Output{i}.target.c_hit{i_off,1}',1, 1, {'nTrials', trialVec}); % use
        fitPre{i,1}{i_off,1} = weibullFitLG(Output{i}.Infor.Orien, Output{i}.pre.c_hit{i_off,1}',1, 1, {'nTrials', trialVpre});
        Thresh(i,i_off) = fitS{i,1}{i_off,1}.thresh;
        Threshpre(i,i_off) = fitPre{i,1}{i_off,1}.thresh;
        
        
        [bootStats{i,1}{i_off,1}]=BootstrapWeibullFit(trialVec, Output{i}.target.c_hit{i_off,1}',NBootstrapReps,Output{i}.Infor.Orien,1, 1);
        [bootStatspre{i,1}{i_off,1}]=BootstrapWeibullFit(trialVpre,Output{i}.pre.c_hit{i_off,1}',NBootstrapReps,Output{i}.Infor.Orien,1, 1);
        Thresh_ci(i,i_off,1:2) =  bootStats{i,1}{i_off,1}.ci95;
        Threshpre_ci(i,i_off,1:2) =  bootStatspre{i,1}{i_off,1}.ci95;
        
    end
    % collapsed all the offs
    Hits = Output{i}.target.all.c_hit;
    trialAll = Output{i}.target.all.HT_num + Output{i}.target.all.Miss_num;
    fitSall{i,1} = weibullFitLG(Output{i}.Infor.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Threshall(i) = fitSall{i,1}.thresh;
    
    [bootStatsall{i,1}]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,Output{i}.Infor.Orien,1, 1);
    Threshall_ci(i,1:2) = bootStatsall{i,1}.ci95;
    
end
%% convert all variables into structure
FitStats = v2struct(fitS,fitSall,fitPre,Thresh,Threshpre,bootStats,bootStatspre,Thresh_ci,Threshpre_ci,Hits,bootStatsall,Threshall_ci);


%% plot LED trials vs non led trials

colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
ledcolor = {[0 0 0] [0.2 0.6 1] [1 0 0.6] [0 0.6 0.2]};
%% first plot the collapsed version of hit (fitting), fa rate,c and d
F=figure;
set(F,'position',[100 50 600 600])
supertitle(['i' num2str(Output{1}.Infor.ID)]);
% supertitle({'  ',['i' num2str(Output{1}.Infor.ID) '-' 'target']});
a=5;

subplot(2,2,1)
for i_led = 1:2
    
    ou = Output{i_led}.target;
    Orien = Output{i_led}.Infor.Orien;
    
    
    maxI = max(Orien);
    minI = min(Orien);
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
    
    %     for i_orien = 1:length(Orien)
    %
    %         line([Orien(i_orien) Orien(i_orien)],ou.all.confi(i_orien,1:2), 'Color',ledcolor{i_led})
    %         hold on
    %     end
    
    if i_led ==1
        text(a+2,1-i_led.*0.1,'Control','Color',ledcolor{i_led})
    else
        text(a+2,1-i_led.*0.1,'LM-0.8mW','Color',ledcolor{i_led})
    end
    
    
    
end
axis([a 100 0 1])

set(gca,'xscale','log','XTick', [10:10:100], 'XTickLabel', {'10';'';'';'';' ';'';' ';' ';'';'100'},'TickDir','out')

xlabel('Orientation change degree')
ylabel('Hit Rate')

subplot(2,2,2)
for i_led = 1:2
    
    scatter(0, Output{i_led}.FA.all.FA,'MarkerEdgeColor',ledcolor{i_led})
    hold on
    plot([0 0], Output{i_led}.FA.all.FA_confi, 'Color',ledcolor{i_led})
    
end
ylim([0 0.2])

set(gca,'XTick', [0], 'XTickLabel', {'Collapsed FA'},'TickDir','out')

ylabel('FA rate')

subplot(2,2,3)
% plot d'
for i_led = 1:2
    Orien = Output{i_led}.Infor.Orien;
    plot(Orien, Output{i_led}.sdt.all.dprime,'Color',ledcolor{i_led})
    hold on
    scatter(Orien,Output{i_led}.sdt.all.dprime,'MarkerEdgeColor',ledcolor{i_led},'MarkerFaceColor',[1 1 1])
    
end
axis([a 100 0 1])

set(gca,'xscale','log','XTick', [10:10:100], 'XTickLabel', {'10';'20';'30';'40';' ';'60';' ';' ';'90';' '},'TickDir','out')

xlabel('Orientation change degree')
ylabel('d prime')
ylim([0 4])

subplot(2,2,4)
% plot criterion
for i_led = 1:2
    Orien = Output{i_led}.Infor.Orien;
    plot(Orien, Output{i_led}.sdt.all.criterion,'Color',ledcolor{i_led})
    hold on
    scatter(Orien,Output{i_led}.sdt.all.criterion,'MarkerEdgeColor',ledcolor{i_led},'MarkerFaceColor',[1 1 1])
    
end
axis([a 100 -0.5 2])

set(gca,'xscale','log','XTick', [10:10:100], 'XTickLabel', {'10';'20';'30';'40';' ';'60';' ';' ';'90';' '},'TickDir','out')

xlabel('Orientation change degree')
ylabel('criterion')


%% plot the threshold difference between conditions and FA rate




F=figure;
set(F,'position',[200 200 600 500])
supertitle(['i' num2str(Output{1}.Infor.ID)]);
h={};
subplot(2,2,1)
for i_led = 1:2
    output = Output{i_led};
    
    Off = output.Infor.Off;
    scatter(Off, FitStats.Thresh(i_led,:),'MarkerEdgeColor',ledcolor{i_led});
    hold on
    
    for i_off=1:length(Off)
        
        line([Off(i_off) Off(i_off)],squeeze(FitStats.Thresh_ci(i_led,i_off,1:2)),'Color',ledcolor{i_led})
        
    end
     if i_led ==1
        text(250,20-i_led.*5,'Control','Color',ledcolor{i_led})
    else
        text(250,20-i_led.*5,'V1-0.3mW','Color',ledcolor{i_led})
    end
end

xlim([200 800])
 ylim([0 50])
xlabel('ISI(ms)')
ylabel('Threshold')
set(gca,'XTick',Off,'TickDir','out')
subplot(2,2,2)

for i_led = 1:2
    output = Output{i_led};
    
    Off = output.Infor.Off;
   scatter(Off, FitStats.Thresh(i_led,:)./FitStats.Thresh(i_led,1),'MarkerEdgeColor',ledcolor{i_led});
    hold on
    
    
    
end

xlim([200 800])
ylim([0 1])
xlabel('ISI(ms)')
ylabel('Threshold-norm.to 250 ms')
set(gca,'XTick',Off,'TickDir','out')



subplot(2,2,3)
for i_led = 1:2
    output = Output{i_led};
    Orien = output.Infor.Orien;
    Off = output.Infor.Off;
    scatter(Off, output.FA.FA,'MarkerEdgeColor',ledcolor{i_led});
    hold on
    
    
    for i_off=1:length(Off)
        
        line([Off(i_off) Off(i_off)],output.FA.FA_confi(i_off,1:2),'Color',ledcolor{i_led})
        
    end
end

xlim([200 800])
ylim([0 0.3])
xlabel('ISI (ms)')
ylabel('FA rate')
set(gca,'XTick',Off,'TickDir','out')


subplot(2,2,4)
for i_led = 1:2
    output = Output{i_led};
    Orien = output.Infor.Orien;
    Off = output.Infor.Off;
    scatter(Off, output.FA.FA./output.FA.FA(1),'MarkerEdgeColor',ledcolor{i_led});
    hold on
    
    
end

xlim([200 800])
ylim([0 5])
xlabel('ISI (ms)')
ylabel('FA rate-norm. to 250ms')
set(gca,'XTick',Off,'TickDir','out')



%% plot the target hit rate

F=figure;
set(F,'position',[50 200 600 500])
supertitle({'  ',['i' num2str(Output{1}.Infor.ID) '-' 'target']});
a=10;
% a=8;


for i_led = 1:2
    subplot(2,2,i_led)
    ou = Output{i_led}.target;
    Orien = Output{i_led}.Infor.Orien;
    Off = Output{i_led}.Infor.Off;
    for i_off =[1 3]%1: length(Off)
        % plot(Orien, S_hit{i_off,1},'Color',colorchoice{i_off})
        h{i_off} = plot(Orien, ou.HT_rate{i_off,1}, 'Color',colorchoice{i_off});
        
        hold on
        
        
    end
    for i_off =[1 3]%1: length(Off)
        %scatter(Orien, S_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
        scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
        hold on
        for i_orien = 1:length(Orien)
            % line([Orien(i_orien) Orien(i_orien)],S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
            line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
            hold on
        end
        
    end
    legend([h{1} h{3}], '250ms','750ms','Location','Southeast')
    legend('boxoff')
    
    axis([a 93 0 1])
    
    set(gca,'xscale','log','XTick', round(Output{1}.Infor.Orien))
    
    xlabel('Orientation change degree')
    ylabel('Hit Rate')
    if i_led==1
        title('Hit Rate-Control')
    else
        title('Hit Rate-AL inhibition')
    end
    
end
subplot(2,2,3)
for i_off = 1
    for i_led = 1:2
        ou = Output{i_led}.target;
        Orien = Output{i_led}.Infor.Orien;
        Off = Output{i_led}.Infor.Off;
        h{i_led} = plot(Orien, ou.HT_rate{i_off,1}, 'Color',ledcolor{i_led});
        hold on
        scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',ledcolor{i_led});
        for i_orien = 1:length(Orien)
            line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',ledcolor{i_led})
            hold on
        end
    end
end
legend([h{1} h{2}], 'Control','AL inhibition','Location','Southeast')
legend('boxoff')

axis([a 93 0 1])
set(gca,'xscale','log','XTick', round(Output{1}.Infor.Orien))
xlabel('Orientation change degree')
ylabel('Hit Rate')
title('Hit Rate-250ms off')

subplot(2,2,4)
for i_off = 3
    for i_led = 1:2
        ou = Output{i_led}.target;
        Orien = Output{i_led}.Infor.Orien;
        Off = Output{i_led}.Infor.Off;
        h{i_led} = plot(Orien, ou.HT_rate{i_off,1}, 'Color',ledcolor{i_led});
        hold on
        scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',ledcolor{i_led});
        for i_orien = 1:length(Orien)
            line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',ledcolor{i_led})
            hold on
        end
    end
end
legend([h{1} h{2}], 'Control','AL inhibition','Location','best')
legend('boxoff')

axis([a 93 0 1])
set(gca,'xscale','log','XTick', round(Output{1}.Infor.Orien))
xlabel('Orientation change degree')
ylabel('Hit Rate')
title('Hit Rate-750ms off')



%% plot FA rate LED vs non LED on FA separated by previous N, FA separated
% by n-1, FA in 100ms window

F=figure;
set(F,'position',[200 200 600 500])
supertitle({['i' num2str(Output{1}.Infor.ID) '-' 'FA rate']});

subplot(2,2,1)
for i_led = 1:2
    output = Output{i_led};
    Orien = output.Infor.Orien;
    Off = output.Infor.Off;
    h{i_led} = scatter(Off, output.FA.FA,'MarkerEdgeColor',ledcolor{i_led});
    hold on
    
    
    for i_off=1:length(Off)
        
        line([Off(i_off) Off(i_off)],output.FA.FA_confi(i_off,1:2),'Color',ledcolor{i_led})
        
    end
end
legend([h{1} h{2}], 'Control','AL-inhibition','Location','Northwest')
legend('boxoff')
xlim([200 800])
ylim([0 0.2])
xlabel('n offs (ms)')
ylabel('FA rate')
set(gca,'XTick',Off)
%  title('FA - All n-1 trials')

subplot(2,2,2)
for i_led = 1:2
    output = Output{i_led};
    Orien = output.Infor.Orien;
    Off = output.Infor.Off;
    h{i_led} = scatter(Off, output.FA.FA./output.FA.FA(1),'MarkerEdgeColor',ledcolor{i_led});
    hold on
    
    
end
legend([h{1} h{2}], 'Control','AL-inhibition','Location','Northwest')
legend('boxoff')
xlim([200 800])
ylim([0 4])
xlabel('n offs (ms)')
ylabel('FA rate-normalized')
set(gca,'XTick',Off)

subplot(2,2,3)

for i_led = 1:2
    output = Output{i_led};
    Orien = output.Infor.Orien;
    Off = output.Infor.Off;
    scatter(Off, output.con.pre2.FA,'MarkerEdgeColor',ledcolor{i_led})
    hold on
    for i_off=1:length(Off)
        
        line([Off(i_off) Off(i_off)],output.con.pre2.FA_confi(i_off,1:2),'Color',ledcolor{i_led})
        
    end
end

ylim([0 0.2])
xlim([200 800])
set(gca,'XTick',Off)
ylabel('FA Rate')
xlabel('n-1 offs (ms)')
%      title('FA - All n trials')

subplot(2,2,4)

for i_led = 1:2
    output = Output{i_led};
    Orien = output.Infor.Orien;
    Off = output.Infor.Off;
    scatter(Off, output.Re.Re,'MarkerEdgeColor',ledcolor{i_led})
    hold on
    for i_off=1:length(Off)
        
        line([Off(i_off) Off(i_off)],output.Re.Re_confi(i_off,1:2),'Color',ledcolor{i_led})
        
    end
end
xlim([200 800])
ylim([0 0.2])
ylabel('Rate within 100ms ')
set(gca,'XTick',Off)
title('release within 100ms stimon')


