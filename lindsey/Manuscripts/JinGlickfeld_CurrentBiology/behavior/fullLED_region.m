% summary of full led inhibition
homeDir = 'Z:\All_Staff\home\miaomiao\Analysis\Behavior\select and analysis';
xlsfile = fullfile(homeDir, 'LED-full.xlsx');
[Gdata, Gtext, Graw] = xlsread(xlsfile,7);
IDs =Gdata(:,1);
Dir = Gtext(2:end,5);
region = Gtext(2:end,2);
power = Gdata(:,3);
background = Gtext(2:end,6);
virus = Gtext(2:end,7);
Thresh = NaN(length(IDs),2); %  first collumn is the control, second is led
FA =  NaN(length(IDs),2);
Fidget = NaN(length(IDs),2);
FA_trial = NaN(length(IDs),2);
dprime = NaN(length(IDs),12,2);
criterion = NaN(length(IDs),12,2);
Oriens = NaN(length(IDs),12,2);
RT = NaN(length(IDs),12,2); % store the reaction time
Pre_thresh = NaN(length(IDs),4,3);% mice * ctrl/ctrl combinations * pre-success, pre-early, pre-miss
Pre_FA = NaN(length(IDs),4,3);
LEDpre_thresh = NaN(length(IDs),4); % all combined: LED/Ctrl; LED/LED; Ctrl/Ctrl; Ctrl/LED
LEDpre_FA = NaN(length(IDs),4); % previouse led trial, followed by control/or followed by led
LEDpre_hit = NaN(length(IDs),4);% for interpolated 22.5 deg
LEDpre_d =  NaN(length(IDs),4);
LEDpre_c = NaN(length(IDs),4);

Pre_fit = {};
TrialID = {};
FA_bintrial = NaN(length(IDs),4,2); % only has four bins
Thresh_bintrial = NaN(length(IDs),4,2);
bootStats = {};
Thresh_pass = NaN(length(IDs),1); % 95 confidence interval lower bound for now 
Thresh_cpass = NaN(length(IDs),1); % combine the success and early to yes catogory
FA_pass = NaN(length(IDs),1);
FA_cpass = NaN(length(IDs),1); % combine the success and early to yes catogory
Hit_22 = NaN(length(IDs),2);
d_22 = NaN(length(IDs),2);
c_22 = NaN(length(IDs),2);

for i_exp = 1:length(IDs)
    load(Dir{i_exp})
    
    % get the threshold for previouse LED, current trial is control (LED/Ctrl)
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.LED_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.LED_pre.all.HT_num + Output{1}.Outcome{1}.target.LED_pre.all.Miss_num;
    LEDpre_fit{i_exp,1} = weibullFitLG(Output{1}.Outcome{1}.target.LED_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,1) = LEDpre_fit{i_exp,1}.thresh;
    LEDpre_FA(i_exp,1) = Output{1}.Outcome{1}.FA.LED_pre.all.FA;
    LEDpre_hit(i_exp,1) = LEDpre_fit{i_exp,1}.modelFun(LEDpre_fit{i_exp,1}.coefEsts, 22.5);
    [LEDpre_d(i_exp,1), LEDpre_c(i_exp,1)] = dprime_simple(LEDpre_hit(i_exp,1),LEDpre_FA(i_exp,1));
    
    % get the threshold for previouse LED, current trial is LED(LED/LED)
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.Same_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.Same_pre.all.HT_num + Output{1}.Outcome{2}.target.Same_pre.all.Miss_num;
    LEDpre_fit{i_exp,2} = weibullFitLG(Output{1}.Outcome{2}.target.Same_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,2) = LEDpre_fit{i_exp,2}.thresh;
    LEDpre_FA(i_exp,2) = Output{1}.Outcome{2}.FA.Same_pre.all.FA;
    LEDpre_hit(i_exp,2) = LEDpre_fit{i_exp,2}.modelFun(LEDpre_fit{i_exp,2}.coefEsts, 22.5); 
    [LEDpre_d(i_exp,2), LEDpre_c(i_exp,2)] = dprime_simple(LEDpre_hit(i_exp,2),LEDpre_FA(i_exp,2));
    
    % get the threshold for previouse control, current trial is control (Ctrl/Ctrl)
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.Same_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.Same_pre.all.HT_num + Output{1}.Outcome{1}.target.Same_pre.all.Miss_num;
    LEDpre_fit{i_exp,3} = weibullFitLG(Output{1}.Outcome{1}.target.Same_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,3) = LEDpre_fit{i_exp,3}.thresh;
    LEDpre_FA(i_exp,3) = Output{1}.Outcome{1}.FA.Same_pre.all.FA;
    LEDpre_hit(i_exp,3) = LEDpre_fit{i_exp,3}.modelFun(LEDpre_fit{i_exp,3}.coefEsts, 22.5);
    [LEDpre_d(i_exp,3), LEDpre_c(i_exp,3)] = dprime_simple(LEDpre_hit(i_exp,3),LEDpre_FA(i_exp,3));
    
    % get the threshold for previouse Control, current trial is LED(Ctrl/LED)
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.LED_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.LED_pre.all.HT_num + Output{1}.Outcome{2}.target.LED_pre.all.Miss_num;
    LEDpre_fit{i_exp,4} = weibullFitLG(Output{1}.Outcome{2}.target.LED_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,4) = LEDpre_fit{i_exp,4}.thresh;
    LEDpre_FA(i_exp,4) = Output{1}.Outcome{2}.FA.LED_pre.all.FA;
    LEDpre_hit(i_exp,4) = LEDpre_fit{i_exp,4}.modelFun(LEDpre_fit{i_exp,4}.coefEsts, 22.5);
    [LEDpre_d(i_exp,4), LEDpre_c(i_exp,4)] = dprime_simple(LEDpre_hit(i_exp,4),LEDpre_FA(i_exp,4));
    
    
    
    %get the threshold for the ctrl/ctrl
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.S_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.S_pre.all.HT_num + Output{1}.Outcome{1}.target.S_pre.all.Miss_num;
    Pre_fit{i_exp,1,1} = weibullFitLG(Output{1}.Outcome{1}.target.S_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,1,1) = Pre_fit{i_exp,1,1}.thresh;
    bootStats{i_exp,1,1}=BootstrapWeibullFit(trialAll',Hits',1000,Output{1}.Outcome{1}.target.S_pre.Orien,1, 1);
    
    
    
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.E_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.E_pre.all.HT_num + Output{1}.Outcome{1}.target.E_pre.all.Miss_num;
    Pre_fit{i_exp,1,2} = weibullFitLG(Output{1}.Outcome{1}.target.E_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,1,2) = Pre_fit{i_exp,1,2}.thresh;
    
    % combine success and early as yes response
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.E_pre.all.HT_num + Output{1}.Outcome{1}.target.E_pre.all.Miss_num + ...
        Output{1}.Outcome{1}.target.S_pre.all.HT_num + Output{1}.Outcome{1}.target.S_pre.all.Miss_num;
    
    Hits = [];
    Hits = (Output{1}.Outcome{1}.target.E_pre.all.HT_num + Output{1}.Outcome{1}.target.S_pre.all.HT_num)./trialAll;
    
    bootStats{i_exp,1,2}=BootstrapWeibullFit(trialAll',Hits',1000,Output{1}.Outcome{1}.target.S_pre.Orien,1, 1);
    
    
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.M_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.M_pre.all.HT_num + Output{1}.Outcome{1}.target.M_pre.all.Miss_num;
    Pre_fit{i_exp,1,3} = weibullFitLG(Output{1}.Outcome{1}.target.M_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,1,3) = Pre_fit{i_exp,1,3}.thresh;
    bootStats{i_exp,1,3}=BootstrapWeibullFit(trialAll',Hits',1000,Output{1}.Outcome{1}.target.M_pre.Orien,1, 1);
    
    
    
    Thresh_pass(i_exp,1) =  bootStats{i_exp,1,3}.ci95(1)> bootStats{i_exp,1,1}.ci95(2);
    Thresh_cpass(i_exp,1) =  bootStats{i_exp,1,3}.ci95(1)> bootStats{i_exp,1,2}.ci95(2);
    
    Pre_FA(i_exp,1,1) = Output{1}.Outcome{1}.FA.S_pre.all.FA;
    Pre_FA(i_exp,1,2) = Output{1}.Outcome{1}.FA.E_pre.all.FA;
    Pre_FA(i_exp,1,3) = Output{1}.Outcome{1}.FA.M_pre.all.FA;
    
    FA_pass(i_exp,1) =  Output{1}.Outcome{1}.FA.M_pre.all.FA_confi(2)<Output{1}.Outcome{1}.FA.S_pre.all.FA_confi(1);
    
    % combine success and early as yes response
    FA_temp = [];
    CR_temp = [];
    confi_temp = [];
    FA_temp = sum(Output{1}.Outcome{1}.FA.E_pre.FA_Num) +  sum(Output{1}.Outcome{1}.FA.S_pre.FA_Num);
    CR_temp =  sum(Output{1}.Outcome{1}.FA.E_pre.CR_Num) +  sum(Output{1}.Outcome{1}.FA.S_pre.CR_Num);
    [~,confi_temp(1:2)] = binofit(FA_temp,FA_temp+CR_temp);
    
    
    FA_cpass(i_exp,1) =  Output{1}.Outcome{1}.FA.M_pre.all.FA_confi(2)<confi_temp(1);
     
     
    TrialID {i_exp,1,1} = Output{1}.Outcome{1}.S_pre_trialID;
    TrialID {i_exp,1,2} = Output{1}.Outcome{1}.E_pre_trialID;
    TrialID {i_exp,1,3} = Output{1}.Outcome{1}.M_pre_trialID;
    
    
    
    % get the threshold for the led/led 
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.S_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.S_pre.all.HT_num + Output{1}.Outcome{2}.target.S_pre.all.Miss_num;
    Pre_fit{i_exp,4,1} = weibullFitLG(Output{1}.Outcome{2}.target.S_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,4,1) = Pre_fit{i_exp,4,1}.thresh;
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.E_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.E_pre.all.HT_num + Output{1}.Outcome{2}.target.E_pre.all.Miss_num;
    Pre_fit{i_exp,4,2} = weibullFitLG(Output{1}.Outcome{2}.target.E_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,4,2) = Pre_fit{i_exp,4,2}.thresh;
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.M_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.M_pre.all.HT_num + Output{1}.Outcome{2}.target.M_pre.all.Miss_num;
    Pre_fit{i_exp,4,3} = weibullFitLG(Output{1}.Outcome{2}.target.M_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,4,3) = Pre_fit{i_exp,4,3}.thresh;
    
    Pre_FA(i_exp,4,1) = Output{1}.Outcome{2}.FA.S_pre.all.FA;
    Pre_FA(i_exp,4,2) = Output{1}.Outcome{2}.FA.E_pre.all.FA;
    Pre_FA(i_exp,4,3) = Output{1}.Outcome{2}.FA.M_pre.all.FA;
    TrialID {i_exp,4,1} = Output{1}.Outcome{2}.S_pre_trialID;
    TrialID {i_exp,4,2} = Output{1}.Outcome{2}.E_pre_trialID;
    TrialID {i_exp,4,3} = Output{1}.Outcome{2}.M_pre_trialID;
    
    % get the threshold for the led/ctrl 
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.LEDS_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.LEDS_pre.all.HT_num + Output{1}.Outcome{1}.target.LEDS_pre.all.Miss_num;
    Pre_fit{i_exp,3,1} = weibullFitLG(Output{1}.Outcome{1}.target.LEDS_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,3,1) = Pre_fit{i_exp,3,1}.thresh;
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.LEDE_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.LEDE_pre.all.HT_num + Output{1}.Outcome{1}.target.LEDE_pre.all.Miss_num;
    Pre_fit{i_exp,3,2} = weibullFitLG(Output{1}.Outcome{1}.target.LEDE_pre.Orien(~isnan(Hits)), Hits(~isnan(Hits))',1, 1, {'nTrials', trialAll(~isnan(Hits))'});
    Pre_thresh(i_exp,3,2) = Pre_fit{i_exp,3,2}.thresh;
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.LEDM_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.LEDM_pre.all.HT_num + Output{1}.Outcome{1}.target.LEDM_pre.all.Miss_num;
    Pre_fit{i_exp,3,3} = weibullFitLG(Output{1}.Outcome{1}.target.LEDM_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,3,3) = Pre_fit{i_exp,3,3}.thresh;
    
    Pre_FA(i_exp,3,1) = Output{1}.Outcome{1}.FA.LEDS_pre.all.FA;
    Pre_FA(i_exp,3,2) = Output{1}.Outcome{1}.FA.LEDE_pre.all.FA;
    Pre_FA(i_exp,3,3) = Output{1}.Outcome{1}.FA.LEDM_pre.all.FA;
    TrialID {i_exp,3,1} = Output{1}.Outcome{1}.LEDS_pre_trialID;
    TrialID {i_exp,3,2} = Output{1}.Outcome{1}.LEDE_pre_trialID;
    TrialID {i_exp,3,3} = Output{1}.Outcome{1}.LEDM_pre_trialID;
    
    % get the threshold for the ctrl/led
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.CtrlS_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.CtrlS_pre.all.HT_num + Output{1}.Outcome{2}.target.CtrlS_pre.all.Miss_num;
    Pre_fit{i_exp,2,1} = weibullFitLG(Output{1}.Outcome{2}.target.CtrlS_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,2,1) = Pre_fit{i_exp,2,1}.thresh;
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.CtrlE_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.CtrlE_pre.all.HT_num + Output{1}.Outcome{2}.target.CtrlE_pre.all.Miss_num;
    Pre_fit{i_exp,2,2} = weibullFitLG(Output{1}.Outcome{2}.target.CtrlE_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    Pre_thresh(i_exp,2,2) = Pre_fit{i_exp,2,2}.thresh;
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.CtrlM_pre.all.c_hit;
    
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.CtrlM_pre.all.HT_num + Output{1}.Outcome{2}.target.CtrlM_pre.all.Miss_num;
    Pre_fit{i_exp,2,3} = weibullFitLG(Output{1}.Outcome{2}.target.CtrlM_pre.Orien(~isnan(Hits)), Hits(~isnan(Hits))',1, 1, {'nTrials', trialAll(~isnan(Hits))'});
    Pre_thresh(i_exp,2,3) = Pre_fit{i_exp,2,3}.thresh;
    
    Pre_FA(i_exp,2,1) = Output{1}.Outcome{2}.FA.CtrlS_pre.all.FA;
    Pre_FA(i_exp,2,2) = Output{1}.Outcome{2}.FA.CtrlE_pre.all.FA;
    Pre_FA(i_exp,2,3) = Output{1}.Outcome{2}.FA.CtrlM_pre.all.FA;
    TrialID {i_exp,2,1} = Output{1}.Outcome{2}.CtrlS_pre_trialID;
    TrialID {i_exp,2,2} = Output{1}.Outcome{2}.CtrlE_pre_trialID;
    TrialID {i_exp,2,3} = Output{1}.Outcome{2}.CtrlM_pre_trialID;
    
    
    for i = 1:2
        Hits = [];
        trialAll = [];
        Hits = Output{i}.target.all.c_hit;
        trialAll = Output{i}.target.all.HT_num + Output{i}.target.all.Miss_num;
        fitSall{i_exp,i} = weibullFitLG(Output{i}.Infor.Orien, Hits',1, 1, {'nTrials', trialAll'});
        Thresh(i_exp,i) = fitSall{i_exp,i}.thresh;
        FA(i_exp,i) = Output{i}.FA.all.FA;
        Fidget(i_exp,i) = Output{i}.FA_percent.fidget;
        FA_trial(i_exp,i) = Output{i}.FA_percent.FA;
        dprime(i_exp,1:length(Output{i}.Infor.Orien),i) = Output{i}.sdt.all.dprime;
        criterion(i_exp,1:length(Output{i}.Infor.Orien),i) = Output{i}.sdt.all.criterion;
        Oriens(i_exp,1:length(Output{i}.Infor.Orien),i) = Output{i}.Infor.Orien;
        RT(i_exp,1:length(Output{i}.Infor.Orien),i) = Output{i}.target.all.RTonHit_mean(1,:);
        FA_bintrial(i_exp,1:4,i)  = Output{i}.bintrial.FA;
        Thresh_bintrial(i_exp,1:4,i)  =  Output{i}.bintrial.thresh;
        % interpolate the 22.5 deg for d' and c
        Hit_22(i_exp,i) = fitSall{i_exp,i}.modelFun(fitSall{i_exp,i}.coefEsts, 22.5); 
        [d_22(i_exp,i), c_22(i_exp,i)] = dprime_simple(Hit_22(i_exp,i),FA(i_exp,i));
%         d_22(i_exp,i) = mean(Output{i}.sdt.all.dprime(Output{i}.Infor.Orien>20 & Output{i}.Infor.Orien<30));
%         c_22(i_exp,i) = mean(Output{i}.sdt.all.criterion(Output{i}.Infor.Orien>20 & Output{i}.Infor.Orien<30));
    end
end
a=squeeze(Oriens(:,:,1));
cbin = Output{1}.bintrial.cbin.Hit; 
%% plot the histogram trial ID difference between pre-success, pre-early, premiss
Colors = [0 0 0;0.2 0.8 0;1 0 0.6]; 
j = 1;
figure;
for i = 1:length(IDs)
    subplot(6,3,i)
    
    histogram(TrialID{i,j,2},14,'FaceColor',Colors(2,:),'EdgeColor',[1 1 1],'FaceAlpha',0.5)
    hold on
    histogram(TrialID{i,j,1},14,'FaceColor',Colors(1,:),'EdgeColor',[1 1 1],'FaceAlpha',0.5)
    histogram(TrialID{i,j,3},14,'FaceColor',Colors(3,:),'EdgeColor',[1 1 1],'FaceAlpha',0.5)
    if i==1
    legend('Pre-early','Pre-Success','Pre-Miss','Location','northeast')
    ylabel('Trial numbers')
    xlabel('X th trial in a session')
    end 
    xlim([0 700])
%     ylim([0 200])
    
    title(['i' num2str(IDs(i)) '  LED area: ' region{i}])
end 

figure
for i = 1:length(IDs)
    subplot(6,3,i)
    
    h=cdfplot(TrialID{i,j,1});
    h.Color = Colors(1,:);
    hold on
    h=cdfplot(TrialID{i,j,2});
    h.Color = Colors(2,:);
    h=cdfplot(TrialID{i,j,3});
    h.Color = Colors(3,:);
   
    if i==1
    legend('Pre-early','Pre-Success','Pre-Miss','Location','northeast')
    
    xlabel('X th trial in a session')
    end 
    xlim([0 700])
    ylim([0 1])
    
    title(['i' num2str(IDs(i)) '  LED area: ' region{i}])
end 

%% summary global average of the normalized value to the previous success condition
Colors = [0 0 0;0.2 0.8 0;1 0 0.6]; % control, early, miss
idx = logical(ones(length(IDs),1));
figure
norm = squeeze(Pre_thresh(idx,1,:))./repmat(squeeze(Pre_thresh(idx,1,1)),1,3);
mean_trl = mean(norm,1);
sem_trl = std(norm,[],1)./sqrt(sum(idx)); 

for i = 1:3
  errorbar(1,mean_trl(i),sem_trl(i),'Color',Colors(i,:))
hold on
scatter([1],mean_trl(i),'MarkerEdgeColor',Colors(i,:),'MarkerFaceColor',[1 1 1],'SizeData',80)
end 


norm = squeeze(Pre_FA(idx,1,:))./repmat(squeeze(Pre_FA(idx,1,1)),1,3);

mean_trl = mean(norm,1);
sem_trl = std(norm,[],1)./sqrt(sum(idx)); 

for i = 1:3
  errorbar(2,mean_trl(i),sem_trl(i),'Color',Colors(i,:))
hold on
scatter([2],mean_trl(i),'MarkerEdgeColor',Colors(i,:),'MarkerFaceColor',[1 1 1],'SizeData',80)
end 


ylim([0 1.6])
xlim([0 3])

% test in one way anova
% [p,t,stats] = anova1(norm);
% multcompare(stats)
% text(1,5,['p = ' num2str(p)])

set(gca,'XTick', 1:1:2,'XTicklabel',{'Norm.Thresh','Norm.FA'},'TickDir','out')
ylabel('Thresh/FA norm. to pre-Success')


%% summary those that pass FA test, check how FA rate change in different LED conditions, seperate by areas
all_region = {'V1','LM','AL','PM'};
figure
i=3;
idx = strcmp(region,all_region(i));%&logical(FA_pass);
Colors = [0 0 0;0.2 0.8 0;1 0 0.6]; % control, early, miss

norm = [];
mean_trl = [];
sem_trl = [];
norm = Pre_FA(idx,:,:)./repmat(squeeze(Pre_FA(idx,1,1)),1,4,3);
mean_trl = squeeze(mean(norm,1));
sem_trl = squeeze(std(norm,[],1))./sqrt(sum(idx)); 

for i_con = 1:4

for j = 1:3
  errorbar(i_con,mean_trl(i_con,j),sem_trl(i_con,j),'Color',Colors(j,:))
hold on
scatter(i_con,mean_trl(i_con,j),'MarkerEdgeColor',Colors(j,:),'MarkerFaceColor',[1 1 1],'SizeData',80)
end 
end 
ylim([0 1.4])
xlim([0.5 4.5])
set(gca,'XTick',1:1:4,'XTicklabel',{'Ctrl/Ctrl','Ctrl/Led','Led/Ctrl','Led/led'},'TickDir','out')
ylabel('FA norm. to pre-success ctrl/ctrl')
title(['area:' all_region{i} '  n = ' num2str(sum(idx)) 'mice'])
%% summary the trial history dependence of FA rate and threshold in ctrl/ctrl
all_region = {'V1','LM','AL','PM'};
i=3;
%idx =strcmp(region,all_region(i));
 idx = logical(ones(length(IDs),1));
 %idx = logical(Thresh_cpass);
Colors = [0 0 0;0.2 0.8 0;1 0 0.6]; 
Color_led = [0 0 0;0.2 0.8 1 ; 0.5 0.5 0.5;0 0.4 0.8];
figure
subplot(2,2,1)
for j=1:4
errorbar([1 2 3],mean(squeeze(Pre_thresh(idx,j,:)),1),std(squeeze(Pre_thresh(idx,j,:)),[],1)./sqrt(length(IDs(idx))),'Color',Color_led(j,:))
hold on
scatter([1 2 3],mean(squeeze(Pre_thresh(idx,j,:)),1),'MarkerEdgeColor',Color_led(j,:),'MarkerFaceColor',[1 1 1])
set(gca,'XTick', 1:1:3,'XTicklabel',{'Success','Early','Miss'},'TickDir','out')

end 
ylim([0 40])
text(1,8,'Ctrl/Ctrl','color',Color_led(1,:))
text(1,5,'Ctrl/Led','color',Color_led(2,:))
text(2,8,'Led/Ctrl','color',Color_led(3,:))
text(2,5,'Led/led','color',Color_led(4,:))
ylabel('Orien thresh(deg)')
xlabel('Previous trial outcome')
% test in one way anova
%[p,t,stats] = anova1(squeeze(Pre_thresh(idx,1,:)));
% text(1,5,['p = ' num2str(p)])

subplot(2,2,2)
for j=1:4
norm = squeeze(Pre_thresh(idx,j,:))./repmat(squeeze(Pre_thresh(idx,1,1)),1,3);
errorbar([1 2 3],mean(norm,1),std(norm,[],1)./sqrt(length(IDs(idx))),'Color',Color_led(j,:))
hold on
scatter([1 2 3],mean(norm,1),'MarkerEdgeColor',Color_led(j,:),'MarkerFaceColor',[1 1 1])
set(gca,'XTick', 1:1:3,'XTicklabel',{'Success','Early','Miss'},'TickDir','out')

end 
ylim([0 2])
ylabel('Norm. orien thresh(deg)')
xlabel('Previous trial outcome')
% [p1,t,stats] = anova1(norm);
% multcompare(stats)
% text(1,0.25,['p = ' num2str(p1)])

subplot(2,2,3)
for j =1:4
errorbar([1 2 3],mean(squeeze(Pre_FA(idx,j,:)),1),std(squeeze(Pre_FA(idx,j,:)),[],1)./sqrt(length(IDs(idx))),'Color',Color_led(j,:))
hold on
scatter([1 2 3],mean(squeeze(Pre_FA(idx,j,:)),1),'MarkerEdgeColor',Color_led(j,:),'MarkerFaceColor',[1 1 1])
set(gca,'XTick', 1:1:3,'XTicklabel',{'Success','Early','Miss'},'TickDir','out')
end 
ylim([0 0.15])
ylabel('FA rate')
xlabel('Previous trial outcome')
% [p2,t,stats] = anova1(Pre_FA);
% multcompare(stats)
% text(1,0.025,['p = ' num2str(p2)])

subplot(2,2,4)
for j=1:4
norm =squeeze(Pre_FA(idx,j,:))./repmat(squeeze(Pre_FA(idx,1,1)),1,3);
errorbar([1 2 3],mean(norm,1),std(norm,[],1)./sqrt(length(IDs(idx))),'Color',Color_led(j,:))
hold on
scatter([1 2 3],mean(norm,1),'MarkerEdgeColor',Color_led(j,:),'MarkerFaceColor',[1 1 1])
set(gca,'XTick', 1:1:3,'XTicklabel',{'Success','Early','Miss'},'TickDir','out')
end 
 ylim([0 1.2])
ylabel('Norm. FA rate')
xlabel('Previous trial outcome')
% [p3,t,stats] = anova1(norm);
% multcompare(stats)
% text(1,0.2,['p = ' num2str(p3)])
%supertitle(['All LED area: ', all_region{i},'   n = ', num2str(sum(idx)), 'mice']) % 
%supertitle(['All areas combined' '  n=' num2str(sum(idx)) ])
supertitle(['areas  '  cell2mat(region(idx)') '  n=' num2str(sum(idx)) ])
%% plot fidget and early release rate
Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Fidget(idx,:);
exp = sum(idx);
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1])
text(i_region,0.07,['n=' num2str(exp)])
ylim([-0.1 0.1])

end 
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Fidget (LED - Ctrl)')

subplot(2,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA_trial(idx,:);
exp = sum(idx);
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1])
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-0.2 0.1])

end 
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta early release rate (LED - Ctrl)')

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Fidget(idx,:);
temp_p = power(idx,:);

scatter(temp_p,diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1])
hold on
 %ylim([-0.25 0.1])
xlim([0 1])
end 
set(gca,'XTick',0:0.2:1,'TickDir','out')
ylabel('\Delta Fidget (LED - Ctrl)')
xlabel('Led/laser Power(mW)')

subplot(2,2,4)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA_trial(idx,:);
temp_p = power(idx,:);

scatter(temp_p,diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1])
hold on
 ylim([-0.25 0.1])
xlim([0 1])
end 
set(gca,'XTick',0:0.2:1,'TickDir','out')
ylabel('\Delta early release rate (LED - Ctrl)')
xlabel('Led/laser Power(mW)')
supertitle('Percentage of fidget release (<3 cycle) and early release(>=3 cycle) of total trials')
%% plot the trial length dependence of only control condition
idx = logical(ones(length(IDs),1));
figure
temp = [];
temp = squeeze(Thresh_bintrial(idx,:,1));
temp = temp./repmat(temp(:,1),1,4); 
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0 0]},0)
hold on
temp = [];
temp = squeeze(FA_bintrial(idx,:,1));
temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0.6 0.4]},0)
% ylim([0 0.4])
xlim([0.5 5.5])
set(gca,'TickDir','out')
xlabel('Trial length (s)')
ylabel('Thresh/FA norm. by first bin')
title('Orientation change detection')


%% plot the overlayed effects of led on threshold and FA normalized by each bin of ctrl
figure
for i =1:4
    subplot(2,2,i)
idx =  strcmp(region,all_region(i));

temp1 = [];
temp1 = squeeze(Thresh_bintrial(idx,:,1));

temp2 = [];
temp2 = squeeze(Thresh_bintrial(idx,:,2));
temp = temp2-temp1; 
% yyaxis left
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0 0]},0)
ylim([-12 20])
% ylabel('Delta Thresh to Ctrl bins')
hold on

temp1 = [];
temp1 = squeeze(FA_bintrial(idx,:,1));

temp2 = [];
temp2 = squeeze(FA_bintrial(idx,:,2));
temp = temp2 - temp1; 
% yyaxis right
% ylim([-0.12 0.2])
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0.6 0.4]},0)
ylabel('Delta FA to Ctrl bins')

%  ylim([-0.05 0.2])
xlim([0.5 5.5])
title(['Area: ', all_region{i},'   n = ', num2str(sum(idx)), 'mice'])
set(gca,'TickDir','out')
xlabel('Trial length (s)')

end 
supertitle('Orientation change detection')

%% plot the trial length dependence of threshold 
figure
i = 4;
idx =  strcmp(region,all_region(i));
subplot(2,3,1)
temp = [];
temp = squeeze(Thresh_bintrial(idx,:,1));
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0 0]},0)
hold on
temp = [];
temp = squeeze(Thresh_bintrial(idx,:,2));
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0.6 1]},0)
ylim([0 45])
xlim([0.5 5.5])
set(gca,'TickDir','out')
xlabel('Trial length (s)')
ylabel('Threshhold (Deg)')

subplot(2,3,2)
temp = [];
temp = squeeze(Thresh_bintrial(idx,:,1));
temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0 0]},0)
hold on
temp = [];
temp = squeeze(Thresh_bintrial(idx,:,2));
temp = temp./repmat(temp(:,1),1,4); 
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0.6 1]},0)
ylim([0 1.1])
xlim([0.5 5.5 ])
set(gca,'TickDir','out')
xlabel('Trial length (s)')
ylabel('Norm. Thresh to each first bin')

subplot(2,3,3)
temp1 = [];
temp1 = squeeze(Thresh_bintrial(idx,:,1));
temp = temp1./temp1; % normalized by corresponding control condition
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0 0]},0)
hold on
temp2 = [];
temp2 = squeeze(Thresh_bintrial(idx,:,2));
temp = temp2./temp1; 
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0.6 1]},0)
ylim([0 2.5])
xlim([0.5 5.5 ])
set(gca,'TickDir','out')
xlabel('Trial length (s)')
ylabel('Norm. Thresh to Ctrl bins')


subplot(2,3,4)
temp = [];
temp = squeeze(FA_bintrial(idx,:,1));
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0 0]},0)
hold on
temp = [];
temp = squeeze(FA_bintrial(idx,:,2));
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0.6 1]},0)
ylim([0 0.25])
xlim([0.5 5.5])
set(gca,'TickDir','out')
xlabel('Trial length (s)')
ylabel('FA rate')

subplot(2,3,5)
temp = [];
temp = squeeze(FA_bintrial(idx,:,1));
temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0 0]},0)
hold on
temp = [];
temp = squeeze(FA_bintrial(idx,:,2));
temp = temp./repmat(temp(:,1),1,4); 
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0.6 1]},0)
ylim([0 3])
xlim([0.5 5.5 ])
set(gca,'TickDir','out')
xlabel('Trial length (s)')
ylabel('Norm. FA to each first bin')

subplot(2,3,6)
temp1 = [];
temp1 = squeeze(FA_bintrial(idx,:,1));
temp = temp1./temp1; % normalized by its control corresponding bin
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0 0]},0)
hold on
temp2 = [];
temp2 = squeeze(FA_bintrial(idx,:,2));
temp = temp2./temp1; 
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',[0 0.6 1]},0)
ylim([0 1.5])
xlim([0.5 5.5 ])
set(gca,'TickDir','out')
xlabel('Trial length (s)')
ylabel('Norm. FA to Ctrl bins')

supertitle(['LED area: ', all_region{i},'   n = ', num2str(sum(idx)), 'mice']) % 




%%  summary of LED effects
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
 text(i_region,19,['n=' num2str(exp)])
ylim([-5 20])

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
ylim([-0.1 0.1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
axis square

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = d_22(idx,:);
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
ylabel('\Delta d prime-22.5 (LED - Ctrl)')
axis square

subplot(2,2,4)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = c_22(idx,:);
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
ylabel('\Delta c-22.5 (LED - Ctrl)')
axis square
%% scatter plot of c22 and d22 for areas involved in orientation change detection
idx = ~strcmp(region,all_region(4));
figure
subplot(1,2,1)
scatter(d_22(idx,1), d_22(idx,2),'MarkerEdgeColor',[0 0 0])
hold on
errorbar(mean(d_22(idx,1)), mean(d_22(idx,2)),std(d_22(idx,1))./sqrt(sum(idx)),std(d_22(idx,1))./sqrt(sum(idx)),std(d_22(idx,2))./sqrt(sum(idx)),std(d_22(idx,2))./sqrt(sum(idx)),'marker','o')
hold on
plot([0 3],[0,3],'k:')
ylim([0 3])
xlim([0 3])
set(gca,'XTick',0:1:3, 'YTick', 0:1:3, 'TickDir','out')
xlabel('d prime (control)')
ylabel('d prime (led)')
axis square
subplot(1,2,2)
scatter(c_22(idx,1), c_22(idx,2),'MarkerEdgeColor',[0 0 0])
hold on
errorbar(mean(c_22(idx,1)), mean(c_22(idx,2)),std(c_22(idx,1))./sqrt(sum(idx)),std(c_22(idx,1))./sqrt(sum(idx)),std(c_22(idx,2))./sqrt(sum(idx)),std(c_22(idx,2))./sqrt(sum(idx)),'marker','o')
plot([0 2],[0,2],'k:')
ylim([-0.1 2])
xlim([-0.1 2])
set(gca,'XTick',0:0.5:2, 'YTick', 0:0.5:2, 'TickDir','out')
xlabel('criterion (control)')
ylabel('criterion (led)')
axis square

%% plot the immedate LED control


Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
all_region = {'V1','LM','AL','PM'};
figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = [LEDpre_thresh(idx,3) LEDpre_thresh(idx,1)];
[h_tp(i_region),p_tp(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
 text(i_region,18,['n=' num2str(exp)])
 ylim([-5 20])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Thresh (LED/Ctrl - Ctrl/Ctrl)')
axis square
subplot(2,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = [LEDpre_FA(idx,3) LEDpre_FA(idx,1)];

[h_fp(i_region),p_fp(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-0.1 0.1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta FA (LED/Ctrl - Ctrl/Ctrl)')
axis square

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = [LEDpre_d(idx,3) LEDpre_d(idx,1)];
[h_dp(i_region),p_dp(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)

  ylim([-1 1])

end 
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta d (LED/Ctrl - Ctrl/Ctrl)')
axis square

subplot(2,2,4)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = [LEDpre_c(idx,3) LEDpre_c(idx,1)];

[h_cp(i_region),p_cp(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-1 1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta c (LED/Ctrl - Ctrl/Ctrl)')
axis square
%% plot LED/LED against Ctrl/LED
figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = [LEDpre_thresh(idx,4) LEDpre_thresh(idx,2)];
[h_tp(i_region),p_tp(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
 text(i_region,18,['n=' num2str(exp)])
 ylim([-5 20])

end 
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Thresh (LED/LED - Ctrl/LED)')
axis square

subplot(2,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = [LEDpre_FA(idx,4) LEDpre_FA(idx,2)];

[h_fp(i_region),p_fp(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-0.1 0.1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta FA (LED/LED - Ctrl/LED)')
axis square

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = [LEDpre_d(idx,4) LEDpre_d(idx,2)];
[h_dp(i_region),p_dp(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)

  ylim([-1 1])

end 
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta d (LED/LED - Ctrl/LED)')
axis square

subplot(2,2,4)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = [LEDpre_c(idx,4) LEDpre_c(idx,2)];

[h_cp(i_region),p_cp(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-1 1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta c (LED/LED - Ctrl/LED)')
axis square
%% plot the immediate LED control
Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
all_region = {'V1','LM','AL','PM'};
figure
subplot(1,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
tempx = LEDpre_thresh(idx,3:4); % previouse control
tempy = LEDpre_thresh(idx,1:2); % previouse LED
exp = sum(idx);
scatter(tempx(:,1),tempy(:,1),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',50)
hold on
scatter(tempx(:,2),tempy(:,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',50)

end 
ylim([0 50])
xlim([0 50])
set(gca,'TickDir','out')
plot([0 50],[0 50],'k:')
xlabel('Threshold after Ctrl trial')
ylabel('Threshold after LED trial')

subplot(1,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
tempx = LEDpre_FA(idx,3:4); % previouse control
tempy = LEDpre_FA(idx,1:2); % previouse LED
exp = sum(idx);
scatter(tempx(:,1),tempy(:,1),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',50)
hold on
scatter(tempx(:,2),tempy(:,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',50)

end 
ylim([0 0.15])
xlim([0 0.15])
set(gca,'TickDir','out')
plot([0 0.15],[0 0.15],'k:')
xlabel('FA rate after Ctrl trial')
ylabel('FA rate after LED trial')


%% plot the RT for 90 deg vs 22.5 deg 20<x<30
for i = 1:size(Oriens,1)
    for j=1:2
        temp =squeeze(RT(i,:,j));
        temp_orien = squeeze(Oriens(i,:,j));
      RT_22(i,j) =  mean(temp(temp_orien>20 & temp_orien<30));
      RT_90(i,j) = temp(temp_orien==max(temp_orien));
      RT_32(i,j) =  mean(temp(temp_orien>30 & temp_orien<=45));
    end 
end 

Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
all_region = {'V1','LM','AL','PM'};
figure
subplot(1,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = RT_22(idx,:);
[h_22(i_region),p_22(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
%  text(i_region,80,['n=' num2str(exp)])
 ylim([-20 80])

end 
hline(0,'k:')
xlim([0 5])
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta RT22 (LED - Ctrl)')

subplot(1,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = RT_90(idx,:);
[h_90(i_region),p_90(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-20 80])

end 
hline(0,'k:')
xlim([0 5])
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta RT90 (LED - Ctrl)')

%% plot the different between RT22 and 90
Diff_RT = RT_32-RT_90; 
figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Diff_RT(idx,:);
[h_22(i_region),p_22(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
%  text(i_region,80,['n=' num2str(exp)])
 ylim([-20 80])

end 
hline(0,'k:')
xlim([0 5])
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('Diff(RT22.5 - RT32)(LED-Ctrl)')
subplot(2,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Diff_RT(idx,1);
[h_22(i_region),p_22(i_region)]= ttest(RT_22(:,1),RT_90(:,1));
exp = sum(idx);
scatter(repmat(i_region,exp,1),temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',10)
hold on
errorbar(i_region,mean(temp),std(temp)./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(temp),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
%  text(i_region,80,['n=' num2str(exp)])
 ylim([-20 80])

end 
hline(0,'k:')
xlim([0 5])
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('(RT22.5 - RT32)Ctrl')

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Diff_RT(idx,2);
[h_22(i_region),p_22(i_region)]= ttest(RT_22(:,1),RT_90(:,1));
exp = sum(idx);
scatter(repmat(i_region,exp,1),temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',10)
hold on
errorbar(i_region,mean(temp),std(temp)./sqrt(exp),'Color',Colors(i_region,:))
hold on
scatter(i_region,mean(temp),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
%  text(i_region,80,['n=' num2str(exp)])
 ylim([-20 80])

end 
hline(0,'k:')
xlim([0 5])
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('(RT22.5 - RT32)LED')

%% scatter plot of power vs threshold change
figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Thresh(idx,:);
temp_p = power(idx,:);
scatter(temp_p,diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1])
hold on
 ylim([-5 20])
xlim([0 1])
end 
set(gca,'XTick',0:0.2:1,'TickDir','out')
xlabel('Led/laser Power(mW)')
ylabel('\Delta Thresh (LED - Ctrl)')

subplot(2,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA(idx,:);
temp_p = power(idx,:);

scatter(temp_p,diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1])
hold on
 ylim([-0.1 0.1])
xlim([0 1])
end 
set(gca,'XTick',0:0.2:1,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
xlabel('Led/laser Power(mW)')

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = d_22(idx,:);
temp_p = power(idx,:);

scatter(temp_p,diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1])
hold on
ylim([-1 1])
xlim([0 1])
end 
set(gca,'XTick',0:0.2:1,'TickDir','out')
ylabel('\Delta d prime-22.5 (LED - Ctrl)')
xlabel('Led/laser Power(mW)')


subplot(2,2,4)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = c_22(idx,:);
temp_p = power(idx,:);

scatter(temp_p,diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1])
hold on
 ylim([-0.5 1])
xlim([0 1])
end 
set(gca,'XTick',0:0.2:1,'TickDir','out')
ylabel('\Delta c-22.5 (LED - Ctrl)')
xlabel('Led/laser Power(mW)')

%% plot delta thresh over delta FA seperate by area
all_region = {'V1','LM','AL','PM'};
Colors(1,:) = [0 0 0];
Colors(2:5,:)= lines(4); % color coded is the virus
shapes = {'square','o','pentagram'}; % vgat, pv cre, gad cre
Virus_unique = unique(virus);
back_unique = unique(background);

figure
for i_region = 1: 4  % four regions
    subplot(2,2,i_region)
    idx_region = strcmp(region,all_region(i_region));
    temp_FA = diff(FA(idx_region,:),1,2);
    temp_thresh = diff(Thresh(idx_region,:),1,2);
    temp_back = background(idx_region);
    temp_virus = virus(idx_region);
    temp_ID = IDs(idx_region);
    for i = 1:sum(idx_region)
        shape_idx = find(strcmp(back_unique,temp_back(i))==1);
        color_idx = find(strcmp(Virus_unique,temp_virus(i))==1);
        h=scatter(temp_FA(i),temp_thresh(i),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(color_idx,:),'SizeData',60);
        h.Marker = char(shapes(shape_idx));
        hold on
        text(temp_FA(i), -4,num2str(temp_ID(i)))
    end 
    ylim([-5 15])
    xlim([-0.1 0.1])
    vline(0,'k:')
    hline(0,'k:')
    title(['Area: ' char(all_region(i_region))])
    ylabel('\Delta Thresh (Deg)')
    xlabel('\Delta FA rate')
    
    if i_region ==1
        for i_shape = 1:length(shapes)
        h = scatter(0.04,15-1.2*i_shape,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
        h.Marker = char(shapes(i_shape));
        text(0.05,15-1.2*i_shape,back_unique(i_shape),'k')
        end
        
        for i_color = 1:size(Colors,1)
        
        text(-0.09,5-1.2*i_color,Virus_unique(i_color),'Color',Colors(i_color,:))    
            
        end 
    end 
    
end 
supertitle('Orientation change detection')



%% plot delta FA over FA 
Colors(1,:) = [0 0 0];
Colors(2:5,:)= lines(4); % color coded is the virus
shapes = {'square','o','pentagram'}; % vgat, pv cre, gad cre
Virus_unique = unique(virus);
back_unique = unique(background);

figure
for i_region = 1: 4  % four regions
    subplot(2,2,i_region)
    idx_region = strcmp(region,all_region(i_region));
    temp_thresh = diff(FA(idx_region,:),1,2);
    temp_FA=FA(idx_region,1);
    temp_back = background(idx_region);
    temp_virus = virus(idx_region);
    for i = 1:sum(idx_region)
        shape_idx = find(strcmp(back_unique,temp_back(i))==1);
        color_idx = find(strcmp(Virus_unique,temp_virus(i))==1);
        h=scatter(temp_FA(i),temp_thresh(i),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(color_idx,:),'SizeData',60);
        h.Marker = char(shapes(shape_idx));
        hold on
    end
    ylim([-0.1 0.1])
    xlim([0 0.3])
    vline(0,'k:')
    hline(0,'k:')
    title(['Area: ' char(all_region(i_region))])
    ylabel('\Delta FA rate')
    xlabel('FA rate')
    
    if i_region ==1
        for i_shape = 1:length(shapes)
        h = scatter(0.2,0.05-0.012*i_shape,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
        h.Marker = char(shapes(i_shape));
        text(0.22,0.05-0.012*i_shape,back_unique(i_shape),'k')
        end
        
        for i_color = 1:size(Colors,1)
        
        text(0,0.1-0.012*i_color,Virus_unique(i_color),'Color',Colors(i_color,:))    
            
        end 
    end 
    
end 
supertitle('Orientation change detection')
%% plot delta Thresh over Thresh
Colors(1,:) = [0 0 0];
Colors(2:5,:)= lines(4); % color coded is the virus
shapes = {'square','o','pentagram'}; % vgat, pv cre, gad cre
Virus_unique = unique(virus);
back_unique = unique(background);

figure
for i_region = 1: 4  % four regions
    subplot(2,2,i_region)
    idx_region = strcmp(region,all_region(i_region));
    temp_thresh = diff(Thresh(idx_region,:),1,2);
    temp_FA=Thresh(idx_region,1);
    temp_back = background(idx_region);
    temp_virus = virus(idx_region);
    for i = 1:sum(idx_region)
        shape_idx = find(strcmp(back_unique,temp_back(i))==1);
        color_idx = find(strcmp(Virus_unique,temp_virus(i))==1);
        h=scatter(temp_FA(i),temp_thresh(i),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(color_idx,:),'SizeData',60);
        h.Marker = char(shapes(shape_idx));
        hold on
    end
    ylim([-5 15])
    xlim([0 40])
    vline(0,'k:')
    hline(0,'k:')
    title(['Area: ' char(all_region(i_region))])
    ylabel('\Delta Thresh (Deg)')
    xlabel('Threshold (Deg)')
    
    if i_region ==4
        for i_shape = 1:length(shapes)
        h = scatter(23,10-1.2*i_shape,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
        h.Marker = char(shapes(i_shape));
        text(25,10-1.2*i_shape,back_unique(i_shape),'k')
        end
        
        for i_color = 1:size(Colors,1)
        
        text(0,15-1.2*i_color,Virus_unique(i_color),'Color',Colors(i_color,:))    
            
        end 
    end 
    
end 
supertitle('Orientation change detection')

%% get the slope the FA trial length and threshold trial length
for i = 1:size(Thresh_bintrial,1)
Thresh_slope(i,1:2) = polyfit(cbin./1000,squeeze(Thresh_bintrial(i,:,1)),1);
FA_slope (i,1:2) = polyfit(cbin./1000,squeeze(FA_bintrial(i,:,1)),1);
end 

figure
for i_region = 1: 4  % four regions
    subplot(2,2,i_region)
    idx_region = strcmp(region,all_region(i_region));
    temp_thresh = Thresh_slope(idx_region,1);
    temp_FA=FA_slope(idx_region,1);
    temp_back = background(idx_region);
    temp_virus = virus(idx_region);
    for i = 1:sum(idx_region)
        shape_idx = find(strcmp(back_unique,temp_back(i))==1);
        color_idx = find(strcmp(Virus_unique,temp_virus(i))==1);
        h=scatter(temp_FA(i),temp_thresh(i),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(color_idx,:),'SizeData',60);
        h.Marker = char(shapes(shape_idx));
        hold on
    end
     ylim([-3 3])
     xlim([-0.02 0.06])
    vline(0,'k:')
    hline(0,'k:')
    title(['Area: ' char(all_region(i_region))])
    ylabel('Ctrl Thresh slope (Deg/s)')
    xlabel('Ctrl FA slope( /s)')
    
%     if i_region ==4
%         for i_shape = 1:length(shapes)
%         h = scatter(23,10-1.2*i_shape,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
%         h.Marker = char(shapes(i_shape));
%         text(25,10-1.2*i_shape,back_unique(i_shape),'k')
%         end
%         
%         for i_color = 1:size(Colors,1)
%         
%         text(0,15-1.2*i_color,Virus_unique(i_color),'Color',Colors(i_color,:))    
%             
%         end 
%     end 
    
end 
supertitle('Orientation change detection: trial length dependence')

%% get the slope the delta FA over slope FA
for i = 1:size(Thresh_bintrial,1)
Thresh_slope(i,1:2) = polyfit(cbin./1000,squeeze(Thresh_bintrial(i,:,1)),1);
FA_slope (i,1:2) = polyfit(cbin./1000,squeeze(FA_bintrial(i,:,1)),1);
end 

figure
for i_region = 1: 4  % four regions
    subplot(2,2,i_region)
    idx_region = strcmp(region,all_region(i_region));
    temp_thresh = diff(FA(idx_region,:),1,2);
    temp_FA=FA_slope(idx_region,1);
    temp_back = background(idx_region);
    temp_virus = virus(idx_region);
    for i = 1:sum(idx_region)
        shape_idx = find(strcmp(back_unique,temp_back(i))==1);
        color_idx = find(strcmp(Virus_unique,temp_virus(i))==1);
        h=scatter(temp_FA(i),temp_thresh(i),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(color_idx,:),'SizeData',60);
        h.Marker = char(shapes(shape_idx));
        hold on
    end
     ylim([-0.1 0.1])
     xlim([-0.02 0.06])
    vline(0,'k:')
    hline(0,'k:')
    title(['Area: ' char(all_region(i_region))])
    ylabel('\Delta FA')
    xlabel('Ctrl FA slope( /s)')
    
%     if i_region ==4
%         for i_shape = 1:length(shapes)
%         h = scatter(23,10-1.2*i_shape,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
%         h.Marker = char(shapes(i_shape));
%         text(25,10-1.2*i_shape,back_unique(i_shape),'k')
%         end
%         
%         for i_color = 1:size(Colors,1)
%         
%         text(0,15-1.2*i_color,Virus_unique(i_color),'Color',Colors(i_color,:))    
%             
%         end 
%     end 
    
end 
supertitle('Orientation change detection')
%% get the slope the delta Thresh over slope FA
for i = 1:size(Thresh_bintrial,1)
Thresh_slope(i,1:2) = polyfit(cbin./1000,squeeze(Thresh_bintrial(i,:,1)),1);
FA_slope (i,1:2) = polyfit(cbin./1000,squeeze(FA_bintrial(i,:,1)),1);
end 

figure
for i_region = 1: 4  % four regions
    subplot(2,2,i_region)
    idx_region = strcmp(region,all_region(i_region));
    temp_thresh = diff(Thresh(idx_region,:),1,2);
    temp_FA=Thresh_slope(idx_region,1);
    temp_back = background(idx_region);
    temp_virus = virus(idx_region);
    for i = 1:sum(idx_region)
        shape_idx = find(strcmp(back_unique,temp_back(i))==1);
        color_idx = find(strcmp(Virus_unique,temp_virus(i))==1);
        h=scatter(temp_FA(i),temp_thresh(i),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(color_idx,:),'SizeData',60);
        h.Marker = char(shapes(shape_idx));
        hold on
    end
     ylim([-5 15])
     xlim([-3 3])
    vline(0,'k:')
    hline(0,'k:')
    title(['Area: ' char(all_region(i_region))])
    ylabel('\Delta Thresh (Deg)')
    xlabel('Ctrl Thresh slope( Deg /s)')
    
%     if i_region ==4
%         for i_shape = 1:length(shapes)
%         h = scatter(23,10-1.2*i_shape,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
%         h.Marker = char(shapes(i_shape));
%         text(25,10-1.2*i_shape,back_unique(i_shape),'k')
%         end
%         
%         for i_color = 1:size(Colors,1)
%         
%         text(0,15-1.2*i_color,Virus_unique(i_color),'Color',Colors(i_color,:))    
%             
%         end 
%     end 
    
end 
supertitle('Orientation change detection')

%% plot all areas together 
for i = 1:size(Thresh_bintrial,1)
Thresh_slope(i,1:2) = polyfit(cbin./1000,squeeze(Thresh_bintrial(i,:,1)),1);
FA_slope (i,1:2) = polyfit(cbin./1000,squeeze(FA_bintrial(i,:,1)),1);
end 
Colors = [];
Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
figure
subplot(2,2,1)
for i_region = 1: 4  
    idx_region = strcmp(region,all_region(i_region));
    temp_FA = diff(FA(idx_region,:),1,2);
    temp_thresh = diff(Thresh(idx_region,:),1,2);
    scatter(temp_FA, temp_thresh,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(i_region,:))
    hold on
    text(0.05,15-2.*i_region,all_region(i_region),'Color',Colors(i_region,:))
end 
ylim([-5 15])
xlim([-0.1 0.1])
vline(0,'k:')
hline(0,'k:')
ylabel('\Delta Thresh (Deg)')
xlabel('\Delta FA')
% [R,P] = corrcoef(diff(FA(~idx_region,:),1,2),diff(Thresh(~idx_region,:),1,2))

subplot(2,2,2)
for i_region = 1: 4  
    idx_region = strcmp(region,all_region(i_region));
    temp_FA = FA(idx_region,1);
    temp_thresh =  diff(FA(idx_region,:),1,2);
    scatter(temp_FA, temp_thresh,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(i_region,:))
    hold on
%     text(0.05,15-2.*i_region,all_region(i_region),'Color',Colors(i_region,:))
end 
ylim([-0.1 0.1])
xlim([0 0.2])
vline(0,'k:')
hline(0,'k:')
ylabel('\Delta FA')
xlabel(' FA')
% [R,P] = corrcoef(FA(:,1),diff(FA,1,2));

subplot(2,2,3)
for i_region = 1: 4  
    idx_region = strcmp(region,all_region(i_region));
    temp_FA = Thresh(idx_region,1);
    temp_thresh =  diff(Thresh(idx_region,:),1,2);
    scatter(temp_FA, temp_thresh,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(i_region,:))
    hold on
%     text(0.05,15-2.*i_region,all_region(i_region),'Color',Colors(i_region,:))
end 
ylim([-5 15])
xlim([0 35])
vline(0,'k:')
hline(0,'k:')
ylabel('\Delta Thresh(Deg)')
xlabel('Threshold(Deg)')
% [R,P] = corrcoef(Thresh(:,1),diff(Thresh,1,2));


subplot(2,2,4)
for i_region = 1: 4  
    idx_region = strcmp(region,all_region(i_region));
    temp_FA = FA_slope(idx_region,1);
    temp_thresh =  diff(FA(idx_region,:),1,2);
    scatter(temp_FA, temp_thresh,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Colors(i_region,:))
    hold on
%     text(0.05,15-2.*i_region,all_region(i_region),'Color',Colors(i_region,:))
end 
ylim([-0.1 0.1])
xlim([-0.02 0.06])
vline(0,'k:')
hline(0,'k:')
ylabel('\Delta FA')
xlabel('FA slope over trial length (/s)')
%  [R,P] = corrcoef(FA_slope(:,1),diff(FA,1,2));
supertitle('Orientation change detection: All areas combined')


%% plot the control threshold and FA over trial length
figure
subplot(1,2,1)
temp = squeeze(Thresh_bintrial(:,:,1));
temp = temp./repmat(temp(:,1),1,4); 
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(size(Thresh_bintrial,1)),{'Color',[0 0 0]},1)
ylim([0 1])
xlim([0.5 5.5])

hold on
temp = squeeze(FA_bintrial(:,:,1));
temp = temp./repmat(temp(:,1),1,4); 
shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(size(Thresh_bintrial,1)),{'Color',[0.1 0.8 0.5]},1)
ylim([0 2.5])
xlim([0.5 5.5])
legend ('Threshold', 'FA')
xlabel('Trial length (s)')
ylabel('Norm. FA/Thresh')
title('Orientation change detection')

