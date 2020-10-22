%% from i563 170411-170414
% modify the code that can load through multiple days 
% need to deal with the the same orientation with different values
 input = inputcombined;
T_Orien = cell2mat(input.tGratingDirectionDeg);
T_Orien = round(mod(T_Orien+180,180).*10,0)./10;
TT_Orien = T_Orien;
TT_Orien(TT_Orien>90) = 180 - TT_Orien(TT_Orien>90); % fold it 



B_Orien =double( cell2mat(input.tBaseGratingDirectionDeg));
B_Orien = mod(B_Orien+180,180);
B_ori = unique(B_Orien);
D_Orien =round((T_Orien - B_Orien).*10,0)./10;
D_Orien = round(mod(D_Orien+180,180).*10,0)./10;
D_ori = unique(D_Orien);
Success = strcmp(input.trialOutcomeCell,'success');
Miss =  strcmp(input.trialOutcomeCell,'ignore');
Early = strcmp(input.trialOutcomeCell,'failure');
ReactTime = cell2mat(input.reactTimesMs);
HoldTime = cell2mat(input.holdTimesMs);
tCycleNum = double(cell2mat(input.tCyclesOn));
Leverdown = double(cell2mat(input.tLeverPressTimeMs));
Leverup = double(cell2mat(input.tLeverReleaseTimeMs)); 
idx_HM = find(Early==0); % index of Hit and Miss trials
% correct the new_RT on hit and misese
    new_RT = cellfun(@minus,input.tLeverReleaseTimeMs,input.targetStimOnMs, 'UniformOutput', false);
    
    
    RT_HM = zeros(1,length(ReactTime));
    RT_HM(~cellfun('isclass',new_RT, 'double')) = cell2mat(new_RT(~cellfun('isclass',new_RT, 'double')));
    
    RT_HM (cellfun('isclass',new_RT, 'double')) = cell2mat(new_RT(cellfun('isclass',new_RT, 'double')));
    
    % adjust the Success, Miss trial conditions
    % throw too fast response into early trials
    trialnum = length(Success);
    SuccessN = zeros(1,trialnum);
    MissN = zeros(1,trialnum);
    EarlyN = Early;
    StimOff = NaN(1,trialnum);
    StimOff_pre = NaN(1,trialnum);
    for i=1:length(idx_HM)
        idx_temp = [];
        idx_temp = idx_HM(i);% trial of hit or miss trial
        cycle = [];
        trialoffs= [];
        cycle = tCycleNum(idx_temp);
        trialoffs = double(input.tStimOffTimes{idx_temp});
        StimOff(idx_temp) = trialoffs(cycle);% interval before target
        StimOff_pre(idx_temp) = trialoffs(cycle-1); % interval before last cycle
        
%         RT_temp = [];        
%         RT_temp = Leverup(idx_temp) - Leverdown(idx_temp) - 100.*cycle -sum(trialoffs(1:cycle));
%         RT_HM(idx_temp) = RT_temp;
        
        if RT_HM(idx_temp)>=200 && RT_HM(idx_temp)<=550
            SuccessN(idx_temp) = 1;   % new success trials
        end
        if RT_HM(idx_temp)<200
            EarlyN(idx_temp) = 1;   % add too fast time into early trials
        end
        if RT_HM(idx_temp)>550
            MissN(idx_temp) = 1;   % new miss trials
        end
        
        
    end

Off = unique(StimOff(idx_HM));

% check if any base orientation impose more difficulties 
for i_orien = 1:length(B_ori)
    
    [output.target{i_orien}] = ISI_HR(StimOff(B_Orien==B_ori(i_orien)),Off,D_ori,SuccessN(B_Orien==B_ori(i_orien)),MissN(B_Orien==B_ori(i_orien)),tCycleNum(B_Orien==B_ori(i_orien)),D_Orien(B_Orien==B_ori(i_orien)),RT_HM(B_Orien==B_ori(i_orien)));
    [output.pre{i_orien}] = ISI_HR(StimOff_pre(B_Orien==B_ori(i_orien)),Off,D_ori,SuccessN(B_Orien==B_ori(i_orien)),MissN(B_Orien==B_ori(i_orien)),tCycleNum(B_Orien==B_ori(i_orien)),D_Orien(B_Orien==B_ori(i_orien)),RT_HM(B_Orien==B_ori(i_orien)));
    % FA rate for different base orientations 
    Input = {};
    Input.stimOnTimeMs =  input.stimOnTimeMs;
    Input.tStimOffTimes = input.tStimOffTimes(B_Orien==B_ori(i_orien));
   
    [output.FA{i_orien}] =  ISI_FA_N(EarlyN(B_Orien==B_ori(i_orien)),Input,Off,tCycleNum(B_Orien==B_ori(i_orien)),sum(B_Orien==B_ori(i_orien)),Leverup(B_Orien==B_ori(i_orien)),Leverdown(B_Orien==B_ori(i_orien)));
    %     for i_d = 1:length(D_ori)
%         Hit_num(i_orien,i_d) = sum(D_Orien==D_ori(i_d)&B_Orien==B_ori(i_orien) & SuccessN);
%         Miss_num(i_orien,i_d) = sum(D_Orien==D_ori(i_d)&B_Orien==B_ori(i_orien) & MissN);
%         Hit(i_orien,i_d) = Hit_num(i_orien,i_d) ./(Hit_num(i_orien,i_d) + Miss_num(i_orien,i_d));
%         [~,confi(i_orien,i_d,1:2)]= binofit(Hit_num(i_orien,i_d),(Hit_num(i_orien,i_d)+ Miss_num(i_orien,i_d)));
%         Hit_RT{i_orien,i_d} = ReactTime(D_Orien==D_ori(i_d)&B_Orien==B_ori(i_orien) & SuccessN);
%         
%     end
%     
    
end

% break up by absolute target value to calculate hit rate
for i_orien = 1:length(B_ori)
    
    T_orien(i_orien,:) = unique(TT_Orien(B_Orien==B_ori(i_orien)));
    
    [output.Abstarget{i_orien}] = ISI_HR(StimOff(B_Orien==B_ori(i_orien)),Off,T_orien(i_orien,:),SuccessN(B_Orien==B_ori(i_orien)),MissN(B_Orien==B_ori(i_orien)),tCycleNum(B_Orien==B_ori(i_orien)),TT_Orien(B_Orien==B_ori(i_orien)),RT_HM(B_Orien==B_ori(i_orien)));
    [output.Abspre{i_orien}] = ISI_HR(StimOff_pre(B_Orien==B_ori(i_orien)),Off,T_orien(i_orien,:),SuccessN(B_Orien==B_ori(i_orien)),MissN(B_Orien==B_ori(i_orien)),tCycleNum(B_Orien==B_ori(i_orien)),TT_Orien(B_Orien==B_ori(i_orien)),RT_HM(B_Orien==B_ori(i_orien)));
    
% for i_t = 1:length(T_orien(i_orien,:))
%    
%     Hit_numT(i_orien,i_t) = sum(TT_Orien==T_orien(i_orien,i_t) &B_Orien==B_ori(i_orien)&SuccessN);
%     Miss_numT(i_orien,i_t) = sum(TT_Orien==T_orien(i_orien,i_t) &B_Orien==B_ori(i_orien)& MissN);
%     Hit_T(i_orien,i_t) = Hit_numT(i_orien,i_t)./(Hit_numT(i_orien,i_t)+Miss_numT(i_orien,i_t));
%     [~,confi_T(i_orien,i_t,1:2)]= binofit(Hit_numT(i_orien,i_t),(Hit_numT(i_orien,i_t)+ Miss_numT(i_orien,i_t)));
%     Hit_TRT{i_orien,i_t} =  ReactTime(TT_Orien==T_orien(i_orien,i_t) & B_Orien==B_ori(i_orien) & SuccessN);
% end
end 

% caclulate FA rate base on the absolute base orien
for i_orien = 1:length(B_ori)
    Early_B(i_orien) = sum(B_Orien==B_ori(i_orien) & EarlyN);
    Hit_numB(i_orien) = sum(B_Orien==B_ori(i_orien) & SuccessN);
    Miss_numB(i_orien) = sum(B_Orien==B_ori(i_orien) & MissN);
    Hit_B(i_orien) = Hit_numB(i_orien)./(Hit_numB(i_orien)+ Miss_numB(i_orien));
    Hit_BRT{i_orien} = ReactTime(B_Orien==B_ori(i_orien) & SuccessN);
    FA_B(i_orien) = Early_B(i_orien)./sum(B_Orien==B_ori(i_orien));
    FA_BRT{i_orien} = HoldTime (B_Orien==B_ori(i_orien) & EarlyN);
end
% get bootstrap fitting for each conditions
fitS = {};
AbsfitS = {};
fitSall = {};
AbsfitSall = {};
Threshall = [];
AbsThreshall = [];
fitPre = {};
AbsfitPre = {};
Thresh = [];
Threshpre = [];
AbsThresh = [];
AbsThreshpre = [];
bootStats={};
AbsbootStats={};
bootStatsall={};
AbsbootStatsall={};
bootStatspre = {};
AbsbootStatspre = {};
NBootstrapReps = 1000;
Thresh_ci = [];
Threshpre_ci = [];
Threshall_ci = [];
AbsThresh_ci = [];
AbsThreshpre_ci = [];
AbsThreshall_ci = [];
for i = 1:length(B_ori)
    for i_off = 1:length(Off)
        % target diffrence
        trialVec = (output.target{i}.HT_num{i_off,1}+output.target{i}.Miss_num{i_off,1})';
        trialVpre =(output.pre{i}.HT_num{i_off,1}+output.pre{i}.Miss_num{i_off,1})';
        fitS{i,1}{i_off,1} = weibullFitLG(D_ori, output.target{i}.c_hit{i_off,1}',1, 1, {'nTrials', trialVec}); % use
        fitPre{i,1}{i_off,1} = weibullFitLG(D_ori, output.pre{i}.c_hit{i_off,1}',1, 1, {'nTrials', trialVpre});
        Thresh(i,i_off) = fitS{i,1}{i_off,1}.thresh;
        Threshpre(i,i_off) = fitPre{i,1}{i_off,1}.thresh;
        
        
        [bootStats{i,1}{i_off,1}]=BootstrapWeibullFit(trialVec,output.target{i}.c_hit{i_off,1}',NBootstrapReps,D_ori,1, 1);
        [bootStatspre{i,1}{i_off,1}]=BootstrapWeibullFit(trialVpre,output.pre{i}.c_hit{i_off,1}',NBootstrapReps,D_ori,1, 1);
        Thresh_ci(i,i_off,1:2) =  bootStats{i,1}{i_off,1}.ci95;
        Threshpre_ci(i,i_off,1:2) =  bootStatspre{i,1}{i_off,1}.ci95;
        
        % absolute target
        trialVec = (output.Abstarget{i}.HT_num{i_off,1}+output.Abstarget{i}.Miss_num{i_off,1})';
        trialVpre =(output.Abspre{i}.HT_num{i_off,1}+output.Abspre{i}.Miss_num{i_off,1})';
        AbsfitS{i,1}{i_off,1} = weibullFitLG(T_orien(i,:), output.Abstarget{i}.c_hit{i_off,1}',1, 1, {'nTrials', trialVec}); % use
        AbsfitPre{i,1}{i_off,1} = weibullFitLG(T_orien(i,:), output.Abspre{i}.c_hit{i_off,1}',1, 1, {'nTrials', trialVpre});
        AbsThresh(i,i_off) = AbsfitS{i,1}{i_off,1}.thresh;
        AbsThreshpre(i,i_off) = AbsfitPre{i,1}{i_off,1}.thresh;
        
        
        [AbsbootStats{i,1}{i_off,1}]=BootstrapWeibullFit(trialVec,output.Abstarget{i}.c_hit{i_off,1}',NBootstrapReps,T_orien(i,:),1, 1);
        [AbsbootStatspre{i,1}{i_off,1}]=BootstrapWeibullFit(trialVpre,output.Abspre{i}.c_hit{i_off,1}',NBootstrapReps,T_orien(i,:),1, 1);
        AbsThresh_ci(i,i_off,1:2) =  AbsbootStats{i,1}{i_off,1}.ci95;
        AbsThreshpre_ci(i,i_off,1:2) =  AbsbootStatspre{i,1}{i_off,1}.ci95;
        
        
        
    end
    % collapsed all the offs for target difference
    Hits = output.target{i}.all.c_hit;
    trialAll =  output.target{i}.all.HT_num + output.target{i}.all.Miss_num;
    fitSall{i,1} = weibullFitLG(D_ori, Hits',1, 1, {'nTrials', trialAll'});
    Threshall(i) = fitSall{i,1}.thresh;
    
    [bootStatsall{i,1}]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,D_ori,1, 1);
    Threshall_ci(i,1:2) = bootStatsall{i,1}.ci95;
    % collapsed all offs for absolute target
    Hits = output.Abstarget{i}.all.c_hit;
    trialAll =  output.Abstarget{i}.all.HT_num + output.Abstarget{i}.all.Miss_num;
    AbsfitSall{i,1} = weibullFitLG(T_orien(i,:), Hits',1, 1, {'nTrials', trialAll'});
    AbsThreshall(i) = AbsfitSall{i,1}.thresh;
    
    [AbsbootStatsall{i,1}]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,T_orien(i,:),1, 1);
    AbsThreshall_ci(i,1:2) = AbsbootStatsall{i,1}.ci95;
    
    
end
%% get the first 1/3 and last 1/3 of training in seperate script
% called randombase_learning
chunk = round(trialnum/3);
[first] = randombase_learning(B_ori,StimOff(1:chunk),B_Orien(1:chunk),D_ori,Off,SuccessN(1:chunk),MissN(1:chunk),tCycleNum(1:chunk),D_Orien(1:chunk),RT_HM(1:chunk),EarlyN(1:chunk),...
    Leverup(1:chunk),Leverdown(1:chunk),TT_Orien(1:chunk),input.tStimOffTimes(1:chunk));

[last] = randombase_learning(B_ori,StimOff(trialnum-chunk+1:trialnum),B_Orien(trialnum-chunk+1:trialnum),D_ori,Off,SuccessN(trialnum-chunk+1:trialnum),MissN(trialnum-chunk+1:trialnum),...
    tCycleNum(trialnum-chunk+1:trialnum),D_Orien(trialnum-chunk+1:trialnum),RT_HM(trialnum-chunk+1:trialnum),EarlyN(trialnum-chunk+1:trialnum),Leverup(trialnum-chunk+1:trialnum),...
    Leverdown(trialnum-chunk+1:trialnum),TT_Orien(trialnum-chunk+1:trialnum),input.tStimOffTimes(trialnum-chunk+1:trialnum));
%% plot the hit rate color coded by different orientations
colors = lines(length(B_ori));
  figure
subplot(1,2,1)

for i_orien = 1:length(B_ori)
    plot(D_ori,output.target{i_orien}.all.HT_rate,'Color',colors(i_orien,:))
    hold on
    scatter(D_ori,output.target{i_orien}.all.HT_rate,'MarkerEdgeColor',colors(i_orien,:))
    for i = 1:length(D_ori)
    plot([D_ori(i),D_ori(i)],output.target{i_orien}.all.confi(i,:),'Color',colors(i_orien,:))
    end 
    text(11,1-0.05*i_orien,[num2str(B_ori(i_orien)) 'Deg'],'Color',colors(i_orien,:))
end


ylim([0 1])
xlim([10 100])
ylabel('Hit rate')
xlabel('\Delta Orientation')
set(gca,'XScale','log','TickDir','out')

subplot(1,2,2)


for i_orien = 1:length(B_ori)
    plot(T_orien(i_orien,:),output.Abstarget{i_orien}.all.HT_rate ,'Color',colors(i_orien,:))
    hold on
    scatter(T_orien(i_orien,:),output.Abstarget{i_orien}.all.HT_rate ,'MarkerEdgeColor',colors(i_orien,:))
    ori_temp = T_orien(i_orien,:);
    for i = 1:length(ori_temp)
    plot([ori_temp(i),ori_temp(i)],output.Abstarget{i_orien}.all.confi(i,:),'Color',colors(i_orien,:))
    end
end

ylim([0 1])
xlim([5 100])
ylabel('Hit rate')
xlabel('Absolute target Orientation')
set(gca,'XScale','log','TickDir','out')
RealFA_B(1) = output.FA{1}.all.FA;
RealFA_B(2) = output.FA{2}.all.FA;
RealFA_B(3) = output.FA{3}.all.FA;
%supertitle([num2str(input.subjectNum) ' - ' input.saveTime(1:6) ' - '  'base orien:' num2str(B_ori) 'Deg' '   Earlys = ' num2str(FA_B)])
supertitle([num2str(input.ID) ' - '  'base orien:' num2str(B_ori) 'Deg' '   Earlys = ' num2str(RealFA_B)])

%% plot the FA rate/threshold separate by differnt offs 
figure
subplot(3,2,1)
for i_orien = 1:length(B_ori)
    for i_off  = 1:length(Off)
        plot([i_off i_off], squeeze(Thresh_ci(i_orien,i_off,:)),'Color',colors(i_orien,:))
        hold on
       
    end 
   
    scatter([1 2 3], Thresh(i_orien,:), 'MarkerEdgeColor',colors(i_orien,:))  
end 
ylim([0 50])
xlim([0.5 3.5])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Thresh delta orien (deg)')
% normalized to its own
subplot(3,2,2)
for i_orien = 1:length(B_ori)
   
    scatter([1 2 3], Thresh(i_orien,:)./Thresh(i_orien,1), 'MarkerEdgeColor',colors(i_orien,:))  
    hold on
    text(0.6,1-0.15*i_orien,[num2str(B_ori(i_orien)) 'Deg'],'Color',colors(i_orien,:))
end 
ylim([0 1])
xlim([0.5 3.5])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Thresh(Deg)')

subplot(3,2,3)
for i_orien = 1:length(B_ori)
    for i_off  = 1:length(Off)
        plot([i_off i_off], squeeze(AbsThresh_ci(i_orien,i_off,:)),'Color',colors(i_orien,:))
        hold on
       
    end 
   
    scatter([1 2 3], AbsThresh(i_orien,:), 'MarkerEdgeColor',colors(i_orien,:))  
end 
ylim([0 50])
xlim([0.5 3.5])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Thresh abs orien (deg)')
% normalized to its own
subplot(3,2,4)
for i_orien = 1:length(B_ori)
   
    scatter([1 2 3], AbsThresh(i_orien,:)./AbsThresh(i_orien,1), 'MarkerEdgeColor',colors(i_orien,:))  
    hold on
    
end 
ylim([0 1])
xlim([0.5 3.5])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Thresh(Deg)')


subplot(3,2,5)
for i_orien = 1:length(B_ori)
    for i_off = 1:length(Off)
        plot([i_off i_off], output.FA{i_orien}.FA_confi(i_off,:),'Color',colors(i_orien,:))
        hold on
    end 
    scatter([1 2 3], output.FA{i_orien}.FA, 'MarkerEdgeColor',colors(i_orien,:))
end 
ylim([0 0.14])
xlim([0.5 3.5])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('FA rate')
% plot normalized by its own version
subplot(3,2,6)
for i_orien = 1:length(B_ori)
    
    scatter([1 2 3], output.FA{i_orien}.FA./output.FA{i_orien}.FA(1), 'MarkerEdgeColor',colors(i_orien,:))
    hold on
end 
ylim([0 6])
xlim([0.5 3.5])
xlabel('ISI (ms)')
ylabel('FA rate')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
supertitle([num2str(input.ID) ' - '  'randombase'])
%% figure
basecolor = {[0 0.6 0],[0 0 0],[1 0.4 0.2]};
subplot(1,2,1)
for i_orien = 1:length(B_ori)
    
    maxI = max(D_ori);
    minI = min(D_ori);
    xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
    h=line(xgrid, fitSall{i_orien,1}.modelFun(fitSall{i_orien,1}.coefEsts, xgrid), 'Color',basecolor{i_orien});
    hold on;
    plot(fitSall{i_orien,1}.intensityX,fitSall{i_orien,1}.fractCorrY, 'o','Color',basecolor{i_orien});
    %thresh = coefEsts(1)*[1 1];
    plot(fitSall{i_orien,1}.thresh*[1 1], [0 fitSall{i_orien,1}.threshY], '--','Color',basecolor{i_orien});
    plot(bootStatsall{i_orien,1}.ci95, fitSall{i_orien,1}.threshY*[1 1], 'Color',basecolor{i_orien});
    
    % set limits correctly
%     xLim = [min(xgrid) max(xgrid)].* [0.75 1.25];
%     xLim = 10.^ceil(log10(xLim) - [1 0]);
   
   
   
    text(7,1-0.05*i_orien,[num2str(B_ori(i_orien)) 'Deg'],'Color',basecolor{i_orien})
end

xlim([6 100])
ylim([0 1])

ylabel('Hit rate')
xlabel('\Delta Orientation')
set(gca,'XScale','log','TickDir','out')
axis square
subplot(1,2,2)
for i_orien = 1:length(B_ori)
    plot([i_orien i_orien], output.FA{i_orien}.all.FA_confi, 'Color',basecolor{i_orien})
    hold on
    scatter(i_orien,output.FA{i_orien}.all.FA,'MarkerEdgeColor',basecolor{i_orien})
end 
xlim([0.8 3.2])
ylim([0 0.1])

ylabel('FA rate')
xlabel('base orientation')
set(gca,'XTick',1:1:3,'XTicklabel',B_ori,'TickDir','out')
axis square



%% plot the threshold on first third and last third
figure

hold on
plot([1 1],first.Threshall_ci(2,:),'color',[0 0 0])
plot([2 2],first.Threshall_ci(1,:),'color',[0 0 0])
plot([3 3],first.Threshall_ci(3,:),'color',[0 0 0])
scatter([2 1 3], first.Threshall,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 0],'SizeData',50)
hold on
plot([1 1],last.Threshall_ci(2,:),'color',[0 0 0])
plot([2 2],last.Threshall_ci(1,:),'color',[0 0 0])
plot([3 3],last.Threshall_ci(3,:),'color',[0 0 0])
scatter([2 1 3], last.Threshall,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 1 1],'SizeData',50)
ylim([0 50])
xlim([0.5 3.5])

set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('Thresh delta orien (deg)')
text(0.7, 47, 'close: first 1/3')
text(0.7, 44, 'open: last 1/3')
title('i563')