function [first] = randombase_learning(B_ori,StimOff,B_Orien,D_ori,Off,SuccessN,MissN,tCycleNum,D_Orien,RT_HM,EarlyN,Leverup,Leverdown,TT_Orien,tStimOffTimes)

% get the first third

for i_orien = 1:length(B_ori)
    
    [first.target{i_orien}] = ISI_HR(StimOff(B_Orien==B_ori(i_orien)),Off,D_ori,SuccessN(B_Orien==B_ori(i_orien)),MissN(B_Orien==B_ori(i_orien)),tCycleNum(B_Orien==B_ori(i_orien)),D_Orien(B_Orien==B_ori(i_orien)),RT_HM(B_Orien==B_ori(i_orien)));

    Input = {};
    Input.stimOnTimeMs =  100;
    Input.tStimOffTimes = tStimOffTimes(B_Orien==B_ori(i_orien));
   
    [first.FA{i_orien}] =  ISI_FA_N(EarlyN(B_Orien==B_ori(i_orien)),Input,Off,tCycleNum(B_Orien==B_ori(i_orien)),sum(B_Orien==B_ori(i_orien)),Leverup(B_Orien==B_ori(i_orien)),Leverdown(B_Orien==B_ori(i_orien)));
  
    
end

% break up by absolute target value to calculate hit rate
for i_orien = 1:length(B_ori)
    
    T_orien(i_orien,:) = unique(TT_Orien(B_Orien==B_ori(i_orien)));
    
    [first.Abstarget{i_orien}] = ISI_HR(StimOff(B_Orien==B_ori(i_orien)),Off,T_orien(i_orien,:),SuccessN(B_Orien==B_ori(i_orien)),MissN(B_Orien==B_ori(i_orien)),tCycleNum(B_Orien==B_ori(i_orien)),TT_Orien(B_Orien==B_ori(i_orien)),RT_HM(B_Orien==B_ori(i_orien)));
   

end 

% get the threshold
NBootstrapReps = 1000;
for i = 1:length(B_ori)
    Hits = first.target{i}.all.c_hit;
    trialAll =  first.target{i}.all.HT_num + first.target{i}.all.Miss_num;
    first.fitSall{i,1} = weibullFitLG(D_ori, Hits',1, 1, {'nTrials', trialAll'});
    first.Threshall(i) = first.fitSall{i,1}.thresh;
    
    [first.bootStatsall{i,1}]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,D_ori,1, 1);
    first.Threshall_ci(i,1:2) = first.bootStatsall{i,1}.ci95;
    % collapsed all offs for absolute target
    Hits = first.Abstarget{i}.all.c_hit;
    trialAll =  first.Abstarget{i}.all.HT_num + first.Abstarget{i}.all.Miss_num;
    first.AbsfitSall{i,1} = weibullFitLG(T_orien(i,:), Hits',1, 1, {'nTrials', trialAll'});
    first.AbsThreshall(i) = first.AbsfitSall{i,1}.thresh;
    
    [first.AbsbootStatsall{i,1}]=BootstrapWeibullFit(trialAll',  Hits',NBootstrapReps,T_orien(i,:),1, 1);
    first.AbsThreshall_ci(i,1:2) = first.AbsbootStatsall{i,1}.ci95;
    
end 

end
