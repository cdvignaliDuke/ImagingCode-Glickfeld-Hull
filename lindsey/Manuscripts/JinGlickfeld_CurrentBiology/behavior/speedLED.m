for led = 1:2
    if led==1
        input = inputcombined;
    end
    if led == 2
        input = inputcombined2;
    end
    Success = strcmp(input.trialOutcomeCell,'success');
    trialnum = length(Success);
    ReactTime = double(cell2mat(input.reactTimesMs));
    Early = strcmp(input.trialOutcomeCell,'failure');
    Miss = strcmp(input.trialOutcomeCell,'ignore');
    HoldTimeMs = double (cell2mat(input.holdTimesMs));
    reqHoldMs =  double (cell2mat(input.tTotalReqHoldTimeMs));
    raw_RT = ReactTime;
    tspeed = double(cell2mat(input.tDotSpeedDPS)) - double(cell2mat(input.tBaseDotSpeedDPS));
    if length(unique(tspeed))>7
        % manually adjust for each mouse
%         tspeed(tspeed>0.4 & tspeed <1) = 0.9;
%         tspeed(tspeed>1 & tspeed <2) = 1.8;
%         
%         tspeed(tspeed>2 & tspeed <4) = 3.5;
%         
%         tspeed(tspeed>4 & tspeed <8) = 7.2;
        
        
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
        Ouput{led}.rawRT{i_speed,1} = raw_RT(idx);
    end
    % get the reaction time distribution on window of 200-700 ms
    [Output{led}.target] = Contrast_HR(speed,SuccessN,MissN,tspeed,RT_HM);
    % get the FA rate
    [Output{led}.FA] = spd_FA(EarlyN,trialnum,HoldTimeMs,600); % 600 is the boundaries where the earliest change detect would occur
    % get the response for short vs long trials
%     % artificially define short: <=2s; long: >=3s
%     short_idx = reqHoldMs<=2000;
%     long_idx = reqHoldMs>=3000;
%     [Output{led}.Tshort] = Contrast_HR(speed,SuccessN(short_idx),MissN(short_idx),tspeed(short_idx),RT_HM(short_idx));
%     [Output{led}.Tlong] = Contrast_HR(speed,SuccessN(long_idx),MissN(long_idx),tspeed(long_idx),RT_HM(long_idx));
    % for binned trials
    [Output{led}.bintrial] = Speed_bin(input,speed,SuccessN,MissN,EarlyN,tspeed,trialnum,HoldTimeMs,RT_HM,600);
    
    [Output{led}.sdt] = Contrast_sdt(speed,Output{led});
    
    
    Output{led}.Infor.speed = speed;
    Output{led}.Infor.ID = input.ID;
    
    sep.tseq{led} = input.tseq;
    sep.EarlyN{led} = EarlyN;
    sep.SuccessN{led} = SuccessN;
    sep.HoldTimeMs{led} = HoldTimeMs;
    sep.tspeed{led} = tspeed;
    sep.RT_HM{led} = RT_HM;
    sep.MissN{led} = MissN;
  
  
    
end

[Output{1}.Outcome] = PreTrlSpeed(sep,600);
save(['i' num2str(input.ID) '_' 'speed-PMLED.mat'],'Output')