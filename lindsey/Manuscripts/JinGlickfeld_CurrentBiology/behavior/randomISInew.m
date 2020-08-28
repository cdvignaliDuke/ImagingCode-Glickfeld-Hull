%% this is corrected for a new RT window analysis
% Hit/FA RT window is ajusted as 200-550 ms after visual onset

if exist('inputcombined2','var')
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
        Leverdown = double(cell2mat(input.tLeverPressTimeMs));
        Leverup = double(cell2mat(input.tLeverReleaseTimeMs)); % time when trial start and trial stop
        tOrientation = double(cell2mat(input.tGratingDirectionDeg));
        CycleNum = double(cell2mat(input.nCyclesOn)); % cycle number allocate to each trial, regardless of how animal react
        tCycleNum = double(cell2mat(input.tCyclesOn)); % cycle number animal actually complete.
        Orien = unique(tOrientation);
        tseq = input.tseq;
        
        % if block2 is on, then tStimOffTimes has extra 250ms in the
        % beginning, get rid of this
        if ~isfield(input,'doLaserstim') && size(input.tStimOffTimes,2)~=0
            input.tStimOffTimes = cellfun(@(x) x(2:end),input.tStimOffTimes, 'UniformOutput', false);
            
            INPUT{led} = input;
            % get the target interval time
            StimOff = NaN(1,trialnum); % interval before target stimulus
            StimOff_pre = NaN(1,trialnum); % interval before last cycle
            idx_HM = find(Early==0); % index of Hit and Miss trials
            Cycle = {};
            Cycle_pre = {};
            for i=1:length(idx_HM)
                idx_temp = [];
                cycle = [];
                trialoffs= [];
                idx_temp = idx_HM(i);% trial of hit or miss trial
                cycle = tCycleNum(idx_temp);
                Triallength = double(input.targetStimOnMs{idx_temp})-Leverdown(idx_temp);
                trialoffs = double(input.tStimOffTimes{idx_temp});
                StimOff(idx_temp) = trialoffs(cycle);% interval before target
                StimOff_pre(idx_temp) = trialoffs(cycle-1); % interval before last cycle
                
                
                a(i)=abs(Triallength-sum(trialoffs(2:cycle+1))- (double(input.stimOnTimeMs).*cycle));
                
            end
            
            % calculate the unique cycle number for each condition, ideally should be
            Off = unique(StimOff(idx_HM)); % unique off times, [250 500 750]
            
            for i_off = 1:length(Off)
                
                Cycle{i_off,1} = unique(CycleNum(StimOff==Off(i_off))); % cycle number for each off before target
                Cycle_pre{i_off,1} = unique(CycleNum(StimOff_pre==Off(i_off))); % cycle number for each off before last cycle
            end
            
           
            
            % correct the new_RT on hit and misese, not work when
            % 110217-120717 when dotargetaftermiss is on
            new_RT = cellfun(@minus,input.tLeverReleaseTimeMs,input.targetStimOnMs, 'UniformOutput', false);
            
            
            RT_HM = zeros(1,length(ReactTime));
            RT_HM(~cellfun('isclass',new_RT, 'double')) = cell2mat(new_RT(~cellfun('isclass',new_RT, 'double')));
            
            RT_HM (cellfun('isclass',new_RT, 'double')) = cell2mat(new_RT(cellfun('isclass',new_RT, 'double')));
%             
            % adjust the Success, Miss trial conditions
            % throw too fast response into early trials
            SuccessN = zeros(1,length(Success));
            MissN = zeros(1,length(Miss));
            EarlyN = Early;
           
            
            for i=1:length(idx_HM)
                idx_temp = [];
                idx_temp = idx_HM(i);% trial of hit or miss trial
%                 RT_temp = [];
%                 cycle_temp = tCycleNum(idx_temp);
%                 trial_offs = double(input.tStimOffTimes{idx_temp});
%                 RT_temp = Leverup(idx_temp) - Leverdown(idx_temp) - 100.*cycle_temp -sum(trial_offs(1:cycle_temp));
%                 RT_HM(idx_temp) = RT_temp;
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
            % calculate the hit/FA rate based on trial length
             [Output{led}.Trl] = ISI_Trl(StimOff,input,Off,Orien,SuccessN,MissN,EarlyN,tCycleNum,tOrientation,trialnum,Leverup,Leverdown);
            % calculate the hit/FA rate matched orientations for each bin?
            [Output{led}.bintrial] = ISI_Samp(input,Off,Orien,SuccessN,MissN,EarlyN,tCycleNum,tOrientation,trialnum,Leverup,Leverdown);
            
            
            % calculate the FA rate based on percentage of actual trials
            
            [Output{led}.FA_percent]  = FA_Trl(EarlyN,tCycleNum,trialnum,Leverup,Leverdown);
            
            % calculate Hit Rate and RT on Hits over different orientation changes separate by offs before the target, need to write into function
            [Output{led}.target] = ISI_HR(StimOff,Off,Orien,SuccessN,MissN,tCycleNum,tOrientation,RT_HM);
            [Output{led}.pre] = ISI_HR(StimOff_pre,Off,Orien,SuccessN,MissN,tCycleNum,tOrientation,RT_HM);
            
            
            % calculate FA rate, should start from cycle 3..
            [Output{led}.FA] = ISI_FA_N(EarlyN,input,Off,tCycleNum,trialnum,Leverup,Leverdown);
            % caculate the release rate during the BS visual presentations
            [Output{led}.Re] = ISI_Re(EarlyN,input,Off,tCycleNum,trialnum,Leverup,Leverdown);
            % caculate the RT distribution for multiple conditions
            [Output{led}.RT] = ISI_RT(EarlyN,input,Off,tCycleNum,Leverup,Leverdown,SuccessN,Orien,tOrientation);
            % calculate dprime and criterion
            [Output{led}.sdt] = ISI_sdt(Off,Orien,Output{led});
            % seperate the
            [Output{led}.con] = ISI_con_N(EarlyN,input,Off,tCycleNum,trialnum,Leverup,Leverdown,SuccessN,tOrientation,Orien,RT_HM,StimOff,StimOff_pre,MissN);
            
            Output{led}.Infor.Orien = Orien;
            Output{led}.Infor.Off = Off;
            Output{led}.Infor.ID = input.ID;
            % get the adjusted reaction time window for both the hit/FA rate.
            sep.tseq{led} = tseq;
            sep.Stimoff{led} = StimOff;
            sep.EarlyN{led} = EarlyN;
            sep.tCycleNum{led} = tCycleNum;
            sep.Leverup{led} = Leverup;
            sep.Leverdown{led} = Leverdown;
            sep.SuccessN{led} = SuccessN;
            sep.tOrientation{led} = tOrientation;
            sep.RT_HM{led} = RT_HM;
            sep.MissN{led} = MissN;
        end
        
        if size(input.tStimOffTimes,2)==0 % seperate consider the fixed 250 ms condition
            INPUT{led} = input;
            idx_HM = find(Early==0); % index of Hit and Miss trials
            % correct the new_RT on hit and misese
            new_RT = cellfun(@minus,input.tLeverReleaseTimeMs,input.targetStimOnMs, 'UniformOutput', false);
            
            
            RT_HM = zeros(1,length(ReactTime));
            RT_HM(~cellfun('isclass',new_RT, 'double')) = cell2mat(new_RT(~cellfun('isclass',new_RT, 'double')));
            
            RT_HM (cellfun('isclass',new_RT, 'double')) = cell2mat(new_RT(cellfun('isclass',new_RT, 'double')));
            
            % adjust the Success, Miss trial conditions
            % throw too fast response into early trials
            SuccessN = zeros(1,length(Success));
            MissN = zeros(1,length(Miss));
            EarlyN = Early;
            
            for i=1:length(idx_HM)
                idx_temp = [];
                idx_temp = idx_HM(i);% trial of hit or miss trial
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
            
            StimOff = repmat(250,1,trialnum);
            %binned the orientation if there are many orientation levels
            if length(unique(cell2mat(inputcombined.tGratingDirectionDeg)))>7
                % manually adjust for each mouse
                
                
                %for i511
%                 tOrientation( tOrientation>=10 &  tOrientation<=20) = 16.6;
%                 tOrientation( tOrientation>=22 &  tOrientation<=28) = 23;
%                 tOrientation( tOrientation>=29 &  tOrientation<=35) = 31.8;
%                 tOrientation( tOrientation>=36 &  tOrientation<=46) = 43;
%                 tOrientation( tOrientation>=48 &  tOrientation<=80) =64;




                Orien = unique( tOrientation);
            end
            
            
            
            [Output{led}.FA_percent]  = FA_Trl(EarlyN,tCycleNum,trialnum,Leverup,Leverdown);
            [Output{led}.bintrial] = ISI_Samp(input,250,Orien,SuccessN,MissN,EarlyN,tCycleNum,tOrientation,trialnum,Leverup,Leverdown);
            
            [Output{led}.target] = ISI_HR(StimOff,250,Orien,SuccessN,MissN,tCycleNum,tOrientation,RT_HM);
            [Output{led}.FA] = Fixed_FA(EarlyN,input,tCycleNum,trialnum,Leverup,Leverdown);
            [Output{led}.sdt] = Fixed_sdt(Orien,Output{led});
            % deal with trial length effects later..
           
            Output{led}.Infor.Orien = Orien;
            Output{led}.Infor.Off = 250;
            Output{led}.Infor.ID = input.ID;
            
            sep.tseq{led} = tseq;
            sep.Stimoff{led} = StimOff;
            sep.EarlyN{led} = EarlyN;
            sep.tCycleNum{led} = tCycleNum;
            sep.Leverup{led} = Leverup;
            sep.Leverdown{led} = Leverdown;
            sep.SuccessN{led} = SuccessN;
            sep.tOrientation{led} = tOrientation;
            sep.RT_HM{led} = RT_HM;
            sep.MissN{led} = MissN;
        end
        
        
        
        
    end
    
    if  size(input.tStimOffTimes,2)==0
        % calculate the trial outcome history dependence
        [Output{1}.Outcome] = PreTrl(sep,INPUT,250,0);
        
        save(['i' num2str(input.ID) '_' 'Fixed-LED.mat'],'Output')
    else
        [Output{1}.Outcome] = PreTrl(sep,INPUT,Off,1);
        
        save(['i' num2str(input.ID) '_' 'ISI-LED.mat'],'Output')
    end
    
end

%% deal with single condition
if ~exist('inputcombined2','var')
    input = inputcombined;
    Success = strcmp(input.trialOutcomeCell,'success');
    trialnum = length(Success);
    ReactTime = double(cell2mat(input.reactTimesMs));
    Early = strcmp(input.trialOutcomeCell,'failure');
    Miss = strcmp(input.trialOutcomeCell,'ignore');
    HoldTimeMs = double (cell2mat(input.holdTimesMs));
    Leverdown = double(cell2mat(input.tLeverPressTimeMs));
    Leverup = double(cell2mat(input.tLeverReleaseTimeMs)); % time when trial start and trial stop
    tOrientation = double(cell2mat(input.tGratingDirectionDeg));
    tbase = double(cell2mat(input.tBaseGratingDirectionDeg));
    if unique(tbase) ==0
        
    else
    tOrientation = tOrientation - unique(tbase);
    end
    CycleNum = double(cell2mat(input.nCyclesOn)); % cycle number allocate to each trial, regardless of how animal react
    tCycleNum = double(cell2mat(input.tCyclesOn)); % cycle number animal actually complete.
    Orien = unique(tOrientation);
    
    % get the target interval time
    StimOff = NaN(1,trialnum); % interval before target stimulus
    StimOff_mean = NaN(1,trialnum);
    StimOff_pre = NaN(1,trialnum); % interval before last cycle
    StimOff_premean = NaN(1,trialnum);
    beforebins = NaN(1,trialnum);
    idx_HM = find(Early==0); % index of Hit and Miss trials
    
    for i=1:length(idx_HM)
        idx_temp = [];
        cycle = [];
        trialoffs= [];
        idx_temp = idx_HM(i);% trial of hit or miss trial
        cycle = tCycleNum(idx_temp);
        Triallength = double(input.targetStimOnMs{idx_temp})-Leverdown(idx_temp);
        trialoffs = double(input.tStimOffTimes{idx_temp});
        StimOff(idx_temp) = trialoffs(cycle);% interval before target
        StimOff_pre(idx_temp) = trialoffs(cycle-1); % interval before last cycle
        % bin the interval acording to the sequence of offs
        % 3 bins,off center 333,500,666, but still call 250,500,750
        %
        %         if mean(trialoffs(1:(cycle-1)))<=400
        %         temp=250;
        %         end
        %         if mean(trialoffs(1:(cycle-1)))>=600
        %         temp=750;
        %         end
        %         if mean(trialoffs(1:(cycle-1)))>400 && mean(trialoffs(1:(cycle-1)))<600
        %         temp=500;
        %         end
        %         beforebins(idx_temp)=mean(trialoffs(1:(cycle-1)));
        %
        %         StimOff_premean(idx_temp) = temp;
        %
        %         % same for the target response
        %          if mean(trialoffs(1:(cycle)))<=400
        %         temp=250;
        %         end
        %         if mean(trialoffs(1:(cycle)))>=600
        %         temp=750;
        %         end
        %         if mean(trialoffs(1:(cycle)))>400 && mean(trialoffs(1:cycle))<600
        %         temp=500;
        %         end
        %
        %         StimOff_mean(idx_temp) = temp;
        
        % or sort out trials that all has 250ms off/ 750ms off/ 500ms off
        temp=0;
        
        if length(unique(trialoffs(1:(cycle-1)))) ==1 % for all the fixed values
            if unique(trialoffs(1:(cycle-1)))==250
                temp=250;
            end
            if unique(trialoffs(1:(cycle-1)))==750
                temp=750;
            end
            if unique(trialoffs(1:(cycle-1)))==500
                temp=500;
            end
        end
        beforebins(idx_temp)=mean(trialoffs(1:(cycle-1)));
        
        StimOff_premean(idx_temp) = temp;
        
        % same for the target response
        temp=0;
        if length(unique(trialoffs(1:(cycle)))) ==1
            if unique(trialoffs(1:(cycle)))==250
                temp=250;
            end
            if unique(trialoffs(1:(cycle)))==750
                temp=750;
            end
            if unique(trialoffs(1:(cycle)))==500
                temp=500;
            end
        end
        StimOff_mean(idx_temp) = temp;
        
        
        a(i)=abs(Triallength-sum(trialoffs(1:cycle))- (double(input.stimOnTimeMs).*cycle));
        
    end
    
    % calculate the unique cycle number for each condition, ideally should be
    Off = unique(StimOff(idx_HM)); % unique off times, [250 500 750]
    
    for i_off = 1:length(Off)
        
        Cycle{i_off,1} = unique(CycleNum(StimOff==Off(i_off))); % cycle number for each off before target
        Cycle_pre{i_off,1} = unique(CycleNum(StimOff_pre==Off(i_off))); % cycle number for each off before last cycle
    end
    
    
    % correct the new_RT on hit and misese
    new_RT = cellfun(@minus,input.tLeverReleaseTimeMs,input.targetStimOnMs, 'UniformOutput', false);
    
    
    RT_HM = zeros(1,length(ReactTime));
    RT_HM(~cellfun('isclass',new_RT, 'double')) = cell2mat(new_RT(~cellfun('isclass',new_RT, 'double')));
    
    RT_HM (cellfun('isclass',new_RT, 'double')) = cell2mat(new_RT(cellfun('isclass',new_RT, 'double')));
    
    % adjust the Success, Miss trial conditions
    % throw too fast response into early trials
    SuccessN = zeros(1,length(Success));
    MissN = zeros(1,length(Miss));
    EarlyN = Early;
    
    for i=1:length(idx_HM)
        idx_temp = [];
        idx_temp = idx_HM(i);% trial of hit or miss trial
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
    
    % calculate the hit/FA rate based on trial length
    [output.Trl] = ISI_Trl(StimOff,input,Off,Orien,SuccessN,MissN,EarlyN,tCycleNum,tOrientation,trialnum,Leverup,Leverdown);
    
    % calculate Hit Rate and RT on Hits over different orientation changes separate by offs before the target, need to write into function
    [output.target] = ISI_HR(StimOff,Off,Orien,SuccessN,MissN,tCycleNum,tOrientation,RT_HM); % for the new time window hit rate
    
    [output.targetbin] = ISI_HR(StimOff_mean,Off,Orien,SuccessN,MissN,tCycleNum,tOrientation,RT_HM); % history dependence
    
    [output.pre] = ISI_HR(StimOff_pre,Off,Orien,SuccessN,MissN,tCycleNum,tOrientation,RT_HM);
    [output.prebin] = ISI_HR(StimOff_premean,Off,Orien,SuccessN,MissN,tCycleNum,tOrientation,RT_HM);
    
    
    % calculate FA rate, should start from cycle 3..and adjust the new
    % window
    [output.FA] =  ISI_FA_N(EarlyN,input,Off,tCycleNum,trialnum,Leverup,Leverdown);
    
    [output.FAbin] =  ISI_FA_bin(EarlyN,input,Off,tCycleNum,trialnum,Leverup,Leverdown);
    
    % caculate the release rate during the BS visual presentations
    [output.Re] = ISI_Re(EarlyN,input,Off,tCycleNum,trialnum,Leverup,Leverdown);
    % caculate the RT distribution for multiple conditions
    [output.RT] = ISI_RT(EarlyN,input,Off,tCycleNum,Leverup,Leverdown,SuccessN,Orien,tOrientation);
    % calculate dprime and criterion
    [output.sdt] = ISI_sdt(Off,Orien,output);
    % seperate the
    [output.con] = ISI_con_N(EarlyN,input,Off,tCycleNum,trialnum,Leverup,Leverdown,SuccessN,tOrientation,Orien,RT_HM,StimOff,StimOff_pre,MissN);
    
    [output.conbin] = ISI_conbin_N(EarlyN,input,Off,tCycleNum,trialnum,Leverup,Leverdown,SuccessN,tOrientation,Orien,RT_HM,StimOff,StimOff_premean,MissN);
    
    
    output.Infor.Orien = Orien;
    output.Infor.Off = Off;
    output.Infor.ID = input.ID;
    save(['i' num2str(input.ID) '_' 'ISI-bin.mat'],'output')
    %save(['i' num2str(input.ID) '_' 'ISI-bad.mat'],'output')
end