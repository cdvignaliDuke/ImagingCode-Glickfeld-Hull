%% contrast detection

for led = 1:2
    if led==1
        input = inputcombined;
    end
    if led == 2
        input = inputcombined2;
    end
    
    if size(input.tCyclesOn,2) >0 % if use FS code
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
        tContrast =  double(cell2mat(input.tGratingContrast));
        Contrast = unique(tContrast);
        
        
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
            if length(unique(cell2mat(inputcombined.tGratingContrast)))>7
                % manually adjust for each mouse
             
            end
            
            
            [Output{led}.FA_percent]  = FA_Contrast(EarlyN,HoldTimeMs,trialnum,900); %700+200 is the boundary between fideget and early release for the FS code
            %[Output{led}.FA_percent]  = FA_Trl(EarlyN,tCycleNum,trialnum,Leverup,Leverdown);
            %[Output{led}.bintrial] = Contrast_bin(input,Contrast,SuccessN,MissN,EarlyN,tContrast,trialnum,HoldTimeMs,RT_HM,900);
            [Output{led}.target] = Contrast_HR(Contrast,SuccessN,MissN,tContrast,RT_HM);
            % needs to modify trial dependence code
            [Output{led}.bintrial] = ISI_Samp(input,250,Contrast,SuccessN,MissN,EarlyN,tCycleNum,tContrast,trialnum,Leverup,Leverdown,1);
            [Output{led}.FA] = Fixed_FA(EarlyN,input,tCycleNum,trialnum,Leverup,Leverdown);

 
            [Output{led}.sdt] = Contrast_sdt(Contrast,Output{led});
            
            
            Output{led}.Infor.Contrast = Contrast;           
            Output{led}.Infor.ID = input.ID;
            
            
            sep.tseq{led} = tseq;
            sep.EarlyN{led} = EarlyN;
            sep.SuccessN{led} = SuccessN;
            sep.HoldTimeMs{led} = HoldTimeMs;
            sep.tContrast{led} = tContrast;
            sep.RT_HM{led} = RT_HM;
            sep.MissN{led} = MissN;
            sep.tCycle{led} = tCycleNum;
            sep.Leverup{led} = Leverup;
            sep.Leverdown{led} = Leverdown;
            
        end
        
    else  % if use HDC code
        
        Success = strcmp(input.trialOutcomeCell,'success');
        trialnum = length(Success);
        ReactTime = double(cell2mat(input.reactTimesMs));
        Early = strcmp(input.trialOutcomeCell,'failure');
        Miss = strcmp(input.trialOutcomeCell,'ignore');
        HoldTimeMs = double (cell2mat(input.holdTimesMs));
        tOrientation = double(cell2mat(input.tGratingDirectionDeg));
        
        Orien = unique(tOrientation);
        tseq = input.tseq;
        tContrast =  double(cell2mat(input.tGratingContrast));
        Contrast = unique(tContrast);
        
        
        
        INPUT{led} = input;
        idx_HM = find(Early==0); % index of Hit and Miss trials
        % just use ReactTime in HDC condtions as new RT
        
        
        RT_HM = zeros(1,length(ReactTime));
        RT_HM(idx_HM) = ReactTime(idx_HM); 
        
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
        if length(unique(cell2mat(inputcombined.tGratingContrast)))>7
            % manually adjust for each mouse
            
            
            
            Contrast = unique( tContrast);
        end
        
        
        
        [Output{led}.FA_percent]  = FA_Contrast(EarlyN,HoldTimeMs,trialnum,600); % 400+200 is the boundary between fideget and early release
       
        [Output{led}.bintrial] = Contrast_bin(input,Contrast,SuccessN,MissN,EarlyN,tContrast,trialnum,HoldTimeMs,RT_HM,600);
        
        [Output{led}.target] = Contrast_HR(Contrast,SuccessN,MissN,tContrast,RT_HM);
        
        [Output{led}.FA] = HDC_FA(EarlyN,trialnum,HoldTimeMs,600); % 400+200  is the boundary between fideget and early release
        
        [Output{led}.sdt] = Contrast_sdt(Contrast,Output{led});
      
        Output{led}.Infor.Contrast = Contrast;
             
        Output{led}.Infor.ID = input.ID;
%         
        sep.tseq{led} = tseq;
        sep.EarlyN{led} = EarlyN;
        sep.SuccessN{led} = SuccessN;
        sep.HoldTimeMs{led} = HoldTimeMs;
        sep.tContrast{led} = tContrast;
        sep.RT_HM{led} = RT_HM;
        sep.MissN{led} = MissN;
%         
        
    end
end





if  size(input.tCyclesOn,2) >0
    % calculate the trial outcome history dependence
    [Output{1}.Outcome] = PreTrlContrast(sep,900);
    
    save(['i' num2str(input.ID) '_' 'contrastFS-PMLED.mat'],'Output')
else
    [Output{1}.Outcome] = PreTrlContrast(sep,600);
    
    save(['i' num2str(input.ID) '_' 'contrastHDC-LED_new.mat'],'Output')
end
