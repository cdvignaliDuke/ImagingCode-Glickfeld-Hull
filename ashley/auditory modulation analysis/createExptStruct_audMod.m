clear all
close all
ds = '_V16s';
doRedOnly = 0;
motionThreshold = 0.05;
%%
rc = behavConstsAV;
awData_temp
eval(['awData_audMod' ds])
ms = struct;
for iexp = 1:size(expt,2)
    subnum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    disp([mouse ' ' expDate])
    
    data = [];
    pupilSizeData = [];
    horizPosPupilData = [];
    vertPosPupilData = [];
    audStimDelayType = [];
    tOri = [];
    tTrStart = [];
    tTargetOn = [];
    wheelTrStartIdx = [];
    speedDataMS = [];
    speedTrStart = [];
    speedTarOn = [];
    doSpeed = 0;
    
    for irun = 1:expt(iexp).nrun
        runFolder = expt(iexp).runs(irun,:);
        runTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];
        
        % load cell timecourses, mworks file, and pupil timecourses
        fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,runFolder);
        if doRedOnly
            load(fullfile(fn,'red_timecourses.mat'))
        else
            load(fullfile(fn,'timecourses.mat'))
        end
        input = Load_SBXdataPlusMWorksData(subnum,expDate,runTime,mouse,runFolder,fName);
        
%         plotExptAudDelay(input,data_tc_subnp,mouse,expDate,irun)
        % combine data
        data = cat(1,data,data_tc_subnp);
        
        fillnan = nan(size(data_tc_subnp,1),1);
        doPupil = false;
%         if ~isempty(expt(iexp).eyeradrange)
%             load(fullfile(fn,'eyeTC.mat'));
%             pupilSizeData = cat(1,pupilSizeData,Area);
%             horizPosPupilData = cat(1,horizPosPupilData,Centroid(:,1));
%             vertPosPupilData = cat(1,vertPosPupilData,Centroid(:,2));
%             doPupil = true;
%         else
%             doPupil = false;
%         end
        
        % trial type
        trType = cell2mat(input.tBlock2TrialNumber);
        ntrials = length(trType);
        stimDelay = NaN(1,ntrials);
        if input.doPairedPulse
            if input.doAuditoryStim
                stimDelay(trType == 0) = 0;
            elseif input.block2DoAuditoryStim
                stimDelay(trType == 1) = 0;
            end
        else
            if input.doAuditoryStim
                if input.baseGratingContrast == 0
                    cyc = cell2mat(input.tCyclesOn(trType == 0));
                    off = cell2mat(input.tStimOffTimeMs(trType == 0));
                    on = cell2mat(input.tStimOnTimeMs(trType == 0));
                    stimDelay(trType == 0) = (off+on).*cyc;
                else
                    stimDelay(trType == 0) = 0;
                end
            end
            if input.block2DoAuditoryStim
                if input.block2BaseGratingContrast == 0
                    cyc = cell2mat(input.tCyclesOn(trType == 1));
                    off = input.block2StimOffTimeMs;
                    on = input.block2StimOnTimeMs;
                    stimDelay(trType == 1) = (off+on).*cyc;
                else                    
                    stimDelay(trType == 1) = 0;
                end
            end            
        end
        audStimDelayType = cat(2,audStimDelayType,stimDelay);

        % grating orientation
        tOri = cat(2,tOri,double(cell2mat(input.tGratingDirectionDeg)));
        
        % running info
        wheelValues = cellfun(@sum, input.wheelSpeedValues);
        if sum(wheelValues) > 0
            doSpeed = 1;
            [speedInTrialMS, speedTrialStartInd,speedTrialTargetInd] = getWheelSpeedFromInput(input);
            speedDataMS = cat(2,speedDataMS,speedInTrialMS);
            speedTrStart = cat(2,speedTrStart,speedTrialStartInd);
            speedTarOn = cat(2,speedTarOn,speedTrialTargetInd);
        end
        
        % variables with frame number offset
        nfr = size(data_tc_subnp,1);

        if irun == 1
            tTrStart = double(cell2mat(input.cFirstStim));
            tTargetOn = double(cell2mat(input.cTargetOn));
            offset = nfr;
            
        else
            tTrStart = cat(2,tTrStart, double(cell2mat(input.cFirstStim))+offset);
            tTargetOn = cat(2,tTargetOn, double(cell2mat(input.cTargetOn))+offset);
            offset = offset+nfr;
        end
        
    end
    
    % mouse info
    ms(iexp).subnum = subnum;
    ms(iexp).date = expDate;
    
    % labeled cells
    if expt(iexp).redLabel
        ms(iexp).redLabel = 1;
        load(fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,'redCellIdentity'))
        ms(iexp).redCellIdentity = redCellIdentity;
        ms(iexp).redCellType = expt(iexp).redCellType;
    else
        ms(iexp).redLabel = 0;
        ms(iexp).redCellIdentity = [];
        ms(iexp).redCellType = [];
    end
    
    % img params
    [nfr,nc] = size(data);
    ntrials = length(tOri);
    frRateHz = expt(iexp).frame_rate;

    % analysis params
    pre_event_frames = frRateHz;
    post_event_frames = frRateHz;
    post_long = frRateHz*10;
    wheelRes = input.speedIntervalMS;
    nWheelRec1S = 1000/wheelRes;
    pre_event_wheel = nWheelRec1S;
    post_event_wheel = nWheelRec1S;
    post_long_wheel = nWheelRec1S*6;
    
    
    % trial type index
    [delayType_ind, delayTypes] = findgroups(audStimDelayType);
    ntrtype = length(delayTypes)+2;
    
    % when does vis stim start?
    tVisStimStart = NaN(1,ntrials);
    tVisStimStart(delayType_ind > 1) = tTargetOn(delayType_ind > 1);
    tVisStimStart(delayType_ind == 1 | isnan(delayType_ind)) = tTrStart(delayType_ind == 1 | isnan(delayType_ind));
        
    % auditory only trials
    audOnly_ind = audStimDelayType >= 1000;
    tAudStimStart = tTrStart(audOnly_ind);
    
    % paired vis + aud trials
    visAudInd = audStimDelayType == 0;
    
    %% running speed 
    if doSpeed
        post_long_wheel = nWheelRec1S*(post_long/frRateHz);
        % aligned to trial start
        speedStartAlign = nan(pre_event_wheel+post_event_wheel);
        speedStartAlignLong = nan(pre_event_wheel+post_long_wheel);
        speedTargetAlign = nan(pre_event_wheel+post_event_wheel);
        for itrial = 1:ntrials
            ind = speedTrStart(itrial)-pre_event_wheel+1:speedTrStart(itrial)+post_event_wheel;
            speedStartAlign(:,itrial) = speedDataMS(ind);
            ind = speedTrStart(itrial)-pre_event_wheel+1:speedTrStart(itrial)+post_long_wheel;
            speedStartAlignLong(:,itrial) = speedDataMS(ind);
            ind = speedTarOn(itrial)-pre_event_wheel+1:speedTarOn(itrial)+post_event_wheel;
            speedTargetAlign(:,itrial) = speedDataMS(ind);
        end
        % aligned to targetcellfun(@(x,y) x(y-pre_event_wheel+1:y+post_event_wheel),speedDataMS,speedTarOn,'unif',0);
    end
    %% dF/F
    
    % aligned to target, F baseline immediately before target, 1s before
    % and after vis onset
    dff_targetAlignAndBL = getDFFEachTrial(data,tTargetOn,tTargetOn,...
    post_event_frames,frRateHz);
    
    % aligned to target, F baseline before trial start
    dff_targetAlignStartBL = getDFFEachTrial(data,tTrStart,tTargetOn,...
    post_event_frames,frRateHz);

    % aligned to start, short time-course
    dff_startAlign = getDFFEachTrial(data,tTrStart,tTrStart,...
        post_event_frames,frRateHz);
    
    % aligned to start, long time-course
    dff_startAlignLong = getDFFEachTrial(data,tTrStart,tTrStart,...
        post_long,frRateHz);

    %% motion trials
    withinTrialMotion = max(diff(squeeze(mean(dff_targetAlignAndBL,2)),1),[],1);
    noMotionInd = withinTrialMotion < motionThreshold;
    
    %% pupil
    if doPupil
        % aligned to target
        pupilSize_targetAlignStartBL = getPupilEachTrial(pupilSizeData,...
            tTrStart,tTargetOn,post_event_frames,frRateHz,'percent');
        pupilHorizPos_targetAlignStartBL = getPupilEachTrial(horizPosPupilData,...
            tTrStart,tTargetOn,post_event_frames,frRateHz,'subtraction');
        pupilVertPos_targetAlignStartBL = getPupilEachTrial(vertPosPupilData,...
            tTrStart,tTargetOn,post_event_frames,frRateHz,'subtraction');
        
        % aligned to start, short time-course
        pupilSize_startAlign = getPupilEachTrial(pupilSizeData,...
            tTrStart,tTrStart,post_event_frames,frRateHz,'percent');
        pupilHorizPos_startAlign = getPupilEachTrial(horizPosPupilData,...
            tTrStart,tTrStart,post_event_frames,frRateHz,'subtraction');
        pupilVertPos_startAlign = getPupilEachTrial(vertPosPupilData,...
            tTrStart,tTrStart,post_event_frames,frRateHz,'subtraction');
        
        % aligned to start, long time-course
        pupilSize_startAlignLong = getPupilEachTrial(pupilSizeData,...
            tTrStart,tTrStart,post_long,frRateHz,'percent');
        pupilHorizPos_startAlignLong = getPupilEachTrial(horizPosPupilData,...
            tTrStart,tTrStart,post_long,frRateHz,'subtraction');
        pupilVertPos_startAlignLong = getPupilEachTrial(vertPosPupilData,...
            tTrStart,tTrStart,post_long,frRateHz,'subtraction');

    end
    
    %% noise correlations
    delayNoiseCorr = cell(1,ntrtype+1);
    ncells = size(data,2);
    itiFrameInd = linspaceNDim(...
        tTrStart-pre_event_frames+1,tTrStart,pre_event_frames);
    itiRespAllTrials = squeeze(mean(reshape(data(itiFrameInd(:),:),...
        ntrials,pre_event_frames,ncells),2));
    itiNoiseCorr = corrcoef(itiRespAllTrials);
    delayNoiseCorr{ntrtype+1} = itiNoiseCorr;
    %% sort trial types
    delayTrials = cell(1,ntrtype);
    noiseCorrVisStim = cell(1,ntrtype);
    noiseCorrBeforeVisStim = cell(1,ntrtype);
    delayLong = cell(1,ntrtype-2);
    delayPupilSize = cell(1,ntrtype);
    delayHorizPosPupil = cell(1,ntrtype);
    delayVertPosPupil = cell(1,ntrtype);
    trType_str = cell(1,ntrtype);
    trOri = cell(1,ntrtype);
    trSpeed = cell(1,ntrtype);
    trSpeedLong = cell(1,ntrtype-2);
    delayLongPsz = cell(1,ntrtype-2);
    delayLongHoriz = cell(1,ntrtype-2);
    delayLongVert = cell(1,ntrtype-2);    
    for idelay = 1:ntrtype
        if idelay == ntrtype % auditory only condition
            ind = audOnly_ind & noMotionInd;
            delayTrials{idelay} = dff_startAlign(:,:,ind);
            noiseCorrVisStim{idelay} = getNoiseCorrSortedData(...
                dff_startAlign(pre_event_frames+1:end,:,ind));
            noiseCorrBeforeVisStim{idelay} = getNoiseCorrSortedData(...
                dff_startAlign(1:pre_event_frames,:,ind));
            trType_str{idelay} = 'aud only';
            trOri{idelay} = tOri(ind);
            if doPupil
                delayPupilSize{idelay} = pupilSize_startAlign(:,ind);
                delayHorizPosPupil{idelay} = pupilHorizPos_startAlign(:,ind);
                delayVertPosPupil{idelay} = pupilVertPos_startAlign(:,ind);                
            end
            if doSpeed
                trSpeed{idelay} = speedStartAlign(:,ind);
            end
        elseif idelay == ntrtype - 1 % visual only condition
            if sum(isnan(delayType_ind)) == 0
                ind = audStimDelayType >= 6000 & noMotionInd;
                delayTrials{idelay} = dff_targetAlignAndBL(:,:,ind);
                noiseCorrVisStim{idelay} = getNoiseCorrSortedData(...
                    dff_targetAlignAndBL(pre_event_frames+1:end,:,ind));
                noiseCorrBeforeVisStim{idelay} = getNoiseCorrSortedData(...
                    dff_targetAlignAndBL(1:pre_event_frames,:,ind));
                if doPupil
                    delayPupilSize{idelay} = pupilSize_targetAlignStartBL(:,ind);
                    delayHorizPosPupil{idelay} = pupilHorizPos_targetAlignStartBL(:,ind);
                    delayVertPosPupil{idelay} = pupilVertPos_targetAlignStartBL(:,ind);                
                end
                if doSpeed
                    trSpeed{idelay} = speedTargetAlign(:,ind);
                end
            else                
                ind = isnan(delayType_ind) & noMotionInd;
                delayTrials{idelay} = dff_startAlign(:,:,ind);
                noiseCorrVisStim{idelay} = getNoiseCorrSortedData(...
                    dff_startAlign(pre_event_frames+1:end,:,ind));
                noiseCorrBeforeVisStim{idelay} = getNoiseCorrSortedData(...
                    dff_startAlign(1:pre_event_frames,:,ind));
                if doPupil
                    delayPupilSize{idelay} = pupilSize_startAlign(:,ind);
                    delayHorizPosPupil{idelay} = pupilHorizPos_startAlign(:,ind);
                    delayVertPosPupil{idelay} = pupilVertPos_startAlign(:,ind);                 
                end
                if doSpeed
                    trSpeed{idelay} = speedStartAlign(:,ind);
                end
            end
            trType_str{idelay} = 'vis only';
            trOri{idelay} = tOri(ind);
        elseif idelay == 1
            ind = delayType_ind == idelay & noMotionInd;
            delayTrials{idelay} = dff_startAlign(:,:,ind);
            noiseCorrVisStim{idelay} = getNoiseCorrSortedData(...
                dff_startAlign(pre_event_frames+1:end,:,ind));
            noiseCorrBeforeVisStim{idelay} = getNoiseCorrSortedData(...
                dff_startAlign(1:pre_event_frames,:,ind));
            delayLong{idelay} = dff_startAlignLong(:,:,ind);
            trType_str{idelay} = num2str(delayTypes(idelay));
            trOri{idelay} = tOri(ind);            
            if doPupil
                delayPupilSize{idelay} = pupilSize_startAlign(:,ind);
                delayHorizPosPupil{idelay} = pupilHorizPos_startAlign(:,ind);
                delayVertPosPupil{idelay} = pupilVertPos_startAlign(:,ind);   
                delayLongPsz{idelay} = pupilSize_startAlignLong(:,ind);
                delayLongHoriz{idelay} = pupilHorizPos_startAlignLong(:,ind);
                delayLongVert{idelay} = pupilVertPos_startAlignLong(:,ind);                
            end
            if doSpeed
                trSpeed{idelay} = speedTargetAlign(:,ind);
                trSpeedLong{idelay} = speedStartAlignLong(:,ind);
            end
        else
            ind = delayType_ind == idelay & noMotionInd;
            delayTrials{idelay} = dff_targetAlignAndBL(:,:,ind);
            noiseCorrVisStim{idelay} = getNoiseCorrSortedData(...
                dff_targetAlignAndBL(pre_event_frames+1:end,:,ind));
            noiseCorrBeforeVisStim{idelay} = getNoiseCorrSortedData(...
                dff_targetAlignAndBL(1:pre_event_frames,:,ind));
            delayLong{idelay} = dff_startAlignLong(:,:,ind);
            trType_str{idelay} = num2str(delayTypes(idelay));
            trOri{idelay} = tOri(ind);            
            if doPupil
                delayPupilSize{idelay} = pupilSize_targetAlignStartBL(:,ind);
                delayHorizPosPupil{idelay} = pupilHorizPos_targetAlignStartBL(:,ind);
                delayVertPosPupil{idelay} = pupilVertPos_targetAlignStartBL(:,ind);   
                delayLongPsz{idelay} = pupilSize_startAlignLong(:,ind);
                delayLongHoriz{idelay} = pupilHorizPos_startAlignLong(:,ind);
                delayLongVert{idelay} = pupilVertPos_startAlignLong(:,ind);         
            end
            if doSpeed
                trSpeed{idelay} = speedTargetAlign(:,ind);
                trSpeedLong{idelay} = speedStartAlignLong(:,ind);
            end
        end
    end
    
    ms(iexp).delayTrials = delayTrials;
    ms(iexp).noiseCorrVisStim = noiseCorrVisStim;
    ms(iexp).noiseCorrBeforeVisStim = noiseCorrBeforeVisStim;
    ms(iexp).delayLong = delayLong;
    ms(iexp).pupil = {delayPupilSize,delayHorizPosPupil,delayVertPosPupil};
    ms(iexp).pupilLong = {delayLongPsz,delayLongHoriz,delayLongVert};
    ms(iexp).audDelayType = trType_str;
    ms(iexp).trOri = trOri;
    ms(iexp).orientations = unique(tOri);
    ms(iexp).frRateHz = frRateHz;
    ms(iexp).pre_event_frames = pre_event_frames;
    if doSpeed
        ms(iexp).trSpeed = trSpeed;
        ms(iexp).trSpeedLong = trSpeedLong;
        ms(iexp).pre_event_wheel = pre_event_wheel;
        ms(iexp).wheelRateHz = nWheelRec1S;
    end
end

%% save 
fnout = fullfile(rc.ashleyAnalysis,'Expt Summaries',['awData_audMod' ds]);
save(fullfile(fnout,['awData_audMod' ds '_CaSummary']),'ms');
