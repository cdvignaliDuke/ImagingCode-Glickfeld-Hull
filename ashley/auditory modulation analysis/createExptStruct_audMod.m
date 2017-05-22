clear all
close all
ds = '_V1';
%%
rc = behavConstsAV;
eval(['awData_audMod' ds])
ms = struct;
for iexp = 1:size(expt,2)
    subNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    
    data = [];
    audStimDelayType = [];
    tOri = [];
    tTrStart = [];
    tTargetOn = [];
    
    for irun = 1:expt(iexp).nrun
        runFolder = expt(iexp).runs(irun,:);
        runTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];
        
        % load cell timecourses and mworks file
        fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,runFolder);
        load(fullfile(fn,'timecourses.mat'))
        input = Load_SBXdataPlusMWorksData(subNum,expDate,runTime,mouse,runFolder,fName);
        
        % combine data
        data = cat(1,data,data_tc_subnp);
        
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
            stimDelay(trType == 0) = input.stimOnTimeMs+input.stimOffTimeMs;
            stimDelay(trType == 1) = input.block2StimOnTimeMs+input.block2StimOffTimeMs;
        end
        audStimDelayType = cat(2,audStimDelayType,stimDelay);

        % grating orientation
        tOri = cat(2,tOri,double(cell2mat(input.tGratingDirectionDeg)));

        % variables with frame number offset
        nfr = size(data_tc_subnp,1);
        if irun == 1
            tTrStart = double(cell2mat(input.cLeverDown));
            tTargetOn = double(cell2mat(input.cTargetOn));
            offset = nfr;
        else
            tTrStart = cat(2,tTrStart, double(cell2mat(input.cLeverDown))+offset);
            tTargetOn = cat(2,tTargetOn, double(cell2mat(input.cTargetOn))+offset);
            offset = offset+nfr;
        end
        
    end
    % img params
    nfr = size(data,1);
    nc = size(data,2);
    ntrials = length(tOri);
    frRateHz = expt(iexp).frame_rate;

    % analysis params
    pre_event_frames = frRateHz;
    post_event_frames = frRateHz;
    
    % mouse info
    ms(iexp).subnum = subNum;
    ms(iexp).date = expDate;
    
    % trial type index
    [delayType_ind, delayTypes] = findgroups(audStimDelayType);
    nTrType = length(delayTypes)+2;
    
    % when does vis stim start?
    tVisStimStart = NaN(1,ntrials);
    tVisStimStart(delayType_ind > 1) = tTargetOn(delayType_ind > 1);
    tVisStimStart(delayType_ind == 1 | isnan(delayType_ind)) = tTrStart(delayType_ind == 1 | isnan(delayType_ind));
        
    % auditory only trials
    audOnly_ind = audStimDelayType >= 1000;
    tAudStimStart = tTrStart(audOnly_ind);
    
    %% dF/F
    
    dff_vis = zeros(pre_event_frames+post_event_frames,nc,ntrials);
    for itrial = 1:ntrials
        f0_ind = tTrStart(itrial)-(pre_event_frames-1):tTrStart(itrial);
        vis_ind = tVisStimStart(itrial)-(pre_event_frames-1):tVisStimStart(itrial)+post_event_frames;
        f0 = mean(data(f0_ind,:),1);
        f1 = data(vis_ind,:);
        dff_vis(:,:,itrial) = bsxfun(@rdivide, (bsxfun(@minus,f1,f0)), f0);
    end
    
    dff_audOnly = zeros(pre_event_frames+post_event_frames,nc,sum(audOnly_ind));
    for itrial = 1:length(tAudStimStart)
        f0_ind = tAudStimStart(itrial)-(pre_event_frames-1):tAudStimStart(itrial);
        aud_ind = tAudStimStart(itrial)-(pre_event_frames-1):tAudStimStart(itrial)+post_event_frames;
        f0 = mean(data(f0_ind,:),1);
        f1 = data(aud_ind,:);
        dff_audOnly(:,:,itrial) = bsxfun(@rdivide, (bsxfun(@minus,f1,f0)), f0);
    end
    
    %% sort trial types
    delayTrials = cell(1,nTrType);
    trType_str = cell(1,nTrType);
    trOri = cell(1,nTrType);
    for idelay = 1:nTrType
        if idelay == nTrType
            delayTrials{idelay} = dff_audOnly;
            trType_str{idelay} = 'aud only';
            trOri{idelay} = tOri(audOnly_ind);
        elseif idelay == nTrType - 1
            ind = isnan(delayType_ind);
            delayTrials{idelay} = dff_vis(:,:,ind);
            trType_str{idelay} = 'vis only';
            trOri{idelay} = tOri(ind);
        else
            ind = delayType_ind == idelay;
            delayTrials{idelay} = dff_vis(:,:,ind);
            trType_str{idelay} = num2str(delayTypes(idelay));
            trOri{idelay} = tOri(ind);
        end
    end
    
    ms(iexp).delayTrials = delayTrials;
    ms(iexp).audDelayType = trType_str;
    ms(iexp).trOri = trOri;
    ms(iexp).orientations = unique(tOri);
    ms(iexp).frRateHz = frRateHz;
    ms(iexp).pre_event_frames = pre_event_frames;

end

%% save 
fnout = fullfile(rc.ashleyAnalysis,'Expt Summaries',['awData_audMod' ds]);
save(fullfile(fnout,['awData_audMod' ds '_CaSummary']),'ms');