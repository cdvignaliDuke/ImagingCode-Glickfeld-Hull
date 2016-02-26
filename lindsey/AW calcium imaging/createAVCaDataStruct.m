function mouse = createAVCaDataStruct(datasetStr);
%collects all data from each experiment according to dataset type

    % set analysis windows
    pre_event_time = 1000;
    post_event_time = 3000; 
    pre_win_time = [-30 70];
    trans_win_time = [150 250];
    minTrialLengthMs = 2500; % time in ms of trial rather than hard-coding cycle number for analysis
    
    %load experiment parameters
    
    eval(['awFSAVdatasets' datasetStr])
    av = behavParamsAV;
    rc = behavConstsAV;
    if strcmp(rc.name,'ashley')
        dataGroup = ['awFSAVdatasets' datasetStr];
    else
        dataGroup = [];
    end

    %create list of mice for file names
    nMice = length(unique({expt.SubNum}));
    str = unique({expt.SubNum});
    values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
    mouse_str = ['i' strjoin(str,'_i')];
    %initialize structure
    s = zeros(1,nMice);
    mouse = struct;
    Dirs_all = [];
    
    %collect data from each experiment
    for iexp = 1:size(expt,2)
        disp([num2str(expt(iexp).date) ' i' num2str(expt(iexp).SubNum)])

        nrun = size(expt(iexp).runs,1);
        dir_run = expt(iexp).dirtuning;
        
        %keep track of which mouse
        imouse = find(values == str2num(expt(iexp).SubNum));
        s(:,imouse) = s(:,imouse)+1;
        mouse(imouse).expt(s(:,imouse)).date = expt(iexp).date;
        mouse(imouse).expt(s(:,imouse)).mouse_name = expt(iexp).SubNum;
        
        %translate time windows into frames 
        pre_event_frames = ceil(pre_event_time*(expt(iexp).frame_rate/1000));
        post_event_frames = ceil(post_event_time*(expt(iexp).frame_rate/1000));
        pre_win_frames = pre_event_frames+round(pre_win_time.*(expt(iexp).frame_rate/1000));
        pre_win = pre_win_frames(1):pre_win_frames(2);
        trans_win_frames = pre_event_frames+round(trans_win_time.*(expt(iexp).frame_rate/1000));
        trans_win = trans_win_frames(1):trans_win_frames(2);
        minTrialLengthFrames = ceil(minTrialLengthMs*(expt(iexp).frame_rate/1000));
                
        %create string for saving mult runs
        runstr = expt(iexp).runs(1,:);
        if nrun>1
            for irun = 2:nrun
                runstr = [runstr '-' expt(iexp).runs(irun,:)];
            end
        end
        fnout = fullfile(rc.caOutputDir, expt(iexp).mouse, expt(iexp).folder, expt(iexp).date, [expt(iexp).date '_' expt(iexp).mouse '_' runstr '_']);
        
        % load and combine mworks data and timecourses
        input = [];
        for irun = 1:nrun
            time = expt(iexp).time_mat(irun,:);
            fn_mworks = [rc.pathStr '\data-i' expt(iexp).SubNum '-' expt(iexp).date '-' time '.mat'];
            if irun == 1
                input = mwLoadData(fn_mworks, [], []);
            else
                try
                    input = [input mwLoadData(fn_mworks, [], [])];
                catch
                    input2 = mwLoadData(fn_mworks, [], []);
                    inpNames1 = fieldnames(input);
                    inpNames2 = fieldnames(input2);
                    inpLong = gt(length(inpNames1),length(inpNames2));
                    if inpLong == 1
                        inpPlusInd = ismember(inpNames1,inpNames2);
                        inpPlus = inpNames1(~inpPlusInd);
                        for i = 1:length(inpPlus)
                            input2.(genvarname(inpPlus(i))) = cell(1,80);
                        end
                    else
                        inpPlusInd = ismember(inpNames2,inpNames1);
                        inpPlus = inpNames2(~inpPlusInd);
                        for i = 1:length(inpPlus)
                            input.(char(genvarname(inpPlus(i)))) = cell(1,80);
                        end
                    end
                    input = [input input2];
                end
            end
        end
        input = concatenateDataBlocks(input);        
        
        %convert important fields to matrices
        run_trials = input.trialsSinceReset;
        cLeverDown = celleqel2mat_padded(input.cLeverDown);
        cLeverUp = celleqel2mat_padded(input.cLeverUp);
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cStimOn = celleqel2mat_padded(input.cStimOn);
        cItiStart = celleqel2mat_padded(input.cItiStart);
        reactTimes = celleqel2mat_padded(input.reactTimesMs);
        if expt(iexp).catch
        cCatchOn = celleqel2mat_padded(input.cCatchOn);
        isFA = celleqel2mat_padded(input.tFalseAlarm);
        tCatchGratingDirectionDeg = chop(celleqel2mat_padded(input.tCatchGratingDirectionDeg),4);
        cDirs = unique(tCatchGratingDirectionDeg(~isnan(tCatchGratingDirectionDeg)));
        end
        tooFastTime = input.nFramesTooFast;
        maxReactTime = input.nFramesReact;
        


        
        %account for accumulation of frames across multiple runs 
        dataTC = [];
        offset = 0;
        for irun = 1:nrun
            ImgFolder = expt(iexp).runs(irun,:);
            fnTC = fullfile(rc.ashleyAnalysis, expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, ImgFolder);
            cd(fnTC);
            load('Timecourses.mat')
            dataTC = cat(1, dataTC, dataTimecourse.dataTCsub);
            offset = offset+size(dataTimecourse.dataTCsub,1);
            if irun < nrun
                startTrial = sum(run_trials(1, 1:irun),2)+1;
                endTrial = sum(run_trials(1,1:irun+1),2);
                cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
                cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
                cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
                cStimOn(1,startTrial:endTrial) = cStimOn(1,startTrial:endTrial)+offset;
                cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
                if expt(iexp).catch
                    cCatchOn(1,startTrial:endTrial) = cCatchOn(1,startTrial:endTrial)+offset;
                end
            end
        end

        ntrials = length(input.trialOutcomeCell);
        tCyclesOn = cell2mat(input.tCyclesOn);
        minCyclesOn = input.minCyclesOn;
        cycles = unique(tCyclesOn);
        V_ind = find(cell2mat(input.tBlock2TrialNumber) == 0);
        AV_ind = find(cell2mat(input.tBlock2TrialNumber) == 1);
        cycTime = input.nFramesOn+input.nFramesOff;
        minTime = minCyclesOn*cycTime;
        frameratems = expt(iexp).frame_rate/1000;
        cycTimeMs = cycTime/frameratems;
        minCyclesAnt = floor(minTrialLengthFrames/cycTime);
        
        % for older datasets that do not have catch trial outcome in input
        disp(sum(cell2mat(input.tShortCatchTrial)))
        if expt(iexp).catch
            catchIndex = find(cell2mat(input.tShortCatchTrial));
            if ~any(strcmp(fieldnames(input),'catchTrialOutcomeCell'))
                catchTrialOutcome = cell(1,ntrials);
                for i = 1:sum(cell2mat(input.tShortCatchTrial))
                if isFA(catchIndex(i)) == 1
                    catchTrialOutcome{catchIndex(i)} = 'FA';
                elseif cCatchOn(catchIndex(i)) == 0
                    catchTrialOutcome{catchIndex(i)} = 'failure';
                elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) < tooFastTime
                    catchTrialOutcome{catchIndex(i)} = 'failure';
                elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) > maxReactTime
                    catchTrialOutcome{catchIndex(i)} = 'CR';
                end
                end
            else
                if length(input.catchTrialOutcomeCell) < ntrials
                    catchTrialOutcome = cell(1,ntrials);
                    for i = 1:sum(cell2mat(input.tShortCatchTrial))
                    if isFA(catchIndex(i)) == 1
                        catchTrialOutcome{catchIndex(i)} = 'FA';
                    elseif isnan(cCatchOn(catchIndex(i)))
                        catchTrialOutcome{catchIndex(i)} = 'failure';
                    elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) < tooFastTime
                        catchTrialOutcome{catchIndex(i)} = 'failure';
                    elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) > maxReactTime
                        catchTrialOutcome{catchIndex(i)} = 'CR';
                    end
                    end
                else
                    catchTrialOutcome = input.catchTrialOutcomeCell;
                end
            end
        end

        tGratingDirectionDeg = chop(celleqel2mat_padded(input.tGratingDirectionDeg),4);
        Dirs = unique(tGratingDirectionDeg);
        Dirs_all = unique([Dirs_all Dirs]);
        clear dataTimecourse
        
        %load direction tuning data
        cellSetsLG
        
        %sort by trial type
        if expt(iexp).catch
            Ix = find(cell2mat(input.tShortCatchTrial)==0);
            CIx = find(cell2mat(input.tShortCatchTrial)==1);
        else
            Ix = 1:length(input.trialOutcomeCell);
        end
        
        FIx = intersect(Ix, find(strcmp(input.trialOutcomeCell, 'failure')));
        SIx = intersect(intersect(Ix, find(reactTimes>200)), find(strcmp(input.trialOutcomeCell, 'success')));
        MIx = intersect(Ix, find(strcmp(input.trialOutcomeCell, 'ignore')));
        FIxlong = intersect(find(tCyclesOn>3), FIx);
        SIxlong = intersect(find(tCyclesOn>3), SIx);
        MIxlong = intersect(find(tCyclesOn>3), MIx);
        Fb1Ix = intersect(V_ind, FIxlong);
        Fb2Ix = intersect(AV_ind, FIxlong);
        Sb1Ix = intersect(V_ind, SIxlong);
        Sb2Ix = intersect(AV_ind, SIxlong);
        Mb1Ix = intersect(V_ind, MIxlong);
        Mb2Ix = intersect(AV_ind, MIxlong);
        Rb1Ix = intersect(Ix,intersect(V_ind, find(tCyclesOn>3)));
        Rb2Ix = intersect(Ix,intersect(AV_ind, find(tCyclesOn>3)));

        %find direction with maximum hits and misses to match datasets
        Sb1IxMatch = [];
        Mb1IxMatch = [];
        nS = zeros(1,length(Dirs));
        nM = zeros(1,length(Dirs));
        for iDir = 1:length(Dirs)
            dirIx = find(tGratingDirectionDeg==Dirs(iDir));
            Sb1IxDir = intersect(dirIx,Sb1Ix);
            Mb1IxDir = intersect(dirIx,Mb1Ix);
            nS(iDir) = length(Sb1IxDir);
            nM(iDir) = length(Mb1IxDir);
            if nS(iDir)<nM(iDir)
                Sb1IxMatch = [Sb1IxMatch Sb1IxDir];
                Mb1IxMatch = [Mb1IxMatch Mb1IxDir(randperm(nM(iDir),nS(iDir)))];
            elseif nS(iDir)>nM(iDir)
                Mb1IxMatch = [Mb1IxMatch Mb1IxDir];
                Sb1IxMatch = [Sb1IxMatch Sb1IxDir(randperm(nS(iDir),nM(iDir)))];
            else
                Sb1IxMatch = [Sb1IxMatch Sb1IxDir];
                Mb1IxMatch = [Mb1IxMatch Mb1IxDir];
            end
        end
        if isfield(input, 'tSoundTargetAmplitude')
            tSoundTargetAmplitude = celleqel2mat_padded(input.tSoundTargetAmplitude)';
            Amps = unique(chop(tSoundTargetAmplitude,2));
            nS = zeros(1,length(Amps));
            nM = zeros(1,length(Amps));
            Sb2IxMatch = [];
            Mb2IxMatch = [];
            if length(Amps)>2
                for iAmp = 2:length(Amps)
                    ampIx = find(tSoundTargetAmplitude==Amps(iAmp));
                    Sb2IxAmp = intersect(ampIx,Sb2Ix);
                    Mb2IxAmp = intersect(ampIx,Mb2Ix);
                    nS(iAmp) = length(Sb2IxAmp);
                    nM(iAmp) = length(Mb2IxAmp);
                    if (nS(iAmp)+nM(iAmp))>0
                        if nS(iAmp)<nM(iAmp)
                            Sb2IxMatch = [Sb2IxMatch; Sb2IxAmp];
                            Mb2IxMatch = [Mb2IxMatch; Mb2IxAmp(randperm(nM(iAmp),nS(iAmp)))];
                        elseif nS(iAmp)>nM(iAmp)
                            Mb2IxMatch = [Mb2IxMatch; Mb2IxAmp];
                            Sb2IxMatch = [Sb2IxMatch; Sb2IxAmp(randperm(nS(iAmp),nM(iAmp)))];
                        else
                            Sb2IxMatch = [Sb2IxMatch; Sb2IxAmp];
                            Mb2IxMatch = [Mb2IxMatch; Mb2IxAmp];
                        end
                    end
                end
            else
                Sb2IxMatch = Sb2Ix;
                Mb2IxMatch = Mb2Ix;
            end
        else
            Amps = [0 input.soundTargetAmplitude];
            Sb2IxMatch = Sb2Ix;
            Mb2IxMatch = Mb2Ix;
        end
                
        
        %names of fields
        mouse(imouse).expt(s(:,imouse)).pre_event_frames = pre_event_frames;
        mouse(imouse).expt(s(:,imouse)).post_event_frames = post_event_frames;
        mouse(imouse).expt(s(:,imouse)).visTargets = Dirs;
        if expt(iexp).catch
            mouse(imouse).expt(s(:,imouse)).catchTargets = cDirs;
        end
        mouse(imouse).expt(s(:,imouse)).audTargets = Amps;
        for iDir = 1:length(Dirs)
            mouse(imouse).expt(s(:,imouse)).target(iDir).name = Dirs(iDir);
            mouse(imouse).expt(s(:,imouse)).target(iDir).ind = find(tGratingDirectionDeg==Dirs(iDir));
            if expt(iexp).catch
            mouse(imouse).expt(s(:,imouse)).target(iDir).cind = find(tCatchGratingDirectionDeg==Dirs(iDir));
            end
        end
        mouse(imouse).expt(s(:,imouse)).cells(1).name = 'resp';
        mouse(imouse).expt(s(:,imouse)).cells(2).name = '0';
        mouse(imouse).expt(s(:,imouse)).cells(3).name = '45';
        mouse(imouse).expt(s(:,imouse)).cells(4).name = '90';
        mouse(imouse).expt(s(:,imouse)).cells(5).name = '135';
        mouse(imouse).expt(s(:,imouse)).cells(6).name = 'base_excit';
        mouse(imouse).expt(s(:,imouse)).cells(7).name = 'base_inhib';
        mouse(imouse).expt(s(:,imouse)).cells(8).name = 'tar_excit';
        mouse(imouse).expt(s(:,imouse)).cells(9).name = 'tar_inhib';
        mouse(imouse).expt(s(:,imouse)).cells(10).name = 'base1_resp';
        mouse(imouse).expt(s(:,imouse)).cells(11).name = 'tar1_resp';
        mouse(imouse).expt(s(:,imouse)).cells(12).name = 'all cells';
        mouse(imouse).expt(s(:,imouse)).cells(1).ind = [];
        mouse(imouse).expt(s(:,imouse)).cells(2).ind = cellsSelect{1};
        mouse(imouse).expt(s(:,imouse)).cells(3).ind = cellsSelect{2};
        mouse(imouse).expt(s(:,imouse)).cells(4).ind = cellsSelect{3};
        mouse(imouse).expt(s(:,imouse)).cells(5).ind = cellsSelect{4};
        mouse(imouse).expt(s(:,imouse)).cells(6).ind = [];
        mouse(imouse).expt(s(:,imouse)).cells(7).ind = [];
        mouse(imouse).expt(s(:,imouse)).cells(8).ind = [];
        mouse(imouse).expt(s(:,imouse)).cells(9).ind = [];
        mouse(imouse).expt(s(:,imouse)).align(1).name = 'press';
        mouse(imouse).expt(s(:,imouse)).align(2).name = 'prevStim';
        mouse(imouse).expt(s(:,imouse)).align(3).name = 'catchStim';
        mouse(imouse).expt(s(:,imouse)).win(1).name = 'pre';
        mouse(imouse).expt(s(:,imouse)).win(2).name = 'trans';
        mouse(imouse).expt(s(:,imouse)).win(1).frames = pre_win;
        mouse(imouse).expt(s(:,imouse)).win(2).frames = trans_win;
        mouse(imouse).expt(s(:,imouse)).info.cyc_time = cycTime;
        mouse(imouse).expt(s(:,imouse)).info.cyc_time_ms = cycTimeMs;
        mouse(imouse).expt(s(:,imouse)).info.nCells = size(dataTC,2);
        mouse(imouse).expt(s(:,imouse)).info.isCatch = expt(iexp).catch;
        mouse(imouse).expt(s(:,imouse)).info.minTrialLengthMs = minTrialLengthMs;
        mouse(imouse).expt(s(:,imouse)).info.minTrialLengthFrames = minTrialLengthFrames;
        mouse(imouse).info.allDirs = Dirs_all;
        
        for i = 1:5
            mouse(imouse).expt(s(:,imouse)).align(i).av(1).name = 'visual';
            mouse(imouse).expt(s(:,imouse)).align(i).av(2).name = 'auditory';
            for ii = 1:2
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).name = 'hit';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).name = 'miss';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).name = 'FA';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).name = 'CR';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).name = 'hit_match';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(6).name = 'miss_match';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).n = eval(['length(Sb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).n = eval(['length(Mb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).n = eval(['length(Fb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).n = eval(['length(Rb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).n = eval(['length(Sb' num2str(ii) 'IxMatch)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(6).n = eval(['length(Mb' num2str(ii) 'IxMatch)']);
            end
        end
        
        
        %% Align data to lever down
        ialign = 1;
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataF = zeros(1,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cLeverDown+post_event_frames-1 <  size(dataTC,1)),[],2)
            Data(:,:,itrial) = dataTC(cLeverDown(itrial)-pre_event_frames:cLeverDown(itrial)+post_event_frames-1,:);
            DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
        end
        
        %identify trials with motion (large peaks in the derivative)
        sz = size(DataDFoverF);
        ind_motion = find(max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>0.05);
        mouse(imouse).expt(s(:,imouse)).align(ialign).ind_motion = ind_motion;
        
        %divide data by trial type and outcome
        for iav = 1:2
            if length(eval(['Sb' num2str(iav) 'Ix']))>0
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));        
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = DataDFoverF(:,:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)); 
            else
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = NaN(sz(1), sz(2));
               
            end
        end
        
        %find trials of a minimum trail length and divide up by outcome
        tCycInd = find(tCyclesOn >= minCyclesAnt);
        for iav = 1:2
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Sb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
        end            
        
        %test for significantly responsive cells
        %to first baseline stimulus (ttest of means)
        [h, p] = ttest(squeeze(mean(DataDFoverF(pre_win,:,setdiff(SIx,ind_motion)),1)), squeeze(mean(DataDFoverF(trans_win,:,setdiff(SIx,ind_motion)),1)), 'dim', 2, 'tail', 'left', 'alpha', 0.05/(length(Dirs)));
        mouse(imouse).expt(s(:,imouse)).align(1).ttest_trans = h;
        %to baseline stimulus before earliest target (ttest of windows across trials)
        [h, p] = ttest(squeeze(mean(DataDFoverF(1:pre_event_frames,:,setdiff(SIx,ind_motion)),1)), squeeze(mean(DataDFoverF(minTime+1:pre_event_frames+minTime,:,setdiff(SIx,ind_motion)),1)), 'dim', 2, 'alpha', 0.05/2);
        mouse(imouse).expt(s(:,imouse)).align(1).ttest_sust = h;
        baseStimRespDiff = squeeze(mean(mean(DataDFoverF(minTime+1:pre_event_frames+minTime,:,setdiff(SIx,ind_motion)),3),1)) - squeeze(mean(mean(DataDFoverF(1:pre_event_frames,:,setdiff(SIx,ind_motion)),3),1));
        %to baseline stimulus before earliest target (vartest of mean response)
        [h, p] = vartest2(squeeze(mean(DataDFoverF(1:pre_event_frames,:,setdiff(SIx,ind_motion)),3)), squeeze(mean(DataDFoverF(minTime+1:pre_event_frames+minTime,:,setdiff(SIx,ind_motion)),3)), 'dim', 2, 'tail', 'left', 'alpha', 0.05/2);
        mouse(imouse).expt(s(:,imouse)).align(1).vartest = h;
        %% Align data to previous stim (could be target or stim before lever release on FA)
        ialign = 2;
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        Data_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cTargetOn+post_event_frames-1 <  size(dataTC,1)),[],2)
            if strcmp(input.trialOutcomeCell(itrial), 'failure')
                Data(:,:,itrial) = dataTC(cStimOn(itrial)-pre_event_frames:cStimOn(itrial)+post_event_frames-1,:);
                Data_CR(:,:,itrial) = dataTC(cStimOn(itrial)-pre_event_frames-cycTime:cStimOn(itrial)+post_event_frames-cycTime-1,:);
            else
                Data(:,:,itrial) = dataTC(cTargetOn(itrial)-pre_event_frames:cTargetOn(itrial)+post_event_frames-1,:);
                Data_CR(:,:,itrial) = dataTC(cTargetOn(itrial)-pre_event_frames-cycTime:cTargetOn(itrial)+post_event_frames-cycTime-1,:);
            end
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), DataF(:,:,itrial));
            DataDF_CR(:,:,itrial) = bsxfun(@minus, Data_CR(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF_CR(:,:,itrial) = bsxfun(@rdivide, DataDF_CR(:,:,itrial), DataF(:,:,itrial));
        end
        
        %identify trials with motion (large peaks in the derivative)
        ind_motion = find(max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>0.05);
        mouse(imouse).expt(s(:,imouse)).align(ialign).ind_motion = ind_motion;
        
        %divide data by trial type and outcome
        for iav = 1:2
            if length(eval(['Sb' num2str(iav) 'Ix']))>0
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = DataDFoverF(:,:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).resp = DataDFoverF_CR(:,:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tcyc = tCyclesOn(:,setdiff(eval(['Sb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = tCyclesOn(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = tCyclesOn(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tcyc = tCyclesOn(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
            else
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tcyc = NaN;
            end
        end
        
        %divide data by target direction
        mouse(imouse).expt(s(:,imouse)).align(ialign).ttest_trans = zeros(size(dataTC,2), length(Dirs));
        tarStimRespDiff = zeros(size(dataTC,2), length(Dirs));
        for iDir = 1:length(Dirs)
            iav = 1;
            vInd = setdiff(find(tGratingDirectionDeg==Dirs(iDir)), ind_motion);
            %if iDir == 1, then the direction is 0 and it's an auditory trial
            if iDir == 1
                indS = intersect(vInd, Sb2Ix);
                indM = intersect(vInd, Mb2Ix);
                indSM = intersect(vInd, Sb2IxMatch);
                indMM = intersect(vInd, Mb2IxMatch);
            else
                indS = intersect(vInd, Sb1Ix);
                indM = intersect(vInd, Mb1Ix);
                indSM = intersect(vInd, Sb1IxMatch);
                indMM = intersect(vInd, Mb1IxMatch);
            end
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iDir} = DataDFoverF(:,:,indS);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).stimResp{iDir} = DataDFoverF(:,:,indM);        
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).stimResp{iDir} = DataDFoverF(:,:,indSM);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).stimResp{iDir} = DataDFoverF(:,:,indMM);
            %test for responsiveness to any one of the presented directions
            [h, p] = ttest(squeeze(mean(DataDFoverF(pre_win,:,[indS indM]),1)), squeeze(mean(DataDFoverF(trans_win,:,[indS indM]),1)), 'dim', 2, 'tail', 'left', 'alpha', 0.05/(length(Dirs)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).ttest_trans(:,iDir) = h;
            tarStimRespDiff(:,iDir) = squeeze(mean(mean(DataDFoverF(trans_win,:,[indS indM]),3),1)) - squeeze(mean(mean(DataDFoverF(pre_win,:,[indS indM]),3),1));
        end
        
        %divide data by target sound amplitude
        for iAmp = 1:length(Amps)
            iav = 2;
            if isfield(input, 'tSoundTargetAmplitude')
                aInd = setdiff(find(tSoundTargetAmplitude==Amps(iAmp)), ind_motion)';
            else
                aInd = setdiff(find(celleqel2mat_padded(input.tDoAuditoryDetect)), ind_motion);
            end
            %if iDir == 1, then the amplitude is 0 and it's a visual trial
            indS = NaN; indM = NaN; indSM = NaN; indMM = NaN;
            if iAmp == 1
                indS = intersect(aInd, Sb1Ix);
                indM = intersect(aInd, Mb1Ix);
                indSM = intersect(aInd, Sb1IxMatch);
                indMM = intersect(aInd, Mb1IxMatch);
            else
                indS = intersect(aInd, Sb2Ix);
                indM = intersect(aInd, Mb2Ix);
                indSM = intersect(aInd, Sb2IxMatch);
                indMM = intersect(aInd, Mb2IxMatch);
            end
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iAmp} = DataDFoverF(:,:,indS);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).stimResp{iAmp} = DataDFoverF(:,:,indM);        
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).stimResp{iAmp} = DataDFoverF(:,:,indSM);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).stimResp{iAmp} = DataDFoverF(:,:,indMM);
        end
        
        %% Align data to vis catch stim
        if expt(iexp).catch
        ialign = 3;
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cCatchOn+post_event_frames-1 <  size(dataTC,1)),[],2)
            if expt(iexp).catch
            if strcmp(catchTrialOutcome{itrial}, 'failure') | isempty(catchTrialOutcome{itrial})
                Data(:,:,itrial) = NaN(pre_event_frames+post_event_frames,size(dataTC,2));
            else
                Data(:,:,itrial) = dataTC(cCatchOn(itrial)-pre_event_frames:cCatchOn(itrial)+post_event_frames-1,:);
            end
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), DataF(:,:,itrial));
            end
        end
        
        %identify trials with motion (large peaks in the derivative)
        ind_motion = find(max(diff(squeeze(nanmean(DataDFoverF,2)),1),[],1)>0.05);
        disp(length(ind_motion));
        mouse(imouse).expt(s(:,imouse)).align(ialign).ind_motion = ind_motion;
        faIx = find(strcmp(catchTrialOutcome,'FA'));
        crIx = find(strcmp(catchTrialOutcome,'CR'));
        
        %divide data by trial type and outcome
            iav = 1;
            hitsMatched2Catch = [];
            missMatched2Catch = [];
            if length([faIx crIx])>0
                hitsMatched2Catch = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(1).stimResp(find(ismember(Dirs,cDirs)));
                missMatched2Catch = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(2).stimResp(find(ismember(Dirs,cDirs)));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = cat(3,hitsMatched2Catch{:});
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = cat(3,missMatched2Catch{:});
                 % FA in this case means responded to catch trial
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = DataDFoverF(:,:,setdiff(faIx,ind_motion));
                 % CR in this case means ignored catch trial
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).resp = DataDFoverF(:,:,setdiff(crIx,ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = tCyclesOn(:,setdiff(faIx,ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = tCyclesOn(:,setdiff(crIx,ind_motion));
            else
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = NaN;
            end
                
        %divide data by catch direction
        
        cind = [];
        for iDir = 1:length(Dirs)
            cind = find(tCatchGratingDirectionDeg == Dirs(iDir));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iDir} = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(1).stimResp{iDir};
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iDir} = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(2).stimResp{iDir};
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).stimResp{iDir} = DataDFoverF(:,:,setdiff(intersect(faIx,cind),ind_motion));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).stimResp{iDir} = DataDFoverF(:,:,setdiff(intersect(crIx,cind),ind_motion));
        end
        end
        %% find cells that are responsive to at least one of the target
        %directions or the first baseline visual stimulus 
        base_ind = find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_trans);
        resp_ind = find(sum(mouse(imouse).expt(s(:,imouse)).align(2).ttest_trans,2)>=1);
        mouse(imouse).expt(s(:,imouse)).cells(1).ind = unique([base_ind; resp_ind]);
        mouse(imouse).expt(s(:,imouse)).cells(6).ind = intersect(find(baseStimRespDiff > 0),(find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_sust)));
        mouse(imouse).expt(s(:,imouse)).cells(7).ind = intersect(find(baseStimRespDiff < 0),(find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_sust)));
        mouse(imouse).expt(s(:,imouse)).cells(8).ind = intersect(find(tarStimRespDiff > 0),resp_ind);
        mouse(imouse).expt(s(:,imouse)).cells(9).ind = intersect(find(tarStimRespDiff < 0),resp_ind);
        mouse(imouse).expt(s(:,imouse)).cells(10).ind = base_ind;
        mouse(imouse).expt(s(:,imouse)).cells(11).ind = resp_ind;
        mouse(imouse).expt(s(:,imouse)).cells(12).ind = unique(cat(1,mouse(imouse).expt(s(:,imouse)).cells(1).ind,mouse(imouse).expt(s(:,imouse)).cells(7).ind,mouse(imouse).expt(s(:,imouse)).cells(7).ind));
    end
    
    save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_CaSummary' datasetStr '.mat']), 'mouse','-v7.3');
end

        
        
        

