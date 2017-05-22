function mouse = createAudControlCaDataStruct(datasetStr,cellsOnly);
%collects all data from each experiment according to dataset type
%set cellsOnly to 1 to grab data from cells only masks (rather than cells
%and dendrites)

    % set analysis windows
    pre_event_time = 1000;
    post_event_time = 4000; 
    pre_win_time = [-30 70];
    trans_win_time = [150 250]; %this is actually ~166-266 ms at 30 Hz
    minTrialLengthMs = 2500; % time in ms of trial rather than hard-coding cycle number for analysis
    
    %load experiment parameters

    eval(['awFSAVdatasets' datasetStr]);
    rc = behavConstsAV;
    if strcmp(rc.name,'ashle')
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
    
    set(0,'defaultfigurepaperorientation','portrait');
    set(0,'defaultfigurepapersize',[8.5 11]);
    set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
    nexp = size(expt,2);
    n = ceil(sqrt(nexp+1));
    if (n^2)-n > nexp+1
        n2 = n-1;
    else
        n2 = n;
    end
    motionHist = figure;
    
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
                            input2.(genvarname(inpPlus{i})) = cell(1,input2.trialSinceReset);
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
        tooFastTime = input.nFramesTooFast;
        maxReactTime = input.nFramesReact;
        


        
        %account for accumulation of frames across multiple runs 
        dataTC = [];
        fnTC = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,expt(iexp).folder, expt(iexp).date);
        if cellsOnly == 2
            load(fullfile(fnTC,'timecourses_dendrites.mat'))
        else
            load(fullfile(fnTC,'timecourses.mat'))
        end
        dataTC = data_tc_subnp;
                
        offset = 0;
        for irun = 1:nrun
            ImgFolder = expt(iexp).runs(irun,:);
            nfr_run = nFramesSbxDataset(expt(iexp).mouse,expt(iexp).date,ImgFolder);
            offset = offset+nfr_run;
            if irun < nrun
                startTrial = sum(run_trials(1, 1:irun),2)+1;
                endTrial = sum(run_trials(1,1:irun+1),2);
                cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
                cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
                cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
                cStimOn(1,startTrial:endTrial) = cStimOn(1,startTrial:endTrial)+offset;
                cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
            end
        end

        ntrials = length(input.trialOutcomeCell);
        tCyclesOn = cell2mat(input.tCyclesOn);
        minCyclesOn = input.minCyclesOn;
        cycles = unique(tCyclesOn);
        V_ind = find(cell2mat(input.tBlock2TrialNumber) == 0);
        AV_ind = find(cell2mat(input.tBlock2TrialNumber) == 1);
        
        %previous trial info
        trType = double(cell2mat(input.tBlock2TrialNumber));
        trType_shift = [NaN trType];
        prevTrType = trType_shift(1:length(trType));
        trType_shift = [NaN NaN trType];
        prev2TrType = trType_shift(1:length(trType));
        trOut = input.trialOutcomeCell;
        trOut_shift = [{NaN} trOut];
        prevTrOut = trOut_shift(1:length(trOut));
        trOut_shift = [{NaN} {NaN} trOut];
        prev2TrOut = trOut_shift(1:length(trOut));
        
        %stim timing
        if iscell(input.nFramesOn)
            cycTime = unique(cell2mat(input.nFramesOn))+unique(cell2mat(input.nFramesOff));
        else
            cycTime = input.nFramesOn+input.nFramesOff;
        end
        minTime = minCyclesOn*cycTime;
        frameratems = expt(iexp).frame_rate/1000;
        cycTimeMs = cycTime/frameratems;
        minCyclesAnt = floor(minTrialLengthFrames/cycTime);

        tGratingDirectionDeg = chop(celleqel2mat_padded(input.tGratingDirectionDeg),4);
        Dirs = unique(tGratingDirectionDeg);
        Dirs_all = unique([Dirs_all Dirs]);
        clear dataTimecourse
        
        %load direction tuning data
        dataPath = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, expt(iexp).dirtuning);
        if cellsOnly == 1
            load(fullfile(dataPath, 'cellsSelect_cellsOnly.mat'));
        elseif cellsOnly == 2
            load(fullfile(dataPath, 'cellsSelect_dendrites.mat'));
        else
            load(fullfile(dataPath, 'cellsSelect.mat'));
        end
        
        %sort by trial type

        Ix = 1:length(input.trialOutcomeCell);
        
        F1Ix = find(tCyclesOn == 1);
        F1oneStim = intersect(V_ind,F1Ix);
        F2oneStim = intersect(AV_ind,F1Ix);
        
        FIx = intersect(Ix, find(strcmp(input.trialOutcomeCell, 'failure')));
        SIx = intersect(intersect(Ix, find(reactTimes>200)), find(strcmp(input.trialOutcomeCell, 'success')));
        SIxAllReact = intersect(Ix, find(strcmp(input.trialOutcomeCell, 'success')));
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
        SIxAllReactlong = intersect(find(tCyclesOn>3), SIxAllReact);
        
% %         %find result of previous trial
% %         if sum(FIx==1)
% %             if sum(intersect(FIx,V_ind)==1)
% %             Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx(2:end)-1,V_ind)));
% %             Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx(2:end)-1,AV_ind)));
% %             Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,AV_ind)));
% %             else
% %             Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx(2:end)-1,V_ind)));
% %             Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx(2:end)-1,AV_ind)));
% %             end
% %             Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
% %             Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
% %             Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,AV_ind)));
% %             Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
% %             SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %             SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %         elseif sum(SIx==1)            
% %             if sum(intersect(SIx,V_ind)==1)
% %             Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx(2:end)-1,V_ind)));
% %             Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx(2:end)-1,AV_ind)));
% %             Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,AV_ind)));
% %             else
% %             Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
% %             Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx(2:end)-1,V_ind)));
% %             Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx(2:end)-1,AV_ind)));
% %             end
% %             Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
% %             Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
% %             SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %             SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %         elseif sum(MIx==1)
% %             if sum(intersect(MIx,V_ind)==1)
% %             Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx(2:end)-1,V_ind)));
% %             Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx(2:end)-1,AV_ind)));
% %             Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
% %             else
% %             Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
% %             Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx(2:end)-1,V_ind)));
% %             Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx(2:end)-1,AV_ind)));
% %             end
% %             Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
% %             Fb5Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb6Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Sb5Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb6Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
% %             SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %             SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %         elseif sum(SIxAllReactlong==1)
% %             if sum(intersect(SIxAllReactlong,V_ind)==1)
% %             SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong(2:end)-1,V_ind)));
% %             SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong(2:end)-1,AV_ind)));
% %             SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %             else
% %             SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %             SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong(2:end)-1,V_ind)));
% %             SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong(2:end)-1,AV_ind)));                
% %             end
% %             Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
% %             Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
% %             Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,AV_ind)));
% %             Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
% %         else
% %             Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
% %             Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
% %             Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,V_ind)));
% %             Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,AV_ind)));
% %             Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,V_ind)));
% %             Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,AV_ind)));
% %             Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
% %             Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
% %             SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %             SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
% %             SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
% %         end
            
        
        SbAR1Ix = intersect(V_ind, SIxAllReactlong);
        SbAR2Ix = intersect(AV_ind, SIxAllReactlong);

%         %find direction with maximum hits and misses to match datasets
%         Sb1IxMatch = [];
%         Mb1IxMatch = [];
%         nS = zeros(1,length(Dirs));
%         nM = zeros(1,length(Dirs));
%         for iDir = 1:length(Dirs)
%             dirIx = find(tGratingDirectionDeg==Dirs(iDir));
%             Sb1IxDir = intersect(dirIx,Sb1Ix);
%             Mb1IxDir = intersect(dirIx,Mb1Ix);
%             nS(iDir) = length(Sb1IxDir);
%             nM(iDir) = length(Mb1IxDir);
%             if nS(iDir)<nM(iDir)
%                 Sb1IxMatch = [Sb1IxMatch Sb1IxDir];
%                 Mb1IxMatch = [Mb1IxMatch Mb1IxDir(randperm(nM(iDir),nS(iDir)))];
%             elseif nS(iDir)>nM(iDir)
%                 Mb1IxMatch = [Mb1IxMatch Mb1IxDir];
%                 Sb1IxMatch = [Sb1IxMatch Sb1IxDir(randperm(nS(iDir),nM(iDir)))];
%             else
%                 Sb1IxMatch = [Sb1IxMatch Sb1IxDir];
%                 Mb1IxMatch = [Mb1IxMatch Mb1IxDir];
%             end
%         end
        
        %names of fields
        mouse(imouse).expt(s(:,imouse)).pre_event_frames = pre_event_frames;
        mouse(imouse).expt(s(:,imouse)).post_event_frames = post_event_frames;
        mouse(imouse).expt(s(:,imouse)).visTargets = Dirs;
        for iDir = 1:length(Dirs)
            mouse(imouse).expt(s(:,imouse)).target(iDir).name = Dirs(iDir);
            mouse(imouse).expt(s(:,imouse)).target(iDir).ind = find(tGratingDirectionDeg==Dirs(iDir));
        end
        mouse(imouse).expt(s(:,imouse)).cells(1).name = 'resp';
        mouse(imouse).expt(s(:,imouse)).cells(2).name = '0';
        mouse(imouse).expt(s(:,imouse)).cells(3).name = '45';
        mouse(imouse).expt(s(:,imouse)).cells(4).name = '90';
        mouse(imouse).expt(s(:,imouse)).cells(5).name = '135';
        mouse(imouse).expt(s(:,imouse)).cells(6).name = 'untuned';
        mouse(imouse).expt(s(:,imouse)).cells(7).name = 'all tuned';
        mouse(imouse).expt(s(:,imouse)).cells(8).name = 'base_excit';
        mouse(imouse).expt(s(:,imouse)).cells(9).name = 'base_inhib';
        mouse(imouse).expt(s(:,imouse)).cells(10).name = 'tar_excit';
        mouse(imouse).expt(s(:,imouse)).cells(11).name = 'tar_inhib';
        mouse(imouse).expt(s(:,imouse)).cells(12).name = 'base1_resp';
        mouse(imouse).expt(s(:,imouse)).cells(13).name = 'tar1_resp';
        mouse(imouse).expt(s(:,imouse)).cells(14).name = 'all cells';
        mouse(imouse).expt(s(:,imouse)).cells(1).ind = [];
        mouse(imouse).expt(s(:,imouse)).cells(2).ind = cellsSelect{1};
        mouse(imouse).expt(s(:,imouse)).cells(3).ind = cellsSelect{2};
        mouse(imouse).expt(s(:,imouse)).cells(4).ind = cellsSelect{3};
        mouse(imouse).expt(s(:,imouse)).cells(5).ind = cellsSelect{4};
        mouse(imouse).expt(s(:,imouse)).cells(6).ind = cellsSelect{5};
        mouse(imouse).expt(s(:,imouse)).cells(7).ind = cellsSelect{6};
        mouse(imouse).expt(s(:,imouse)).align(1).name = 'press';
        mouse(imouse).expt(s(:,imouse)).align(2).name = 'targetStim';
        mouse(imouse).expt(s(:,imouse)).align(3).name = 'catchStim';
        mouse(imouse).expt(s(:,imouse)).win(1).name = 'pre';
        mouse(imouse).expt(s(:,imouse)).win(2).name = 'trans';
        mouse(imouse).expt(s(:,imouse)).win(1).frames = pre_win;
        mouse(imouse).expt(s(:,imouse)).win(2).frames = trans_win;
        mouse(imouse).expt(s(:,imouse)).info.cyc_time = cycTime;
        mouse(imouse).expt(s(:,imouse)).info.cyc_time_ms = cycTimeMs;
        mouse(imouse).expt(s(:,imouse)).info.nCells = size(dataTC,2);
        mouse(imouse).expt(s(:,imouse)).info.minTrialLengthMs = minTrialLengthMs;
        mouse(imouse).expt(s(:,imouse)).info.minTrialLengthFrames = minTrialLengthFrames;
        mouse(imouse).expt(s(:,imouse)).info.dirs = Dirs;
        mouse(imouse).info.allDirs = Dirs_all;
        mouse(imouse).expt(s(:,imouse)).tuning(1).name = 'ori preference';
        mouse(imouse).expt(s(:,imouse)).tuning(2).name = 'OSI';
        mouse(imouse).expt(s(:,imouse)).tuning(3).name = 'DSI';
        mouse(imouse).expt(s(:,imouse)).tuning(1).outcome = ori_ind_all;
        mouse(imouse).expt(s(:,imouse)).tuning(2).outcome = OSI;
        mouse(imouse).expt(s(:,imouse)).tuning(3).outcome = DSI;
        
        for i = 1:2
            mouse(imouse).expt(s(:,imouse)).align(i).av(1).name = 'visual';
            mouse(imouse).expt(s(:,imouse)).align(i).av(2).name = 'auditory';
            for ii = 1:2
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).name = 'hit';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).name = 'miss';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).name = 'FA';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).name = 'CR';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).name = 'hit_match';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(6).name = 'miss_match';
                if i == 1
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).name = 'unmatched hits - dirs';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).name = 'unmatched miss - dirs';        
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).name = 'miss match FA';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).name = 'miss match CR';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).name = 'FA match hits';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(6).name = 'FA match miss';
                end
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
        motionThreshold = 0.1;
        sz = size(DataDFoverF);
        ind_motion = find(max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold);
        mouse(imouse).expt(s(:,imouse)).align(ialign).ind_motion = ind_motion;
        
        figure(motionHist);
        subplot(n,n2,iexp)
        hist(max(diff(squeeze(mean(DataDFoverF,2)),1),[],1));
        hold on
        vline(motionThreshold,'k:')
        xlim([0 0.1])
        xlabel('max diff')
        ylabel('n trials')
               
        %divide data by trial type and outcome
        for iav = 1:2
            if length(eval(['SbAR' num2str(iav) 'Ix']))>0
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = DataDFoverF(:,:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));        
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = DataDFoverF(:,:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)); 
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).oneStim = DataDFoverF(:,:,setdiff(eval(['F' num2str(iav) 'oneStim']),ind_motion));                 
                 
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tcyc = tCyclesOn(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = tCyclesOn(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = tCyclesOn(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tcyc = tCyclesOn(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).prevTrType = prevTrType(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).prevTrType = prevTrType(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).prevTrType = prevTrType(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).prevTrType = prevTrType(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).prevTrType = prevTrType(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).prevTrType = prevTrType(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));                
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).prev2TrType = prev2TrType(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).prev2TrType = prev2TrType(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).prev2TrType = prev2TrType(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).prev2TrType = prev2TrType(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).prev2TrType = prev2TrType(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).prev2TrType = prev2TrType(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).prevTrOut = prevTrOut(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).prevTrOut = prevTrOut(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).prevTrOut = prevTrOut(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).prevTrOut = prevTrOut(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).prevTrOut = prevTrOut(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).prevTrOut = prevTrOut(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));                
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).prev2TrOut = prev2TrOut(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).prev2TrOut = prev2TrOut(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).prev2TrOut = prev2TrOut(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).prev2TrOut = prev2TrOut(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).prev2TrOut = prev2TrOut(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).prev2TrOut = prev2TrOut(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
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
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion)));
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
        end
        
        %% find cumulative average for each cycle length
        cycs = 1:max(tCyclesOn);
        for icyc = 1:length(cycs)
            minCyc = cycs(icyc);
            tCycInd = find(tCyclesOn >= minCyc);
            oneCycInd = find(tCyclesOn == minCyc);
            trLengthFrames = cycTime*(icyc);

            for iav = 1:2
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion)));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion)));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
            end

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
        ind_motion = find(max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold);
        mouse(imouse).expt(s(:,imouse)).align(ialign).ind_motion = ind_motion;
        
        %divide data by trial type and outcome
        for iav = 1:2
            if length(eval(['Sb' num2str(iav) 'Ix']))>0
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = DataDFoverF(:,:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).reactTimes = reactTimes(setdiff(eval(['Sb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).reactTimes = reactTimes(setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).reactTimes = reactTimes(setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tcyc = tCyclesOn(:,setdiff(eval(['Sb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = tCyclesOn(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).resp = DataDFoverF_CR(:,:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                %react times for each trial
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).reactTimes = reactTimes(setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).reactTimes = reactTimes(setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).reactTimes = reactTimes(setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                %trial length 
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = tCyclesOn(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tcyc = tCyclesOn(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
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
            for iav = 1:2;
            vInd = setdiff(find(tGratingDirectionDeg==Dirs(iDir)), ind_motion);

            indS = intersect(vInd, eval(['Sb' num2str(iav) 'Ix']));%Sb1Ix);
            indM = intersect(vInd, eval(['Mb' num2str(iav) 'Ix']));%Mb1Ix);
%             if iav < 3
%             indSM = intersect(vInd, Sb1IxMatch);
%             indMM = intersect(vInd, Mb1IxMatch);
%             end
            
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iDir} = DataDFoverF(:,:,indS);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).stimResp{iDir} = DataDFoverF(:,:,indM);  
%             if iav < 3
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).stimResp{iDir} = DataDFoverF(:,:,indSM);
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).stimResp{iDir} = DataDFoverF(:,:,indMM);
%             end
            
             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).trialL{iDir} = tCyclesOn(:,indS);
             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).trialL{iDir} = tCyclesOn(:,indM);
             
            %test for responsiveness to any one of the presented directions
            if iav == 1
            [h, p] = ttest(squeeze(mean(DataDFoverF(pre_win,:,[indS indM]),1)), squeeze(mean(DataDFoverF(trans_win,:,[indS indM]),1)), 'dim', 2, 'tail', 'left', 'alpha', 0.05/(length(Dirs)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).ttest_trans(:,iDir) = h;
            tarStimRespDiff(:,iDir) = squeeze(mean(mean(DataDFoverF(trans_win,:,[indS indM]),3),1)) - squeeze(mean(mean(DataDFoverF(pre_win,:,[indS indM]),3),1));
            end
            end
        end
        
        %% Align data to lever down, include target
        ialign = 3;
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
        
        cycs = 1:max(tCyclesOn);
        for icyc = 1:length(cycs)
            minCyc = cycs(icyc);
            oneCycInd = find(tCyclesOn == minCyc);
        for iav = 1:2
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));              
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
            
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tTargetOn{icyc} = cTargetOn(intersect(oneCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tTargetOn{icyc} = cTargetOn(intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tTargetOn{icyc} = cTargetOn(intersect(oneCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tTargetOn{icyc} = cTargetOn(intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
            
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).direction{icyc} = tGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).direction{icyc} = tGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).direction{icyc} = tGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).direction{icyc} = tGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
           
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).name = 'hit';
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).name = 'miss';             
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).name = 'hit - matched';
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).name = 'miss - matched';
            
        end
        end
        %% find cells that are responsive to at least one of the target
        %directions or the first baseline visual stimulus 
        base_ind = find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_trans);
        resp_ind = find(sum(mouse(imouse).expt(s(:,imouse)).align(2).ttest_trans,2)>=1);
        mouse(imouse).expt(s(:,imouse)).cells(1).ind = unique([base_ind; resp_ind]);
        mouse(imouse).expt(s(:,imouse)).cells(8).ind = intersect(find(baseStimRespDiff > 0),(find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_sust)));
        mouse(imouse).expt(s(:,imouse)).cells(9).ind = intersect(find(baseStimRespDiff < 0),(find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_sust)));
        mouse(imouse).expt(s(:,imouse)).cells(10).ind = intersect(find(mean(tarStimRespDiff,2) > 0),resp_ind);
        mouse(imouse).expt(s(:,imouse)).cells(11).ind = intersect(find(mean(tarStimRespDiff,2) < 0),resp_ind);
        mouse(imouse).expt(s(:,imouse)).cells(12).ind = base_ind;
        mouse(imouse).expt(s(:,imouse)).cells(13).ind = resp_ind;
        mouse(imouse).expt(s(:,imouse)).cells(14).ind = unique(cat(1,mouse(imouse).expt(s(:,imouse)).cells(1).ind,mouse(imouse).expt(s(:,imouse)).cells(8).ind));
    end
    
    figure(motionHist)
    if cellsOnly == 1
        save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_CaSummary_cells' datasetStr '.mat']), 'mouse','-v7.3');
    elseif cellsOnly == 2
        save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_CaSummary_dendrites' datasetStr '.mat']), 'mouse','-v7.3');
    else
    print(fullfile(rc.caOutputDir, dataGroup, [datasetStr '_motionHist.pdf']), '-dpdf')
    save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_CaSummary' datasetStr '.mat']), 'mouse','-v7.3');
    end
end
        
        
        

