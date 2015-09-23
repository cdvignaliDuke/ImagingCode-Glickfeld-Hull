function mouse = createAVCaDataStruct(doPlot);
    close all
    awFSAVdatasets
    av = behavParamsAV;
    rc = behavConstsAV;
    pre_event_time = 1000;
    post_event_time = 3000;
    trans_win_time = [100 200];
    mid_win_time = [400 500];
    late_win_time = [1000 1100];
    s = zeros(1,2);
    mouse = struct;
    for iexp = 1:size(expt,2)
        disp(num2str(iexp))
        SubNum = expt(iexp).SubNum;
        date_name = expt(iexp).date;
        runs = expt(iexp).runs;
        time_mat = expt(iexp).time_mat;
        mouse_name = expt(iexp).mouse;
        folder = expt(iexp).folder;
        frame_rate = expt(iexp).frame_rate;
        nrun = size(runs,1);
        dir_run = expt(iexp).dirtuning;
        doCatch = expt(iexp).catch;
        
        str = sprintf('%f,', av.mouse);
        values = textscan(str, '%f', 'delimiter', ',', 'EmptyValue', NaN);
        imouse = find(values{1} == str2num(SubNum));
        s(:,imouse) = s(:,imouse)+1;
        pre_event_frames = ceil(pre_event_time*(frame_rate/1000));
        post_event_frames = ceil(post_event_time*(frame_rate/1000));
        trans_win_frames = pre_event_frames+round(trans_win_time.*(frame_rate/1000));
        trans_win = trans_win_frames(1):trans_win_frames(2);
        mid_win_frames = pre_event_frames+round(mid_win_time.*(frame_rate/1000));
        mid_win = mid_win_frames(1):mid_win_frames(2);
        late_win_frames = pre_event_frames+round(late_win_time.*(frame_rate/1000));
        late_win = late_win_frames(1):late_win_frames(2);
        
        %create string for saving mult runs
        runstr = runs(1,:);
        if nrun>1
            for irun = 2:nrun
                runstr = [runstr '-' runs(irun,:)];
            end
        end
        fnout = fullfile(rc.caOutputDir, mouse_name, date_name, [date_name '_' mouse_name '_' runstr '_']);
        
        % load and combine mworks data and timecourses
        input = [];
        for irun = 1:nrun
            time = time_mat(irun,:);
            fn_mworks = [rc.pathStr '\data-i' SubNum '-' date_name '-' time '.mat'];
            if irun == 1
                input = mwLoadData(fn_mworks, [], []);
            else
                input = [input mwLoadData(fn_mworks, [], [])];
            end
        end
        input = concatenateDataBlocks(input);
        
        run_trials = input.trialsSinceReset;
        cLeverDown = cell2mat(input.cLeverDown);
        cLeverUp = cell2mat(input.cLeverUp);
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cStimOn = celleqel2mat_padded(input.cStimOn);
        cItiStart = cell2mat(input.cItiStart);

        dataTC = [];
        offset = 0;
        for irun = 1:nrun
            ImgFolder = runs(irun,:);
            fnTC = fullfile(rc.ashleyAnalysis, mouse_name,folder, date_name, ImgFolder);
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
            end
        end

        ntrials = length(input.trialOutcomeCell);
        trialOutcome = cell2mat(input.trialOutcomeCell);
        tCyclesOn = cell2mat(input.tCyclesOn);
        cycles = unique(tCyclesOn);
        V_ind = find(cell2mat(input.tBlock2TrialNumber) == 0);
        AV_ind = find(cell2mat(input.tBlock2TrialNumber) == 1);
        cycTime = input.nFramesOn+input.nFramesOff;

        tGratingDirectionDeg = chop(cell2mat(input.tGratingDirectionDeg),4);
        Dirs = unique(tGratingDirectionDeg);
        clear dataTimecourse
        
        %load direction tuning data
        cellSetsLG
        
        %sort by trial type
        if doCatch
            Ix = find(cell2mat(input.isCatchTrial)==0);
        else
            Ix = 1:length(input.trialOutcomeCell);
        end
        
        FIx = intersect(Ix, find(strcmp(input.trialOutcomeCell, 'failure')));
        SIx = intersect(Ix, find(strcmp(input.trialOutcomeCell, 'success')));
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
        
        %find direction with maximum hits and misses
        Sb1IxMatch = [];
        Mb1IxMatch = [];
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
            tSoundTargetAmplitude = cell2mat_padded(input.tSoundTargetAmplitude)';
            Amps = unique(tSoundTargetAmplitude);
            Sb1IxMatch = [];
            Mb1IxMatch = [];
            if length(Amps)>2
                for iAmp = 1:length(Amps)
                    ampIx = find(tSoundTargetAmplitude==Amps(iAmp));
                    Sb2IxAmp = intersect(ampIx,Sb2Ix);
                    Mb2IxAmp = intersect(ampIx,Mb2Ix);
                    nS(iAmp) = length(Sb2IxAmp);
                    nM(iAmp) = length(Mb2IxAmp);
                    if nS(iAmp)<nM(iAmp)
                        Sb2IxMatch = [Sb2IxMatch Sb2IxAmp];
                        Mb2IxMatch = [Mb2IxMatch Mb2IxAmp(randperm(nM(iAmp),nS(iAmp)))];
                    elseif nS(iAmp)>nM(iAmp)
                        Mb2IxMatch = [Mb2IxMatch Mb2IxAmp];
                        Sb2IxMatch = [Sb2IxMatch Sb2IxAmp(randperm(nS(iAmp),nM(iAmp)))];
                    else
                        Sb2IxMatch = [Sb2IxMatch Sb2IxAmp];
                        Mb2IxMatch = [Mb2IxMatch Mb2IxAmp];
                    end
                end
            else
                Sb2IxMatch = Sb2Ix;
                Mb2IxMatch = Mb2Ix;
            end
        else
            Sb2IxMatch = Sb2Ix;
            Mb2IxMatch = Mb2Ix;
        end
                
        
        %names of fields
        mouse(imouse).expt(s(:,imouse)).cells(1).name = 'all';
        mouse(imouse).expt(s(:,imouse)).cells(2).name = '0';
        mouse(imouse).expt(s(:,imouse)).cells(3).name = '45';
        mouse(imouse).expt(s(:,imouse)).cells(4).name = '90';
        mouse(imouse).expt(s(:,imouse)).cells(5).name = '135';
        mouse(imouse).expt(s(:,imouse)).cells(1).ind = 1:size(dataTC,2);
        mouse(imouse).expt(s(:,imouse)).cells(2).ind = cellsSelect{1};
        mouse(imouse).expt(s(:,imouse)).cells(3).ind = cellsSelect{2};
        mouse(imouse).expt(s(:,imouse)).cells(4).ind = cellsSelect{3};
        mouse(imouse).expt(s(:,imouse)).cells(5).ind = cellsSelect{4};
        mouse(imouse).expt(s(:,imouse)).align(1).name = 'press';
        mouse(imouse).expt(s(:,imouse)).align(2).name = 'release';
        mouse(imouse).expt(s(:,imouse)).align(3).name = 'target';
        mouse(imouse).expt(s(:,imouse)).align(4).name = 'prevStim';
        mouse(imouse).expt(s(:,imouse)).align(5).name = 'prevBaseStim';
        for i = 1:5
            mouse(imouse).expt(s(:,imouse)).align(i).av(1).name = 'visual';
            mouse(imouse).expt(s(:,imouse)).align(i).av(2).name = 'auditory';
            for ii = 1:2
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).name = 'hit';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).name = 'miss';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).name = 'FA';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).name = 'CR';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).name = 'hit_match';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).n = eval(['length(Sb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).n = eval(['length(Mb' num2str(ii) 'IxMatch)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).n = eval(['length(Fb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).n = eval(['length(Rb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).n = eval(['length(Sb' num2str(ii) 'IxMatch)']);                    
            end
        end
        
        
        
        %% Align data to lever up
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cLeverUp+post_event_frames-1 <  size(dataTC,1)),[],2)
            Data(:,:,itrial) = dataTC(cLeverUp(itrial)-pre_event_frames:cLeverUp(itrial)+post_event_frames-1,:);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
        end
        
        DataDFoverFavg = squeeze(mean(DataDFoverF(:,:,:),2));
        
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(trans_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(3).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(trans_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(5).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(trans_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(3).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(trans_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(5).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(mid_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(3).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(mid_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(5).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(mid_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(3).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(mid_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(5).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(late_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).late_resp = [mean(mean(DataDFoverF(late_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(3).late_resp = [mean(mean(DataDFoverF(late_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(late_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(5).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(late_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).late_resp = [mean(mean(DataDFoverF(late_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(3).late_resp = [mean(mean(DataDFoverF(late_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(late_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(5).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
        
        if doPlot
            figure;
            tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
            subplot(2,2,1)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
            subplot(2,2,2)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
            subplot(2,2,3)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
            subplot(2,2,4)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            alignYaxes
            title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
            suptitle([date_name ' ' mouse_name ' ' runstr '- align to lever release'])
            print([fnout 'release_align_SE_AV.pdf'], '-dpdf')
        end
        
        %% Align data to previous stim
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
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
            DataDF_CR(:,:,itrial) = bsxfun(@minus, Data_CR(:,:,itrial), mean(Data_CR(1:pre_event_frames,:,itrial),1));
            DataDFoverF_CR(:,:,itrial) = bsxfun(@rdivide, DataDF_CR(:,:,itrial), mean(Data_CR(1:pre_event_frames,:,itrial),1));
        end
        
        DataDFoverFavg = squeeze(mean(DataDFoverF(:,:,:),2));
        DataDFoverFavg_CR = squeeze(mean(DataDFoverF_CR(:,:,:),2));
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(1).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(trans_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(2).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(3).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(trans_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(4).trans_resp = [mean(mean(DataDFoverF_CR(trans_win,:,Rb1Ix),1),3) std(mean(DataDFoverF_CR(trans_win,:,Rb1Ix),1),[],3)./sqrt(length(Rb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(5).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(1).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(trans_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(2).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(3).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(trans_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(4).trans_resp = [mean(mean(DataDFoverF_CR(trans_win,:,Rb2Ix),1),3) std(mean(DataDFoverF_CR(trans_win,:,Rb2Ix),1),[],3)./sqrt(length(Rb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(5).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(1).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(mid_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(2).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(3).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(mid_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(4).mid_resp = [mean(mean(DataDFoverF_CR(mid_win,:,Rb1Ix),1),3) std(mean(DataDFoverF_CR(mid_win,:,Rb1Ix),1),[],3)./sqrt(length(Rb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(5).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(1).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(mid_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(2).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(3).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(mid_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(4).mid_resp = [mean(mean(DataDFoverF_CR(mid_win,:,Rb2Ix),1),3) std(mean(DataDFoverF_CR(mid_win,:,Rb2Ix),1),[],3)./sqrt(length(Rb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(5).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(1).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(late_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(2).late_resp = [mean(mean(DataDFoverF(late_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(3).late_resp = [mean(mean(DataDFoverF(late_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(late_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(4).late_resp = [mean(mean(DataDFoverF_CR(late_win,:,Rb1Ix),1),3) std(mean(DataDFoverF_CR(late_win,:,Rb1Ix),1),[],3)./sqrt(length(Rb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(5).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(1).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(late_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(2).late_resp = [mean(mean(DataDFoverF(late_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(3).late_resp = [mean(mean(DataDFoverF(late_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(late_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(4).late_resp = [mean(mean(DataDFoverF_CR(late_win,:,Rb2Ix),1),3) std(mean(DataDFoverF_CR(late_win,:,Rb2Ix),1),[],3)./sqrt(length(Rb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(5).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
        
        
        if doPlot 
            %Hit vs FA
            figure;
            tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
            subplot(2,2,1)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
            subplot(2,2,2)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
            subplot(2,2,3)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
            subplot(2,2,4)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            alignYaxes
            title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
            suptitle([date_name ' ' mouse_name ' ' runstr '- align to prev stim'])
            print([fnout 'prevStim_align_SE_AV.pdf'], '-dpdf')
           
            %Hit vs Miss
            figure;
            tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
            subplot(2,2,1)
            if length(Mb1IxMatch)>2
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), std(DataDFoverFavg(:,Mb1IxMatch),[],2)/sqrt(length(Mb1IxMatch)), 'r');
            else
                plot(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), 'r');
            end
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1IxMatch),2), std(DataDFoverFavg(:,Sb1IxMatch),[],2)/sqrt(length(Sb1IxMatch)), 'k');
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Visual trials: ' num2str(length(Sb1IxMatch)) ' Successes; ' num2str(length(Mb1IxMatch)) ' Misses']) 
            subplot(2,2,2)
            if length(Mb2Ix)>2
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb2Ix),2), std(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'r');
            else
                plot(tt, mean(DataDFoverFavg(:,Mb2Ix),2), 'r');
            end
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Mb2Ix)) ' Misses']) 
            subplot(2,2,3)
            if length(Mb1IxMatch)>2
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), std(DataDFoverFavg(:,Mb1IxMatch),[],2)/sqrt(length(Mb1IxMatch)), 'g');
            else
                plot(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), 'g');
            end
            hold on
            if length(Mb2Ix)>2
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb2Ix),2), std(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'k'); 
            else
                plot(tt, mean(DataDFoverFavg(:,Mb2Ix),2), 'r');
            end
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Missed trials: ' num2str(length(Mb1IxMatch)) ' Visual; ' num2str(length(Mb2Ix)) ' Auditory']) 
            subplot(2,2,4)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1IxMatch),2), std(DataDFoverFavg(:,Sb1IxMatch),[],2)/sqrt(length(Sb1IxMatch)), 'g');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            alignYaxes
            title(['Success trials: ' num2str(length(Sb1IxMatch)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
            suptitle([date_name ' ' mouse_name ' ' runstr '- align to previous stim'])
            print([fnout 'prevStim_align_SM_AV.pdf'], '-dpdf')
        end
        
        if doPlot
            %Hits vs FAs
            for iOri = 1:nOri
                DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
                figure;
                tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
                subplot(2,2,1)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
                subplot(2,2,2)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
                subplot(2,2,3)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
                subplot(2,2,4)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                alignYaxes
                title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
                suptitle([date_name ' ' mouse_name ' ' runstr '- align to previous stim on- ' num2str(Oris(iOri)) ' deg select: n = ' num2str(length(cellsSelect{iOri}))])
                print([fnout 'prevStimOn_align_SE_AV_' num2str(Oris(iOri)) 'pref.pdf'], '-dpdf')
            end
            for iOri = 1:nOri
                DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
                figure;
                tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
                subplot(2,2,1)
                if length(Mb1IxMatch)>2
                    shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), std(DataDFoverFavg(:,Mb1IxMatch),[],2)/sqrt(length(Mb1IxMatch)), 'r');
                else
                    plot(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), 'r');
                end
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1IxMatch),2), std(DataDFoverFavg(:,Sb1IxMatch),[],2)/sqrt(length(Sb1IxMatch)), 'k');
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Visual trials: ' num2str(length(Sb1IxMatch)) ' Successes; ' num2str(length(Mb1IxMatch)) ' Misses']) 
                subplot(2,2,2)
                if length(Mb2Ix)>2
                    shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb2Ix),2), std(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'r');
                else
                    plot(tt, mean(DataDFoverFavg(:,Mb2Ix),2), 'r');
                end
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Mb2Ix)) ' Misses']) 
                subplot(2,2,3)
                if length(Mb1IxMatch)>2
                    shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), std(DataDFoverFavg(:,Mb1IxMatch),[],2)/sqrt(length(Mb1IxMatch)), 'g');
                else
                    plot(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), 'g');
                end
                hold on
                if length(Mb2Ix)>2
                    shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb2Ix),2), std(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'k'); 
                else
                    plot(tt, mean(DataDFoverFavg(:,Mb2Ix),2), 'k');
                end
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Missed trials: ' num2str(length(Mb1IxMatch)) ' Visual; ' num2str(length(Mb2Ix)) ' Auditory']) 
                subplot(2,2,4)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1IxMatch),2), std(DataDFoverFavg(:,Sb1IxMatch),[],2)/sqrt(length(Sb1IxMatch)), 'g');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                alignYaxes
                title(['Success trials: ' num2str(length(Sb1IxMatch)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
                suptitle([date_name ' ' mouse_name ' ' runstr '- align to previous stim on- ' num2str(Oris(iOri)) ' deg select: n = ' num2str(length(cellsSelect{iOri}))])
                print([fnout 'prevStimOn_align_SM_AV_' num2str(Oris(iOri)) 'pref.pdf'], '-dpdf')
            end
        end
        
        %% Align data to previous baseline stim
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        Data_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cStimOn+post_event_frames-1 <  size(dataTC,1)),[],2)
            Data(:,:,itrial) = dataTC(cStimOn(itrial)-pre_event_frames:cStimOn(itrial)+post_event_frames-1,:);
            Data_CR(:,:,itrial) = dataTC(cStimOn(itrial)-pre_event_frames-cycTime:cStimOn(itrial)+post_event_frames-cycTime-1,:);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
            DataDF_CR(:,:,itrial) = bsxfun(@minus, Data_CR(:,:,itrial), mean(Data_CR(1:pre_event_frames,:,itrial),1));
            DataDFoverF_CR(:,:,itrial) = bsxfun(@rdivide, DataDF_CR(:,:,itrial), mean(Data_CR(1:pre_event_frames,:,itrial),1));
        end
        
        DataDFoverFavg = squeeze(mean(DataDFoverF(:,:,:),2));
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(1).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(trans_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(2).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(3).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(trans_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(4).trans_resp = [mean(mean(DataDFoverF_CR(trans_win,:,Rb1Ix),1),3) std(mean(DataDFoverF_CR(trans_win,:,Rb1Ix),1),[],3)./sqrt(length(Rb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(5).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(1).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(trans_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(2).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(3).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(trans_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(4).trans_resp = [mean(mean(DataDFoverF_CR(trans_win,:,Rb2Ix),1),3) std(mean(DataDFoverF_CR(trans_win,:,Rb2Ix),1),[],3)./sqrt(length(Rb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(5).trans_resp = [mean(mean(DataDFoverF(trans_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(trans_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(1).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(mid_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(2).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(3).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(mid_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(4).mid_resp = [mean(mean(DataDFoverF_CR(mid_win,:,Rb1Ix),1),3) std(mean(DataDFoverF_CR(mid_win,:,Rb1Ix),1),[],3)./sqrt(length(Rb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(5).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(1).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(mid_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(2).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(3).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(mid_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(4).mid_resp = [mean(mean(DataDFoverF_CR(mid_win,:,Rb2Ix),1),3) std(mean(DataDFoverF_CR(mid_win,:,Rb2Ix),1),[],3)./sqrt(length(Rb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(5).mid_resp = [mean(mean(DataDFoverF(mid_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(mid_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(1).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb1Ix),1),3) std(mean(DataDFoverF(late_win,:,Sb1Ix),1),[],3)./sqrt(length(Sb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(2).late_resp = [mean(mean(DataDFoverF(late_win,:,Mb1IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Mb1IxMatch),1),[],3)./sqrt(length(Mb1IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(3).late_resp = [mean(mean(DataDFoverF(late_win,:,Fb1Ix),1),3) std(mean(DataDFoverF(late_win,:,Fb1Ix),1),[],3)./sqrt(length(Fb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(4).late_resp = [mean(mean(DataDFoverF_CR(late_win,:,Rb1Ix),1),3) std(mean(DataDFoverF_CR(late_win,:,Rb1Ix),1),[],3)./sqrt(length(Rb1Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(1).outcome(5).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb1IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Sb1IxMatch),1),[],3)./sqrt(length(Sb1IxMatch))]; 
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(1).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb2Ix),1),3) std(mean(DataDFoverF(late_win,:,Sb2Ix),1),[],3)./sqrt(length(Sb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(2).late_resp = [mean(mean(DataDFoverF(late_win,:,Mb2IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Mb2IxMatch),1),[],3)./sqrt(length(Mb2IxMatch))];        
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(3).late_resp = [mean(mean(DataDFoverF(late_win,:,Fb2Ix),1),3) std(mean(DataDFoverF(late_win,:,Fb2Ix),1),[],3)./sqrt(length(Fb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(4).late_resp = [mean(mean(DataDFoverF_CR(late_win,:,Rb2Ix),1),3) std(mean(DataDFoverF_CR(late_win,:,Rb2Ix),1),[],3)./sqrt(length(Rb2Ix))];
        mouse(imouse).expt(s(:,imouse)).align(4).av(2).outcome(5).late_resp = [mean(mean(DataDFoverF(late_win,:,Sb2IxMatch),1),3) std(mean(DataDFoverF(late_win,:,Sb2IxMatch),1),[],3)./sqrt(length(Sb2IxMatch))]; 
       
            %Hit vs FA
        if doPlot    
            figure;
            tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
            subplot(2,2,1)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
            subplot(2,2,2)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
            subplot(2,2,3)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
            subplot(2,2,4)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            alignYaxes
            title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
            suptitle([date_name ' ' mouse_name ' ' runstr '- align to prev base stim'])
            print([fnout 'prevStim_align_SE_AV.pdf'], '-dpdf')
            
            %Hit vs Miss
            figure;
            tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
            subplot(2,2,1)
            if length(Mb1IxMatch)>2
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), std(DataDFoverFavg(:,Mb1IxMatch),[],2)/sqrt(length(Mb1IxMatch)), 'r');
            else
                plot(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), 'r');
            end
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1IxMatch),2), std(DataDFoverFavg(:,Sb1IxMatch),[],2)/sqrt(length(Sb1IxMatch)), 'k');
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Visual trials: ' num2str(length(Sb1IxMatch)) ' Successes; ' num2str(length(Mb1IxMatch)) ' Misses']) 
            subplot(2,2,2)
            if length(Mb2Ix)>2
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb2Ix),2), std(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'r');
            else
                plot(tt, mean(DataDFoverFavg(:,Mb2Ix),2), 'r');
            end
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Mb2Ix)) ' Misses']) 
            subplot(2,2,3)
            if length(Mb1IxMatch)>2
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), std(DataDFoverFavg(:,Mb1IxMatch),[],2)/sqrt(length(Mb1IxMatch)), 'g');
            else
                plot(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), 'g');
            end
            hold on
            if length(Mb2Ix)>2
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb2Ix),2), std(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'k'); 
            else
                plot(tt, mean(DataDFoverFavg(:,Mb2Ix),2), 'r');
            end
            xlim([-500 1500])
            hold on
            vline(0,'k')
            title(['Missed trials: ' num2str(length(Mb1IxMatch)) ' Visual; ' num2str(length(Mb2Ix)) ' Auditory']) 
            subplot(2,2,4)
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1IxMatch),2), std(DataDFoverFavg(:,Sb1IxMatch),[],2)/sqrt(length(Sb1IxMatch)), 'g');
            hold on
            shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
            xlim([-500 1500])
            hold on
            vline(0,'k')
            alignYaxes
            title(['Success trials: ' num2str(length(Sb1IxMatch)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
            suptitle([date_name ' ' mouse_name ' ' runstr '- align to previous base stim'])
            print([fnout 'prevStim_align_SM_AV.pdf'], '-dpdf')
        end
        
        if doPlot
            %Hits vs FAs
            for iOri = 1:nOri
                DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
                figure;
                tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
                subplot(2,2,1)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
                subplot(2,2,2)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
                subplot(2,2,3)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb1Ix),2), std(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Fb2Ix),2), std(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
                subplot(2,2,4)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1Ix),2), std(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                alignYaxes
                title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
                suptitle([date_name ' ' mouse_name ' ' runstr '- align to previous base stim on- ' num2str(Oris(iOri)) ' deg select: n = ' num2str(length(cellsSelect{iOri}))])
                print([fnout 'prevBaseStimOn_align_SE_AV_' num2str(Oris(iOri)) 'pref.pdf'], '-dpdf')
            end
            for iOri = 1:nOri
                DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
                figure;
                tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
                subplot(2,2,1)
                if length(Mb1IxMatch)>2
                    shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), std(DataDFoverFavg(:,Mb1IxMatch),[],2)/sqrt(length(Mb1IxMatch)), 'r');
                else
                    plot(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), 'r');
                end
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1IxMatch),2), std(DataDFoverFavg(:,Sb1IxMatch),[],2)/sqrt(length(Sb1IxMatch)), 'k');
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Visual trials: ' num2str(length(Sb1IxMatch)) ' Successes; ' num2str(length(Mb1IxMatch)) ' Misses']) 
                subplot(2,2,2)
                if length(Mb2Ix)>2
                    shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb2Ix),2), std(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'r');
                else
                    plot(tt, mean(DataDFoverFavg(:,Mb2Ix),2), 'r');
                end
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Mb2Ix)) ' Misses']) 
                subplot(2,2,3)
                if length(Mb1IxMatch)>2
                    shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), std(DataDFoverFavg(:,Mb1IxMatch),[],2)/sqrt(length(Mb1IxMatch)), 'g');
                else
                    plot(tt, mean(DataDFoverFavg(:,Mb1IxMatch),2), 'g');
                end
                hold on
                if length(Mb2Ix)>2
                    shadedErrorBar(tt, mean(DataDFoverFavg(:,Mb2Ix),2), std(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'k'); 
                else
                    plot(tt, mean(DataDFoverFavg(:,Mb2Ix),2), 'k');
                end
                xlim([-500 1500])
                hold on
                vline(0,'k')
                title(['Missed trials: ' num2str(length(Mb1IxMatch)) ' Visual; ' num2str(length(Mb2Ix)) ' Auditory']) 
                subplot(2,2,4)
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb1IxMatch),2), std(DataDFoverFavg(:,Sb1IxMatch),[],2)/sqrt(length(Sb1IxMatch)), 'g');
                hold on
                shadedErrorBar(tt, mean(DataDFoverFavg(:,Sb2Ix),2), std(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                xlim([-500 1500])
                hold on
                vline(0,'k')
                alignYaxes
                title(['Success trials: ' num2str(length(Sb1IxMatch)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
                suptitle([date_name ' ' mouse_name ' ' runstr '- align to previous base stim on- ' num2str(Oris(iOri)) ' deg select: n = ' num2str(length(cellsSelect{iOri}))])
                print([fnout 'prevBaseStimOn_align_SM_AV_' num2str(Oris(iOri)) 'pref.pdf'], '-dpdf')
            end
        end
    end
    mouse_str = [];
    for imouse = 1:size(av,2)
        mouse_str = [mouse_str 'i' num2str(av(imouse).mouse) '_'];  
    end
   save(fullfile(rc.caOutputDir, [date '_' mouse_str 'CaSummary.mat']), 'mouse');
end

        
        
        

