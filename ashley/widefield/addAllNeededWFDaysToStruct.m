%% run one at a time
function addAllNeededWFDaysToStruct(mouse)
while true
    close all;
    rc = behavConstsWF(mouse);
    event_buffer_time = 2000;
    f_win = [event_buffer_time/3:2*(event_buffer_time/3)];
    if ~exist(fullfile(rc.structOutput, [mouse '_dSStruct.mat']))
        dS = struct;
        dayN = 1;
    else
        load(fullfile(rc.structOutput, [mouse '_dSStruct.mat']))
        dayN = length(dS)+1;
    end
    [xd indexRowN] = findNextIndexForDataStructWF(rc, dS);

    if isempty(indexRowN)
        disp('No more records to do');
        dateStr = [];
        return
    else

        % do the one we found
        dateStr = xd.DateStr{indexRowN};
        fprintf('Starting subject %s, date %s\n', ...
        mouse, dateStr);

        % use imaging data to determine number of runs
        eval(['i' mouse '_paths'])
        if mouse == '633'
            data_folder = fullfile(data_folder, [dateStr '_' mouse]);
        else
            data_folder = fullfile(data_folder, dateStr);
        end
        nrun = str2num(xd.Runs{indexRowN});
        runs = [];
        for irun = 1:nrun
            runs = [runs eval(['xd.ImgDataRun' num2str(irun) '(indexRowN)'])];
        end

        % load behavioral data for each run 
        for irun = 1:nrun
            timeStr = num2str(eval(['xd.MatFileRun' num2str(irun) '(indexRowN)']));
            fn_mworks = fullfile(behav_folder,['data-i' mouse '-' dateStr '-' timeStr '.mat']);
            if irun == 1
                input = mwLoadData(fn_mworks, [], []);
            else
                input = [input mwLoadData(fn_mworks, [], [])];
            end
        end
        input = concatenateDataBlocks(input);

        %extract counter info
        run_trials = input.trialsSinceReset;
        cLeverDown = cell2mat(input.cLeverDown);
        cLeverUp = cell2mat(input.cLeverUp);
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cStimOn = celleqel2mat_padded(input.cStimOn);
        cItiStart = cell2mat(input.cItiStart);

        %if timecourses don't exist load imaging data and account for counter offset
        if exist(fullfile(rc.structOutput, [mouse '_' dateStr '_roi_TC.mat']), 'file')
            load(fullfile(rc.structOutput, [mouse '_' dateStr '_roi_TC.mat']))
        else
            img_data = [];   
            offset = 0;
            for irun = 1:nrun
                temp_data = readtiff(fullfile(data_folder,[behav_run num2str(runs(irun))], [behav_run num2str(runs(irun)) suffix '.tif']));
                img_data = cat(3,img_data,temp_data); 
                offset = offset+size(temp_data,3);
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
            clear temp_data

            %register to retinotopy dataset
            load(fullfile(rc.structOutput, [mouse '_' dateStr '_reg_data.mat']))
            sz_target  = size(roi_avg);
            mytform    = maketform('affine',input_points(1:3,:), base_points(1:3,:));
            
            numWorkers = 10;
            if isempty(gcp)
                pool = parpool(numWorkers);
            end
            sz_orig = size(img_data);
            sz_orig(3) = (floor(sz_orig(3)./numWorkers)).*numWorkers;
            registered = zeros(sz_orig(1),sz_orig(2), sz_orig(3)./numWorkers,numWorkers);
            tic
            parfor i = 1:numWorkers
                registered(:,:,:,i) = imtransform(img_data(:,:,1+((sz_orig(3)./numWorkers)*(i-1)):(sz_orig(3)./numWorkers)*i),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);   
            end
            toc
            delete(gcp)
            registered = reshape(registered, sz_orig);
            
            %sort by alignment
            
            ntrials = length(cLeverDown);
            frame_rate = str2num(xd.FrameRate{indexRowN});
            event_buffer_frames = ceil(event_buffer_time*(frame_rate./1000));
            f_win_frames = ceil(f_win.*(frame_rate./1000));
            sz = size(registered);
            pressAlign = zeros(sz(1), sz(2), 2*event_buffer_frames,ntrials);
            targetAlign = zeros(sz(1), sz(2), 2*event_buffer_frames,ntrials);
            releaseAlign = zeros(sz(1), sz(2), 2*event_buffer_frames,ntrials);
            for itrial = 1:ntrials
                if cLeverDown(itrial)-event_buffer_frames>0 & cLeverDown(itrial)+event_buffer_frames<sz(3)
                    pressAlign(:,:,:,itrial) = registered(:,:,cLeverDown(itrial)-event_buffer_frames:cLeverDown(itrial)+event_buffer_frames-1);
                else
                    pressAlign(:,:,:,itrial) = nan(sz(1), sz(2), 2*event_buffer_frames);
                end
                if ~isnan(cTargetOn(itrial)) & cTargetOn(itrial)+event_buffer_frames<sz(3)
                    targetAlign(:,:,:,itrial) = registered(:,:,cTargetOn(itrial)-event_buffer_frames:cTargetOn(itrial)+event_buffer_frames-1);
                else
                    targetAlign(:,:,:,itrial) = nan(sz(1), sz(2), 2*event_buffer_frames);
                end
                if cLeverUp(itrial)+event_buffer_frames<sz(3)
                    releaseAlign(:,:,:,itrial) = registered(:,:,cLeverUp(itrial)-event_buffer_frames:cLeverUp(itrial)+event_buffer_frames-1);
                else
                    releaseAlign(:,:,:,itrial) = nan(sz(1), sz(2), 2*event_buffer_frames);
                end
            end
            
            clear img_data
            
            %transform to dF/F
            F = mean(pressAlign(:,:,f_win_frames,:),3);
            pressAlign_dF = bsxfun(@minus,pressAlign,F);
            pressAlign_dFoverF = bsxfun(@rdivide,pressAlign_dF,F);
            targetAlign_dF = bsxfun(@minus,targetAlign,F);
            targetAlign_dFoverF = bsxfun(@rdivide,targetAlign_dF,F);
            releaseAlign_dF = bsxfun(@minus,releaseAlign,F);
            releaseAlign_dFoverF = bsxfun(@rdivide,releaseAlign_dF,F);
            clear pressAlign releaseAlign targetAlign pressAlign_dF releaseAlign_dF targetAlign_dF
            
            %sort by trial type
            successIx = find(strcmp(input.trialOutcomeCell,'success'));
            failureIx = find(strcmp(input.trialOutcomeCell,'failure'));
            missedIx = find(strcmp(input.trialOutcomeCell,'ignore'));
            
            if ~exist(fullfile(rc.structOutput, dateStr))
                mkdir(fullfile(rc.structOutput, dateStr));
            end
            save(fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_outcomeIx.mat']),'successIx', 'failureIx', 'missedIx');
            
            pressAlign_dFoverF_SIx = squeeze(mean(pressAlign_dFoverF(:,:,:,successIx),4));
            pressAlign_dFoverF_FIx = squeeze(mean(pressAlign_dFoverF(:,:,:,failureIx),4));
            pressAlign_dFoverF_MIx = squeeze(mean(pressAlign_dFoverF(:,:,:,missedIx),4));
            targetAlign_dFoverF_SIx = squeeze(mean(targetAlign_dFoverF(:,:,:,successIx),4));
            targetAlign_dFoverF_MIx = squeeze(mean(targetAlign_dFoverF(:,:,:,missedIx),4));
            releaseAlign_dFoverF_SIx = squeeze(mean(releaseAlign_dFoverF(:,:,:,successIx),4));
            releaseAlign_dFoverF_FIx = squeeze(mean(releaseAlign_dFoverF(:,:,:,failureIx),4));
            clear pressAlign_dFoverF releaseAlign_dFoverF targetAlign_dFoverF
            
            writetiff(pressAlign_dFoverF_SIx, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignSIx.tif']));
            writetiff(pressAlign_dFoverF_MIx, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignMIx.tif']));
            writetiff(pressAlign_dFoverF_FIx, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignFIx.tif'])); 
            writetiff(targetAlign_dFoverF_SIx, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_targetAlignSIx.tif']));   
            writetiff(targetAlign_dFoverF_MIx, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_targetAlignMIx.tif'])); 
            writetiff(releaseAlign_dFoverF_SIx, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_releaseAlignSIx.tif']));    
            writetiff(releaseAlign_dFoverF_FIx, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_releaseAlignFIx.tif'])); 
            
            %extract timecourses
            load(fullfile(rc.structOutput, [mouse '_' roi_date '_roi_masks.mat']))
            reg_avg = mean(registered,3);
            mask_cell_V(find(reg_avg==0)) = 0;
            roiTC = stackGetTimeCourses(registered, mask_cell_V);
            roiTC_np = zeros(size(roiTC));
            nROI = size(roiTC,2);
            for i = 1:nROI
                np = np_V(:,:,i);
                np(find(reg_avg==0)) = 0;
                roiTC_np(:,i) = stackGetTimeCourses(registered, np);
            end
            save(fullfile(rc.structOutput, [mouse '_' dateStr '_roi_TC.mat']), 'roiTC', 'roiTC_np')
            clear registered
        end
        
        %align ROIs to events
        nROI = size(roiTC,2);
        ntrials = length(cLeverDown);
        frame_rate = str2num(xd.FrameRate{indexRowN});
        event_buffer_frames = ceil(event_buffer_time*(frame_rate./1000));
        pressAlign = zeros(2*event_buffer_frames,nROI,2,ntrials);
        targetAlign = zeros(2*event_buffer_frames,nROI,2,ntrials);
        releaseAlign = zeros(2*event_buffer_frames,nROI,2,ntrials);
        for itrial = 1:ntrials
            if cLeverDown(itrial)-event_buffer_frames>0 & cLeverDown(itrial)+event_buffer_frames<size(roiTC,1)
                pressAlign(:,:,1,itrial) = roiTC(cLeverDown(itrial)-event_buffer_frames:cLeverDown(itrial)+event_buffer_frames-1,:);
                pressAlign(:,:,2,itrial) = roiTC_np(cLeverDown(itrial)-event_buffer_frames:cLeverDown(itrial)+event_buffer_frames-1,:);
            else
                pressAlign(:,:,:,itrial) = nan(2*event_buffer_frames,nROI,2);
            end
            if ~isnan(cTargetOn(itrial)) & cTargetOn(itrial)+event_buffer_frames<size(roiTC,1)
                targetAlign(:,:,1,itrial) = roiTC(cTargetOn(itrial)-event_buffer_frames:cTargetOn(itrial)+event_buffer_frames-1,:);
                targetAlign(:,:,2,itrial) = roiTC_np(cTargetOn(itrial)-event_buffer_frames:cTargetOn(itrial)+event_buffer_frames-1,:);
            else
                targetAlign(:,:,:,itrial) = nan(2*event_buffer_frames,nROI,2);
            end
            if cLeverUp(itrial)+event_buffer_frames<size(roiTC,1)
                releaseAlign(:,:,1,itrial) = roiTC(cLeverUp(itrial)-event_buffer_frames:cLeverUp(itrial)+event_buffer_frames-1,:);
                releaseAlign(:,:,2,itrial) = roiTC_np(cLeverUp(itrial)-event_buffer_frames:cLeverUp(itrial)+event_buffer_frames-1,:);
            else
                releaseAlign(:,:,:,itrial) = nan(2*event_buffer_frames,nROI,2);
            end
        end
        
        dS(dayN).date = dateStr;
        dS(dayN).pressAlign = pressAlign;
        dS(dayN).targetAlign = targetAlign;
        dS(dayN).releaseAlign = releaseAlign;
        dS(dayN).trialOutcomeCell = input.trialOutcomeCell;
        dS(dayN).holdTimesMs = cell2mat_padded(input.holdTimesMs);
        dS(dayN).reactTimesMs = cell2mat_padded(input.reactTimesMs);
        dS(dayN).reqHoldTimesMs = cell2mat_padded(input.tTotalReqHoldTimeMs);
        dS(dayN).reactTimeMs = input.reactTimeMs;
        dS(dayN).stimOnTimeMs = input.stimOnTimeMs;
        dS(dayN).eventBufferFrames = event_buffer_frames;
        dS(dayN).frameRate = frame_rate;
        
        save(fullfile(rc.structOutput, [mouse '_dSStruct.mat']), 'dS'); 
    end 
            
    if isempty(dateStr)
        break; 
    end
end

