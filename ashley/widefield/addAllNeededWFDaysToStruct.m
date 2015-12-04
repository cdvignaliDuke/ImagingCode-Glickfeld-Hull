%% run one at a time
function addAllNeededWFDaysToStruct(mouse)
while true
    close all;
    rc = behavConstsWF(mouse);
    event_buffer_time = 2000;
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
            load(fullfile(rc.structOutput, [mouse '_' dateStr '_roi_TC.mat']), 'roi_TC')
        else
            img_data = [];   
            offset = 0;
            for irun = 1:nrun
                temp_data = readtiff(fullfile(data_folder,[dateStr '_' mouse],[mouse behav_run num2str(runs(irun))], [mouse behav_run num2str(runs(irun)) '_MMStack.ome.tif']));
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
            load(fullfile(rc.structOutput, [mouse '_' dateStr '_reg_data.mat']), 'img_avg', 'roi_avg', 'input_points', 'base_points')
            sz_target  = size(roi_avg);
            mytform    = maketform('affine',input_points(1:3,:), base_points(1:3,:));
            
            numWorkers = 10;
            pool = parpool(numWorkers);
            sz_orig = size(img_data);
            registered = zeros(sz_orig(1),sz_orig(2), sz_orig(3)./numWorkers,numWorkers);
            tic
            parfor i = 1:numWorkers
                registered(:,:,:,i) = imtransform(img_data(:,:,1+((sz_orig(3)./numWorkers)*(i-1)):(sz_orig(3)./numWorkers)*i),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);   
            end
            toc
            delete(gcp)
            registered = reshape(registered, sz_orig);
            clear img_data;

            %extract timecourses
            load(fullfile(rc.structOutput, [mouse '_' roi_date '_roi_masks.mat']))
            roi_TC = stackGetTimeCourses(registered, mask_cell_V);
            save(fullfile(rc.structOutput, [mouse '_' dateStr '_roi_TC.mat']), 'roi_TC')
            clear registered
        end
        
        %align ROIs to events
        nROI = size(roi_TC,2);
        ntrials = length(cLeverDown);
        frame_rate = str2num(xd.FrameRate{indexRowN});
        event_buffer_frames = ceil(event_buffer_time*(frame_rate./1000));
        pressAlign = zeros(2*event_buffer_frames,nROI,ntrials);
        targetAlign = zeros(2*event_buffer_frames,nROI,ntrials);
        releaseAlign = zeros(2*event_buffer_frames,nROI,ntrials);
        for itrial = 1:ntrials
            if cLeverDown(itrial)-event_buffer_frames>0 & cLeverDown(itrial)+event_buffer_frames<size(roi_TC,1)
                pressAlign(:,:,itrial) = roi_TC(cLeverDown(itrial)-event_buffer_frames:cLeverDown(itrial)+event_buffer_frames-1,:);
            else
                pressAlign(:,:,itrial) = nan(2*event_buffer_frames,nROI);
            end
            if ~isnan(cTargetOn(itrial)) & cTargetOn(itrial)+event_buffer_frames<size(roi_TC,1)
                targetAlign(:,:,itrial) = roi_TC(cTargetOn(itrial)-event_buffer_frames:cTargetOn(itrial)+event_buffer_frames-1,:);
            else
                targetAlign(:,:,itrial) = nan(2*event_buffer_frames,nROI);
            end
            if cLeverUp(itrial)+event_buffer_frames<size(roi_TC,1)
                releaseAlign(:,:,itrial) = roi_TC(cLeverUp(itrial)-event_buffer_frames:cLeverUp(itrial)+event_buffer_frames-1,:);
            else
                releaseAlign(:,:,itrial) = nan(2*event_buffer_frames,nROI);
            end
        end
        
        dS(dayN).date = dateStr;
        dS(dayN).pressAlign = pressAlign;
        dS(dayN).targetAlign = targetAlign;
        dS(dayN).releaseAlign = releaseAlign;
        dS(dayN).trialOutcomeCell = input.trialOutcomeCell;
        dS(dayN).holdTimesMs = cell2mat_padded(input.holdTimesMs);
        dS(dayN).reactTimesMs = cell2mat_padded(input.reactTimesMs);
        dS(dayN).reactTimesMs = cell2mat_padded(input.tTotalReqHoldTimeMs);
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

