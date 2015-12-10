function makeAlignedMovies(mouse)
while true
    rc = behavConstsWF(mouse);
    event_buffer_time = 2000;
    f_win = [event_buffer_time/3:2*(event_buffer_time/3)];
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols); 
    eval(['i' mouse '_paths'])
    load(fullfile(rc.structOutput, [mouse '_' roi_date '_roi_masks.mat']))
    sz = size(mask_cell);
    pre_press_SIx_cat = zeros(sz(1), sz(2), xd.nRows);
    post_press_SIx_cat= zeros(sz(1), sz(2), xd.nRows);
    pre_target_SIx_cat= zeros(sz(1), sz(2), xd.nRows);
    post_target_SIx_cat= zeros(sz(1), sz(2), xd.nRows);
    post_release_SIx_cat= zeros(sz(1), sz(2), xd.nRows);
    pre_press_MIx_cat= zeros(sz(1), sz(2), xd.nRows);
    post_press_MIx_cat= zeros(sz(1), sz(2), xd.nRows);
    pre_target_MIx_cat= zeros(sz(1), sz(2), xd.nRows);
    post_target_MIx_cat= zeros(sz(1), sz(2), xd.nRows);
    pre_press_FIx_cat= zeros(sz(1), sz(2), xd.nRows);
    post_press_FIx_cat= zeros(sz(1), sz(2), xd.nRows);
    post_release_FIx_cat= zeros(sz(1), sz(2), xd.nRows);
    for iD = 1:xd.nRows
        dateStr = xd.DateStr{iD};
        if ~exist(fullfile(rc.structOutput,dateStr, [mouse '_' dateStr '_pressAlignSIx.tif']))
            if ~exist(fullfile(rc.structOutput,dateStr))
                mkdir(fullfile(rc.structOutput,dateStr))
            end
            
            fprintf('Starting subject %s, date %s\n', ...
            mouse, dateStr);
            indexRowN = iD;
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
            
            %load imaging data
            img_data = [];   
            offset = 0;
            for irun = 1:nrun
                temp_data = readtiff(fullfile(data_folder,[dateStr '_' mouse_name],[behav_run num2str(runs(irun))], [behav_run num2str(runs(irun)) suffix '.tif']));
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
            
            %sort by alignment
            
            ntrials = length(cLeverDown);
            frame_rate = str2num(xd.FrameRate{indexRowN});
            event_buffer_frames = ceil(event_buffer_time*(frame_rate./1000));
            f_win_frames = ceil(f_win.*(frame_rate./1000));
            sz = size(img_data);
            pressAlign = zeros(sz(1), sz(2), 2*event_buffer_frames,ntrials);
            targetAlign = zeros(sz(1), sz(2), 2*event_buffer_frames,ntrials);
            releaseAlign = zeros(sz(1), sz(2), 2*event_buffer_frames,ntrials);
            for itrial = 1:ntrials
                if cLeverDown(itrial)-event_buffer_frames>0 & cLeverDown(itrial)+event_buffer_frames<sz(3)
                    pressAlign(:,:,:,itrial) = img_data(:,:,cLeverDown(itrial)-event_buffer_frames:cLeverDown(itrial)+event_buffer_frames-1);
                else
                    pressAlign(:,:,:,itrial) = nan(sz(1), sz(2), 2*event_buffer_frames);
                end
                if ~isnan(cTargetOn(itrial)) & cTargetOn(itrial)+event_buffer_frames<sz(3)
                    targetAlign(:,:,:,itrial) = img_data(:,:,cTargetOn(itrial)-event_buffer_frames:cTargetOn(itrial)+event_buffer_frames-1);
                else
                    targetAlign(:,:,:,itrial) = nan(sz(1), sz(2), 2*event_buffer_frames);
                end
                if cLeverUp(itrial)+event_buffer_frames<sz(3)
                    releaseAlign(:,:,:,itrial) = img_data(:,:,cLeverUp(itrial)-event_buffer_frames:cLeverUp(itrial)+event_buffer_frames-1);
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
            
            save(fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_outcomeIx.mat']),'successIx', 'failureIx', 'missedIx');
            
            pressAlign_dFoverF_SIx = squeeze(mean(pressAlign_dFoverF(:,:,:,successIx),4));
            pressAlign_dFoverF_FIx = squeeze(mean(pressAlign_dFoverF(:,:,:,failureIx),4));
            pressAlign_dFoverF_MIx = squeeze(mean(pressAlign_dFoverF(:,:,:,missedIx),4));
            targetAlign_dFoverF_SIx = squeeze(mean(targetAlign_dFoverF(:,:,:,successIx),4));
            targetAlign_dFoverF_MIx = squeeze(mean(targetAlign_dFoverF(:,:,:,missedIx),4));
            releaseAlign_dFoverF_SIx = squeeze(mean(releaseAlign_dFoverF(:,:,:,successIx),4));
            releaseAlign_dFoverF_FIx = squeeze(mean(releaseAlign_dFoverF(:,:,:,failureIx),4));
            clear pressAlign_dFoverF releaseAlign_dFoverF targetAlign_dFoverF
            
            %align to roi
            load(fullfile(rc.structOutput, [mouse '_' dateStr '_reg_data.mat']))
            sz_target  = size(roi_avg);
            mytform    = maketform('affine',input_points(1:3,:), base_points(1:3,:));
            sz_orig = size(pressAlign_dFoverF_SIx);
            
            pressAlign_dFoverF_SIx_reg = imtransform(pressAlign_dFoverF_SIx,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);   
            pressAlign_dFoverF_MIx_reg = imtransform(pressAlign_dFoverF_MIx,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]); 
            pressAlign_dFoverF_FIx_reg = imtransform(pressAlign_dFoverF_FIx,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]); 
            targetAlign_dFoverF_SIx_reg = imtransform(targetAlign_dFoverF_SIx,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);   
            targetAlign_dFoverF_MIx_reg = imtransform(targetAlign_dFoverF_MIx,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]); 
            releaseAlign_dFoverF_SIx_reg = imtransform(releaseAlign_dFoverF_SIx,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);   
            releaseAlign_dFoverF_FIx_reg = imtransform(releaseAlign_dFoverF_FIx,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
            clear pressAlign_dFoverF_SIx pressAlign_dFoverF_MIx pressAlign_dFoverF_FIx releaseAlign_dFoverF_SIx releaseAlign_dFoverF_FIx targetAlign_dFoverF_SIx targetAlign_dFoverF_MIx
            
            writetiff(pressAlign_dFoverF_SIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignSIx.tif']));
            writetiff(pressAlign_dFoverF_MIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignMIx.tif']));
            writetiff(pressAlign_dFoverF_FIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignFIx.tif'])); 
            writetiff(targetAlign_dFoverF_SIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_targetAlignSIx.tif']));   
            writetiff(targetAlign_dFoverF_MIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_targetAlignMIx.tif'])); 
            writetiff(releaseAlign_dFoverF_SIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_releaseAlignSIx.tif']));    
            writetiff(releaseAlign_dFoverF_FIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_releaseAlignFIx.tif'])); 
        else
            load(fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_outcomeIx.mat']),'successIx', 'failureIx', 'missedIx');
            pressAlign_dFoverF_SIx_reg = readtiff(pressAlign_dFoverF_SIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignSIx.tif']));
            pressAlign_dFoverF_MIx_reg = readtiff(pressAlign_dFoverF_MIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignMIx.tif']));
            pressAlign_dFoverF_FIx_reg = readtiff(pressAlign_dFoverF_FIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_pressAlignFIx.tif'])); 
            targetAlign_dFoverF_SIx_reg = readtiff(targetAlign_dFoverF_SIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_targetAlignSIx.tif']));   
            targetAlign_dFoverF_MIx_reg = readtiff(targetAlign_dFoverF_MIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_targetAlignMIx.tif'])); 
            releaseAlign_dFoverF_SIx_reg = readtiff(releaseAlign_dFoverF_SIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_releaseAlignSIx.tif']));    
            releaseAlign_dFoverF_FIx_reg = readtiff(releaseAlign_dFoverF_FIx_reg, fullfile(rc.structOutput, dateStr,[mouse '_' dateStr '_releaseAlignFIx.tif'])); 
        end
        pre_press_SIx_cat(:,:,iD) = mean(pressAlign_dFoverF_SIx_reg(:,:,51:60),3);
        post_press_SIx_cat(:,:,iD) = mean(pressAlign_dFoverF_SIx_reg(:,:,61:70),3);
        pre_target_SIx_cat(:,:,iD) = mean(targetAlign_dFoverF_SIx_reg(:,:,51:60),3);
        post_target_SIx_cat(:,:,iD) = mean(targetAlign_dFoverF_SIx_reg(:,:,65:74),3);
        post_release_SIx_cat(:,:,iD) = mean(releaseAlign_dFoverF_SIx_reg(:,:,65:74),3);
        if length(missedIx>0)
            pre_press_MIx_cat(:,:,iD) = mean(pressAlign_dFoverF_MIx_reg(:,:,51:60),3);
            post_press_MIx_cat(:,:,iD) = mean(pressAlign_dFoverF_MIx_reg(:,:,61:70),3);
            pre_target_MIx_cat(:,:,iD) = mean(targetAlign_dFoverF_MIx_reg(:,:,51:60),3);
            post_target_MIx_cat(:,:,iD) = mean(targetAlign_dFoverF_MIx_reg(:,:,65:74),3);
        end
        pre_press_FIx_cat(:,:,iD) = mean(pressAlign_dFoverF_FIx_reg(:,:,51:60),3);
        post_press_FIx_cat(:,:,iD) = mean(pressAlign_dFoverF_FIx_reg(:,:,61:70),3);
        post_release_FIx_cat(:,:,iD) = mean(releaseAlign_dFoverF_FIx_reg(:,:,65:74),3);
        clear pressAlign_dFoverF_SIx_reg pressAlign_dFoverF_MIx_reg pressAlign_dFoverF_FIx_reg releaseAlign_dFoverF_SIx_reg releaseAlign_dFoverF_FIx_reg targetAlign_dFoverF_SIx_reg targetAlign_dFoverF_MIx_reg
    end
    if ~exist(fullfile(rc.structOutput,'summary'))
        mkdir(fullfile(rc.structOutput,'summary'));
    end
    writetiff(pre_press_SIx_cat, fullfile(rc.structOutput,'summary',[mouse '_prepressAlignSIx.tif']));
    writetiff(pre_press_MIx_cat, fullfile(rc.structOutput,'summary',[mouse '_prepressAlignMIx.tif']));
    writetiff(pre_press_FIx_cat, fullfile(rc.structOutput,'summary',[mouse '_prepressAlignFIx.tif']));
    writetiff(post_press_SIx_cat, fullfile(rc.structOutput,'summary',[mouse '_postpressAlignSIx.tif']));
    writetiff(post_press_MIx_cat, fullfile(rc.structOutput,'summary',[mouse '_postpressAlignMIx.tif']));
    writetiff(post_press_FIx_cat, fullfile(rc.structOutput,'summary',[mouse '_postpressAlignFIx.tif']));
    writetiff(pre_target_SIx_cat, fullfile(rc.structOutput,'summary',[mouse '_pretargetAlignSIx.tif']));
    writetiff(pre_target_MIx_cat, fullfile(rc.structOutput,'summary',[mouse '_pretargetAlignMIx.tif']));
    writetiff(post_target_SIx_cat, fullfile(rc.structOutput,'summary',[mouse '_posttargetAlignSIx.tif']));
    writetiff(post_target_MIx_cat, fullfile(rc.structOutput,'summary',[mouse '_posttargetAlignMIx.tif']));
    writetiff(post_release_SIx_cat, fullfile(rc.structOutput,'summary',[mouse '_postreleaseAlignSIx.tif']));
    writetiff(post_release_FIx_cat, fullfile(rc.structOutput,'summary',[mouse '_postreleaseAlignFIx.tif']));
    return
end
            
            
            
            
            
            
            
