function extractAllNeededWFTCs(mouse)
while true
    rc = behavConstsWF(mouse);
    [xd indexRowN] = findNextIndexForWFTC(rc);
    
    if isempty(indexRowN)
        disp('No more records to do');
        dateStr = [];
        return
    else
        % do the one we found
        dateStr = xd.DateStr{indexRowN};
        fprintf('Starting subject %s, date %s\n', ...
        mouse, dateStr);
        data_folder = fullfile(rc.dataRootDir, [dateStr '_' mouse]);
        nrun = str2num(xd.Runs{indexRowN});
        data = [];
        for irun = 1:nrun
            data_temp = readtiff(fullfile(data_folder,[dateStr '_' mouse '_' num2str(irun)], [dateStr '_' mouse '_' num2str(irun) '_MMStack.ome.tif']));
            data = cat(3, data, data_temp);
            if irun == 1
                timeStr = num2str(eval(['xd.MatFileRun' num2str(irun) '(indexRowN)']));
                fn_mworks = fullfile(rc.pathStr,['data-' mouse '-' dateStr '-' timeStr '.mat']);
                input = mwLoadData(fn_mworks, [], []);
            end
        end
        clear data_temp
        
        if indexRowN == 1
            tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
            tGratingContrast = celleqel2mat_padded(input.tGratingContrast);
            cStimOn = celleqel2mat_padded(input.cStimOn);
            cDecision = celleqel2mat_padded(input.cDecision);
            SIx = strcmp(input.trialOutcomeCell, 'success');
            nTrials = length(tLeftTrial);
            sz = size(data);
            data_f = nan(sz(1),sz(2),nTrials);
            data_targ = nan(sz(1),sz(2),nTrials);
            for itrial = 1:nTrials
                if cStimOn(itrial)+25<sz(3)
                    data_f(:,:,itrial) = mean(data(:,:,cStimOn(itrial)-20:cStimOn(itrial)-1),3);
                    data_targ(:,:,itrial) = mean(data(:,:,cStimOn(itrial)+5:cStimOn(itrial)+25),3);
                end
            end
            data_targ_dfof = (data_targ-data_f)./data_f;
            indL = intersect(find(SIx),intersect(find(tGratingContrast==1),find(tLeftTrial)));
            data_dfof_L = mean(data_targ_dfof(:,:,indL),3);
            indR = intersect(find(SIx),intersect(find(tGratingContrast==1),find(~tLeftTrial)));
            data_dfof_R = mean(data_targ_dfof(:,:,indR),3);
            avg_resp = data_dfof_R-data_dfof_L;
            figure; imagesc(avg_resp); colormap(gray);
            clim([-0.02 .05])
            %adjust number of ROIs to choose
            nROI = 4;
            for i = 1:nROI
                roi(i) = impoly;
            end

            siz = size(avg_resp);
            mask = zeros(siz(1),siz(2),nROI);
            for i = 1:nROI
                mask(:,:,i) = createMask(roi(i));
            end

            roi_cluster = sum(mask,3);
            mask_cell = bwlabel(roi_cluster);
            figure; imagesc(mask_cell);
            
            area_list = strvcat('V1','LM','AL','PM');
            save(fullfile(rc.structOutput, [mouse '_roi_masks.mat']), 'roi_cluster', 'mask_cell', 'area_list', 'data_dfof_R', 'data_dfof_L', 'avg_resp');
            roiTC = stackGetTimeCourses(data, mask_cell);
            save(fullfile(rc.structOutput, [mouse '_' dateStr '_tc.mat']), 'roiTC');
        else
            load(fullfile(rc.structOutput, [mouse '_' xd.DateStr{1} '_roi_masks.mat']))
            load(fullfile(rc.structOutput, [mouse '_' xd.DateStr{indexRowN} '_reg_data.mat']))
            %register
            sz_target  = size(roi_avg);
            mytform    = maketform('affine',input_points(1:3,:), base_points(1:3,:));
            registered = imtransform(data,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]); 
            roiTC = stackGetTimeCourses(registered, mask_cell);
            save(fullfile(rc.structOutput, [mouse '_' dateStr '_tc.mat']), 'roiTC');
        end
    end
end            
            
