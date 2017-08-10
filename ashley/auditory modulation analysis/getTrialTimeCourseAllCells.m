function neuronTimeCourse = getTrialTimeCourseAllCells(ms,fileSavePath)

%% loop through each mouse
delays = uniqueItemsInStructField(ms,'audDelayType');
ndelay = length(delays);
oris = uniqueItemsInStructField(ms,'orientations');
nori = length(oris);

nexp = size(ms,2);
neuronTimeCourse = struct;
for idelay = 1:ndelay
    neuronTimeCourse(idelay).trType = delays{idelay};
    neuronTimeCourse(idelay).resp2VisStim = [];
    neuronTimeCourse(idelay).resp2VisStimErr = [];
    neuronTimeCourse(idelay).resp2VisStimTimeS = [];
    neuronTimeCourse(idelay).longResp2TrialStart = [];
    neuronTimeCourse(idelay).longResp2TrialStartErr = [];
    neuronTimeCourse(idelay).longTimeCourseTimeS = [];   
    for iexp = 1:nexp
        audDelayType_ind = strcmp(ms(iexp).audDelayType, delays{idelay});
        [nFramesResp,nCells,~] = size(ms(iexp).delayTrials{1}); 
        nFramesFromTrialStart = size(ms(iexp).delayLong{1},1);
        if ~any(audDelayType_ind)
            neuronTimeCourse(idelay).resp2VisStim = cat(2, ...
                neuronTimeCourse(idelay).resp2VisStim, nan(nFramesResp,nCells,nori));
            neuronTimeCourse(idelay).resp2VisStimErr = cat(2, ...
                neuronTimeCourse(idelay).resp2VisStimErr, nan(nFramesResp,nCells,nori));
            neuronTimeCourse(idelay).longResp2TrialStart = cat(2, ...
                neuronTimeCourse(idelay).longResp2TrialStart, ...
                nan(nFramesFromTrialStart,nCells,nori));
            neuronTimeCourse(idelay).longResp2TrialStartErr = cat(2, ...
                neuronTimeCourse(idelay).longResp2TrialStartErr, ...
                nan(nFramesFromTrialStart,nCells,nori));
        else
            % get imaging data for this trial type          
            neuronResp2VisStim = ms(iexp).delayTrials{audDelayType_ind};

            % imaging params for this experiment
            frRateHz = ms(iexp).frRateHz;
            preVisStimFrames = ms(iexp).pre_event_frames;
            nTimePointsResp2VisStim = size(neuronResp2VisStim,1);         

            % trial info
            trOris = ms(iexp).trOri{audDelayType_ind};
            ori_ind = findgroups(trOris);
            nori = length(oris);

            %response to each orientation
            resp2VisStim = nan(nTimePointsResp2VisStim,nCells,nori);
            resp2VisStimErr = nan(nTimePointsResp2VisStim,nCells,nori);
            for iori = 1:nori                
                resp2VisStim(:,:,iori) = mean(neuronResp2VisStim(...
                    :,:,ori_ind == iori),3);
                resp2VisStimErr(:,:,iori) = ste(neuronResp2VisStim(...
                    :,:,ori_ind == iori),3);
            end

            neuronTimeCourse(idelay).resp2VisStim = cat(2, ...
                neuronTimeCourse(idelay).resp2VisStim,resp2VisStim);
            neuronTimeCourse(idelay).resp2VisStimErr = cat(2, ...
                neuronTimeCourse(idelay).resp2VisStimErr,resp2VisStimErr);

            if isempty(neuronTimeCourse(idelay).resp2VisStimTimeS)
                resp2VisStimTimeS = ...
                    [-preVisStimFrames+1:nTimePointsResp2VisStim-preVisStimFrames]...
                    /frRateHz;
                neuronTimeCourse(idelay).resp2VisStimTimeS = cat(2,...
                    neuronTimeCourse(idelay).resp2VisStimTimeS,...
                    resp2VisStimTimeS);
            end 

            if idelay <= ndelay-2
                neuronLongResp2TrialStart = ms(iexp).delayLong{audDelayType_ind};
                nTimePointsResp2TrialStart = size(neuronLongResp2TrialStart,1);
                longResp2TrialStart = nan(nTimePointsResp2TrialStart,nCells,nori);
                longResp2TrialStartErr = nan(nTimePointsResp2TrialStart,nCells,nori);
                for iori = 1:nori               
                    longResp2TrialStart(:,:,iori) = mean(...
                        neuronLongResp2TrialStart(:,:,ori_ind == iori),3);
                    longResp2TrialStartErr(:,:,iori) = ste(...
                        neuronLongResp2TrialStart(:,:,ori_ind == iori),3);
                end
                neuronTimeCourse(idelay).longResp2TrialStart = cat(2, ...
                    neuronTimeCourse(idelay).longResp2TrialStart,longResp2TrialStart);
                neuronTimeCourse(idelay).longResp2TrialStartErr = cat(2, ...
                    neuronTimeCourse(idelay).longResp2TrialStartErr,longResp2TrialStartErr);
                if isempty(neuronTimeCourse(idelay).longTimeCourseTimeS)            
                    longTimeCourseTimeS = ...
                        [-preVisStimFrames+1:nTimePointsResp2TrialStart-preVisStimFrames]...
                        /frRateHz;
                    neuronTimeCourse(idelay).longTimeCourseTimeS = cat(2,...
                        neuronTimeCourse(idelay).longTimeCourseTimeS,...
                        longTimeCourseTimeS);                
                end
            end
        end
    end
end
% save(fullfile(fileSavePath,'neuronTrialTimeCourses.mat'),'neuronTimeCourse')
