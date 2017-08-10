function speedDataAllExpt = getSpeedDataAcrossExpt(ms,fileSavePath)
%% Params - should put this in its own script
respWinDelayS = 0.067;
respWinDelayFr = 2;
respWinS = 0.25;

%%
nexp = size(ms,2);
delays = uniqueItemsInStructField(ms,'audDelayType');
ndelay = length(delays);
speedDataAllExpt = struct;

for idelay = 1:ndelay
        speedDataAllExpt(idelay).trType = delays{idelay};
        speedDataAllExpt(idelay).exptsWithSpeed = [];
        speedDataAllExpt(idelay).speedAlignVisStimEaExpt = [];
        speedDataAllExpt(idelay).speedAlignVisStimErrEaExpt = [];
        speedDataAllExpt(idelay).speedAlignStartEaExpt = [];
        speedDataAllExpt(idelay).speedAlignStartErrEaExpt = [];
        speedDataAllExpt(idelay).neuronCorr2Speed = [];
        speedDataAllExpt(idelay).speedResp = [];   
        speedDataAllExpt(idelay).speedRespErr = [];
        speedDataAllExpt(idelay).speedAlignVisStimTimeStampS = [];
        speedDataAllExpt(idelay).speedAlignStartTimeStampS = [];
    for iexp = 1:nexp
        audDelayType_ind = strcmp(ms(iexp).audDelayType, delays{idelay}); 
        if ~any(audDelayType_ind) | isempty(ms(iexp).trSpeed)
            nCells = size(ms(iexp).delayTrials{1},2);
            nTimePointsAlignVisStim = 20;
            nTimePointsAlignStart = 90;
            speedDataAllExpt(idelay).speedAlignVisStimEaExpt = cat(2,...
                speedDataAllExpt(idelay).speedAlignVisStimEaExpt,...
                nan(nTimePointsAlignVisStim,1));
            speedDataAllExpt(idelay).speedAlignVisStimErrEaExpt = cat(2,...
                speedDataAllExpt(idelay).speedAlignVisStimErrEaExpt,...
                nan(nTimePointsAlignVisStim,1));
            speedDataAllExpt(idelay).speedAlignStartEaExpt = cat(2,...
                speedDataAllExpt(idelay).speedAlignStartEaExpt,...
                nan(nTimePointsAlignStart,1));
            speedDataAllExpt(idelay).speedAlignStartErrEaExpt = cat(2,...
                speedDataAllExpt(idelay).speedAlignStartErrEaExpt,...
                nan(nTimePointsAlignStart,1));
            speedDataAllExpt(idelay).neuronCorr2Speed = cat(2,...
                speedDataAllExpt(idelay).neuronCorr2Speed,...
                nan(1,nCells));
            speedDataAllExpt(idelay).speedResp = cat(2,...
                speedDataAllExpt(idelay).speedResp, NaN);
            speedDataAllExpt(idelay).speedRespErr = cat(2,...
                speedDataAllExpt(idelay).speedRespErr,NaN);
        else
            exptIdentityName = {[ms(iexp).subnum '-' ms(iexp).date]};
            speedDataAllExpt(idelay).exptsWithSpeed = cat(2,...
                speedDataAllExpt(idelay).exptsWithSpeed, exptIdentityName);
           
            speedAlignVisStim = getNormedTimeCourse(...
                ms(iexp).trSpeed{audDelayType_ind},ms(iexp).pre_event_wheel,...
                'subtraction');
            speedDataAllExpt(idelay).speedAlignVisStimEaExpt = ...
                cat(2,speedDataAllExpt(idelay).speedAlignVisStimEaExpt,...
                mean(speedAlignVisStim,2));
            speedDataAllExpt(idelay).speedAlignVisStimErrEaExpt = ...
                cat(2,speedDataAllExpt(idelay).speedAlignVisStimErrEaExpt,...
                ste(speedAlignVisStim,2));  
            if isempty(speedDataAllExpt(idelay).speedAlignVisStimTimeStampS)
                speedAlignVisStimTimeStampS = double(-ms(iexp).pre_event_wheel+1:...
                    size(speedAlignVisStim,1)-ms(iexp).pre_event_wheel)...
                        /double(ms(iexp).wheelRateHz);
                speedDataAllExpt(idelay).speedAlignVisStimTimeStampS = cat(2,...
                    speedDataAllExpt(idelay).speedAlignVisStimTimeStampS,...
                    speedAlignVisStimTimeStampS);
            end
            
            if idelay <= ndelay-2
                speedAlignStart = getNormedTimeCourse(...
                    ms(iexp).trSpeedLong{audDelayType_ind},ms(iexp).pre_event_wheel,...
                    'subtraction');
                speedDataAllExpt(idelay).speedAlignStartEaExpt = ...
                    cat(2,speedDataAllExpt(idelay).speedAlignStartEaExpt,...
                    mean(speedAlignStart,2));
                speedDataAllExpt(idelay).speedAlignStartErrEaExpt = ...
                    cat(2,speedDataAllExpt(idelay).speedAlignStartErrEaExpt,...
                    ste(speedAlignStart,2));
                if isempty(speedDataAllExpt(idelay).speedAlignStartTimeStampS)
                    speedAlignStartTimeStampS = double(-ms(iexp).pre_event_wheel+1:...
                    size(speedAlignStart,1)-ms(iexp).pre_event_wheel)...
                        /double(ms(iexp).wheelRateHz);
                    speedDataAllExpt(idelay).speedAlignStartTimeStampS = cat(2,...
                        speedDataAllExpt(idelay).speedAlignStartTimeStampS,...
                        speedAlignStartTimeStampS);
                end
            end 
            
            respWindowFr_wheel = getResponseWheelFrames(ms(iexp),respWinDelayS,respWinS);
            respWindowFr_image = getResponseImgFrames(ms(iexp),respWinDelayFr,respWinS);
            speedRespEaTrial = mean(speedAlignVisStim(respWindowFr_wheel,:),1)';
            neuronRespEaTrial = squeeze(mean(ms(iexp).delayTrials...
                {audDelayType_ind}(respWindowFr_image,:,:),1))';
            neuronCorr2Speed = corr(neuronRespEaTrial,speedRespEaTrial)';
            speedDataAllExpt(idelay).neuronCorr2Speed = cat(2,...
                speedDataAllExpt(idelay).neuronCorr2Speed,neuronCorr2Speed);
            
            speedResp = mean(mean(speedAlignVisStim(respWindowFr_wheel,:)));
            speedRespErr = ste(mean(speedAlignVisStim(respWindowFr_wheel,:),1),2);
            
            speedDataAllExpt(idelay).speedResp = cat(2,...
                speedDataAllExpt(idelay).speedResp, speedResp);
            speedDataAllExpt(idelay).speedRespErr = cat(2,...
                speedDataAllExpt(idelay).speedRespErr,speedRespErr);
        end
        
        
    end
end


end

%% local functions used above
function respWindowFrames = getResponseWheelFrames(ms,respWinDelayS,respWinS)
    wheelRateHz = double(ms.wheelRateHz);
    delayFrames = round(respWinDelayS*wheelRateHz);
    preFrames = double(ms.pre_event_wheel);
    respWindowFrames = preFrames + delayFrames:preFrames + ... 
        delayFrames + round(respWinS*wheelRateHz);
end

function respWindowFrames = getResponseImgFrames(ms,respWinDelayFr,respWinS)
    frRateHz = ms.frRateHz;
    preFrames = ms.pre_event_frames;
    respWindowFrames = preFrames + respWinDelayFr:preFrames + ... 
        respWinDelayFr + round(respWinS*frRateHz);
end

