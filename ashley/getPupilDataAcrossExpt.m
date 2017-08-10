function pupilDataAllExpt = getPupilDataAcrossExpt(ms,fileSavePath)
%% Params - should put this in its own script
respWinDelayFr = 2;
respWinS = 0.25;

%%
nexp = size(ms,2);
delays = uniqueItemsInStructField(ms,'audDelayType');
ndelay = length(delays);
pupilDataAllExpt = struct;

for idelay = 1:ndelay
        pupilDataAllExpt(idelay).trType = delays{idelay};
        pupilDataAllExpt(idelay).exptsWithPupil = [];
        pupilDataAllExpt(idelay).pupilSizeAlignVisStimEaExpt = [];
        pupilDataAllExpt(idelay).pupilHorizPosAlignVisStimEaExpt = [];
        pupilDataAllExpt(idelay).pupilVertPosAlignVisStimEaExpt = [];
        pupilDataAllExpt(idelay).pupilSizeAlignStartEaExpt = [];
        pupilDataAllExpt(idelay).pupilHorizPosAlignStartEaExpt = [];
        pupilDataAllExpt(idelay).pupilVertPosAlignStartEaExpt = [];
        pupilDataAllExpt(idelay).neuronCorr2PupilSize = [];
        pupilDataAllExpt(idelay).pupilSizeResp = [];
        pupilDataAllExpt(idelay).pupilPupilHorizPosResp = [];
        pupilDataAllExpt(idelay).pupilPupilVertPosResp = [];
    for iexp = 1:nexp
        audDelayType_ind = strcmp(ms(iexp).audDelayType, delays{idelay});
        
        if ~any(audDelayType_ind) | isempty(ms(iexp).pupil{1}{1})
            continue
        else
            exptIdentityName = {[ms(iexp).subnum '-' ms(iexp).date]};
            pupilDataAllExpt(idelay).exptsWithPupil = cat(2,...
                pupilDataAllExpt(idelay).exptsWithPupil, exptIdentityName);
            
            pupilSizeAlignVisStim = ms(iexp).pupil{1}{audDelayType_ind};
            pupilHorizPosAlignVisStim = ms(iexp).pupil{2}{audDelayType_ind};
            pupilVertPosAlignVisStim = ms(iexp).pupil{3}{audDelayType_ind};
            
            pupilDataAllExpt(idelay).pupilSizeAlignVisStimEaExpt = ...
                cat(2,pupilDataAllExpt(idelay).pupilSizeAlignVisStimEaExpt,...
                mean(pupilSizeAlignVisStim,2));
            pupilDataAllExpt(idelay).pupilHorizPosAlignVisStimEaExpt = ...
                cat(2,pupilDataAllExpt(idelay).pupilHorizPosAlignVisStimEaExpt,...
                mean(pupilHorizPosAlignVisStim,2));
            pupilDataAllExpt(idelay).pupilVertPosAlignVisStimEaExpt = ...
                cat(2,pupilDataAllExpt(idelay).pupilVertPosAlignVisStimEaExpt,...
                mean(pupilVertPosAlignVisStim,2));
            
            if idelay <= ndelay-2
                pupilSizeAlignStart = ms(iexp).pupilLong{1}{audDelayType_ind};
                pupilHorizPosAlignStart = ms(iexp).pupilLong{2}{audDelayType_ind};
                pupilVertPosAlignStart = ms(iexp).pupilLong{3}{audDelayType_ind};

                pupilDataAllExpt(idelay).pupilSizeAlignStartEaExpt = ...
                    cat(2,pupilDataAllExpt(idelay).pupilSizeAlignStartEaExpt,...
                    mean(pupilSizeAlignStart,2));
                pupilDataAllExpt(idelay).pupilHorizPosAlignStartEaExpt = ...
                    cat(2,pupilDataAllExpt(idelay).pupilHorizPosAlignStartEaExpt,...
                    mean(pupilHorizPosAlignStart,2));            
                pupilDataAllExpt(idelay).pupilVertPosAlignStartEaExpt = ...
                    cat(2,pupilDataAllExpt(idelay).pupilVertPosAlignStartEaExpt,...
                    mean(pupilVertPosAlignStart,2));
            end
            
            respWindowFr = getResponseFrames(ms(iexp),respWinDelayFr,respWinS);
            
            pupilSizeRespEaTrial = mean(pupilSizeAlignVisStim(respWindowFr,:),1)';
            neuronRespEaTrial = squeeze(mean(ms(iexp).delayTrials...
                {audDelayType_ind}(respWindowFr,:,:),1))';
            neuronCorr2PupilSize = corr(neuronRespEaTrial,pupilSizeRespEaTrial)';
            pupilDataAllExpt(idelay).neuronCorr2PupilSize = cat(2,...
                pupilDataAllExpt(idelay).neuronCorr2PupilSize,neuronCorr2PupilSize);
            
            pupilSizeResp = nanmean(nanmean(pupilSizeAlignVisStim(respWindowFr,:)));
            pupilPupilHorizPosResp = nanmean(nanmean(pupilHorizPosAlignVisStim(respWindowFr,:)));
            pupilPupilVertPosResp = nanmean(nanmean(pupilVertPosAlignVisStim(respWindowFr,:)));
            
            pupilDataAllExpt(idelay).pupilSizeResp = cat(2,...
                pupilDataAllExpt(idelay).pupilSizeResp, pupilSizeResp);
            pupilDataAllExpt(idelay).pupilPupilHorizPosResp = cat(2,...
                pupilDataAllExpt(idelay).pupilPupilHorizPosResp, pupilPupilHorizPosResp);
            pupilDataAllExpt(idelay).pupilPupilVertPosResp = cat(2,...
                pupilDataAllExpt(idelay).pupilPupilVertPosResp, pupilPupilVertPosResp);
            
        end
        
        
    end
end


end

%% local functions used above
function respWindowFrames = getResponseFrames(ms,respWinDelayFr,respWinS)
    frRateHz = ms.frRateHz;
    preFrames = ms.pre_event_frames;
    respWindowFrames = preFrames + respWinDelayFr:preFrames + ... 
        respWinDelayFr + round(respWinS*frRateHz);
end
