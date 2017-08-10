function neuronTuning = getOriTuningAllCells(ms,fileSavePath)
% function to summarize orientation tuning for each delay for each
% neuron, output as a structure, save the structure.

%% Params
respWinDelayFr = 2;
respWinS = 0.25;
nboot = 10;
reliableFitCutoff = 30;

%% loop through each mouse
delays = uniqueItemsInStructField(ms,'audDelayType');
ndelay = length(delays);

nexp = size(ms,2);
neuronTuning = struct;
for idelay = 1:ndelay
    % define all fields to be concatenated across experiments
    neuronTuning(idelay).trType = delays{idelay};
    neuronTuning(idelay).exptIdentity = [];
    neuronTuning(idelay).relaibleResponseCurve = [];
    neuronTuning(idelay).rawTuningCurve = [];
    neuronTuning(idelay).rawTuningCurveErr = [];
    neuronTuning(idelay).rawPrefInd = [];
    neuronTuning(idelay).rawOrthInd = [];
    neuronTuning(idelay).reliableFits = [];
    neuronTuning(idelay).fitTuningCurve = [];
    neuronTuning(idelay).fitPrefInd = [];
    neuronTuning(idelay).fitOrthInd = [];
    neuronTuning(idelay).isLabelled = [];
    neuronTuning(idelay).cellType = [];
    for iexp = 1:nexp
        % get imaging data for this trial type
        nCells = size(ms(iexp).delayTrials{1},2);
        nOrientations = length(ms(iexp).orientations);
        orisSmooth = 0:1:180;
        % save experiment identity
        exptIdentityName = {[ms(iexp).subnum '-' ms(iexp).date]};
        neuronTuning(idelay).exptIdentity = cat(2,...
            neuronTuning(idelay).exptIdentity, repmat(exptIdentityName,1,nCells));

        audDelayType_ind = strcmp(ms(iexp).audDelayType, delays{idelay});
        if ~any(audDelayType_ind)
            neuronTuning(idelay).relaibleResponseCurve = cat(2, ...
                neuronTuning(idelay).relaibleResponseCurve, ...
                nan(nOrientations,nCells));
            neuronTuning(idelay).rawTuningCurve = cat(2,... 
                neuronTuning(idelay).rawTuningCurve,nan(nOrientations,nCells));
            neuronTuning(idelay).rawTuningCurveErr = cat(2, ...
                neuronTuning(idelay).rawTuningCurveErr, nan(nOrientations,nCells));
            neuronTuning(idelay).rawPrefInd = cat(2,... 
                neuronTuning(idelay).rawPrefInd, nan(1,nCells));
            neuronTuning(idelay).rawOrthInd = cat(2,... 
                neuronTuning(idelay).rawOrthInd, nan(1,nCells));
            neuronTuning(idelay).reliableFits = cat(2,... 
                neuronTuning(idelay).reliableFits, nan(1,nCells));
            neuronTuning(idelay).fitTuningCurve = cat(2,... 
                neuronTuning(idelay).fitTuningCurve, nan(length(orisSmooth),nCells));
            neuronTuning(idelay).fitPrefInd = cat(2,... 
                neuronTuning(idelay).fitPrefInd, nan(1,nCells));
            neuronTuning(idelay).fitOrthInd = cat(2,... 
                neuronTuning(idelay).fitOrthInd, nan(1,nCells));
            neuronTuning(idelay).isLabelled = cat(2,... 
                neuronTuning(idelay).isLabelled, nan(1,nCells));
            neuronTuning(idelay).cellType = cat(2,... 
                neuronTuning(idelay).cellType, nan(1,nCells));        
        else
            neuronActivity = ms(iexp).delayTrials{audDelayType_ind};

            % orientation tuning
            respWinFr = getResponseFrames(ms(iexp),respWinDelayFr,respWinS);
            trOris = ms(iexp).trOri{audDelayType_ind};
            [ori_ind, oris] = findgroups(trOris);
            nOri = length(oris);
            orisFit = [oris 180];
            orisSmooth = 0:1:180;

                % raw data
            neuronTrialResp = squeeze(mean(neuronActivity(respWinFr,:,:),1));
            neuronTuningExpt = splitapply(@(x) mean(x,2),neuronTrialResp,ori_ind)';
            neuronTuningErrExpt = splitapply(@(x) ste(x,2),neuronTrialResp,ori_ind)';
            [~, neuronMaxRespOri_ind] = max(neuronTuningExpt,[],1);
            neuronOrth2Pref_ind = getOrthogonalInd(neuronMaxRespOri_ind,oris);
            
                % reliable repsonses
            neuronReliableResponses = zeros(nOri,nCells);
            for iOri = 1:nOri
                thisOriInd = ori_ind == iOri;
                resp2Ori = neuronTrialResp(:,thisOriInd)';
                resp2OriTtest = ttest(resp2Ori,[],'tail','right','alpha',0.05/nOri);
                neuronReliableResponses(iOri,:) = resp2OriTtest;
            end
            reliableResponsesVector = neuronReliableResponses(:);
            rawTuningData4FittingVector = neuronTuningExpt(:);
            rawTuningData4FittingVector(reliableResponsesVector == 0) = 0;
            rawTuningDataRect = reshape(rawTuningData4FittingVector,nOri,nCells);
            
            neuronTuning(idelay).relaibleResponseCurve = cat(2, ...
                neuronTuning(idelay).relaibleResponseCurve, ...
                neuronReliableResponses);
            neuronTuning(idelay).rawTuningCurve = cat(2, ...
                neuronTuning(idelay).rawTuningCurve, neuronTuningExpt);
            neuronTuning(idelay).rawTuningCurveErr = cat(2, ...
                neuronTuning(idelay).rawTuningCurveErr, neuronTuningErrExpt);
            neuronTuning(idelay).rawPrefInd = cat(2, ...
                neuronTuning(idelay).rawPrefInd, neuronMaxRespOri_ind);
            neuronTuning(idelay).rawOrthInd = cat(2, ...
                neuronTuning(idelay).rawOrthInd, neuronOrth2Pref_ind);

                % fit data
            neuronTuningCircular = cat(1, neuronTuningExpt, neuronTuningExpt(1,:))';
            if strcmp(delays{idelay},'vis only')
                neuronTuningBoot = getBootstrappedTuningCurves(neuronTrialResp,...
                    ori_ind, nboot);
                neuronTuningCircularBoot = permute(...
                    cat(1,neuronTuningBoot,neuronTuningBoot(1,:,:)),[2,1,3]);
                [~, ~, distanceOf90PctFits] = vonmisesReliableFit(...
                    neuronTuningCircular, neuronTuningCircularBoot,orisFit,nboot);
                reliableFits = distanceOf90PctFits < reliableFitCutoff;
                neuronTuning(idelay).reliableFits = cat(2,...
                    neuronTuning(idelay).reliableFits, reliableFits);
            else                
                neuronTuning(idelay).reliableFits = cat(2,...
                    neuronTuning(idelay).reliableFits, nan(1,nCells));
            end
            fitTuningCurves = getVonMisesFits(neuronTuningCircular,orisFit,orisSmooth);
            neuronTuning(idelay).fitTuningCurve = cat(2,...
                neuronTuning(idelay).fitTuningCurve, fitTuningCurves);
            [~, neuronPeakOfFit_ind] = max(fitTuningCurves,[],1);
            neuronOrth2PeakFit_ind = getOrthogonalInd(neuronPeakOfFit_ind,orisSmooth);
            neuronTuning(idelay).fitPrefInd = cat(2, ...
                neuronTuning(idelay).fitPrefInd, neuronPeakOfFit_ind);
            neuronTuning(idelay).fitOrthInd = cat(2, ...
                neuronTuning(idelay).fitOrthInd, neuronOrth2PeakFit_ind);

            % make a logical vector indicating whether a cell is labelled &
            % make a cell array of strings indicating label type
            labelledCells = false(1,nCells);
            cellType = cell(1,nCells);
            cellType(:) = {nan};
            if ms(iexp).redLabel
                labelledCells_ind = ms(iexp).redCellIdentity(...
                    ms(iexp).redCellIdentity ~= 0);
                labelledCells(labelledCells_ind) = true;

                cellType(labelledCells) = {ms(iexp).redCellType};  
            end
            neuronTuning(idelay).isLabelled = cat(2, ... 
                neuronTuning(idelay).isLabelled, labelledCells);        
            neuronTuning(idelay).cellType = cat(2, ... 
                neuronTuning(idelay).cellType, cellType);
        end
    end
end
% save(fullfile(fileSavePath,'oriTuningAndFits.mat'),'neuronTuning')
%% local functions used above
function respWindowFrames = getResponseFrames(ms,respWinDelayFr,respWinS)
    frRateHz = ms.frRateHz;
    preFrames = ms.pre_event_frames;
    respWindowFrames = preFrames + respWinDelayFr:preFrames + ... 
        respWinDelayFr + round(respWinS*frRateHz);
end

function orthogonal_ind = getOrthogonalInd(prefferred_ind,oris)
    nori = length(oris);
    orthogonal_ind = prefferred_ind + (nori/2);
    orthogonal_ind(orthogonal_ind > nori) = ... 
        prefferred_ind(orthogonal_ind > nori) - (nori/2);
end

function tuningBoots = getBootstrappedTuningCurves(neuronTrialResp, trOri_ind, nBoot)
    nori = length(unique(trOri_ind));
    tuningBoots = nan(nori,size(neuronTrialResp,1),nBoot);
    for iori = 1:nori
        ind = trOri_ind == iori;
        for iboot = 1:nBoot
            boot_ind = randsample(find(ind), sum(ind),1);
            neuronBootTuning = mean(neuronTrialResp(:,boot_ind),2)';
            tuningBoots(iori,:,iboot) = neuronBootTuning;
        end
    end
end

function y_fit = getVonMisesFits(neuronTuningCurves,orisFit,orisSmooth)
    neuronTuningCurves = neuronTuningCurves';
    nc = size(neuronTuningCurves,2);
    y_fit = zeros(length(orisSmooth),nc);
    for icell = 1:nc
        [b, k, R, u] = miaovonmisesfit_ori(deg2rad(orisFit),neuronTuningCurves(:,icell));
        y_fit(:,icell) = b+R.*exp(k.*(cos(2.*(deg2rad(orisSmooth)-u))-1));
    end
end

end