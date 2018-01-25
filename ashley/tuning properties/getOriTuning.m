function [avgResponseEaOri,semResponseEaOri,vonMisesFitAllCells,fitReliability,R_square,tuningTC] = ...
    getOriTuning(tc,mworks,downSampleFactor)
basewin = 8:12;
respwin = 13:16;
nBoot = 1000;


[nFrames,nCells] = size(tc);

%params
nOn = (mworks.nScansOn)./downSampleFactor;
nOff = (mworks.nScansOff)./downSampleFactor;
nTrials = mworks.trialSinceReset;
if mod(nFrames,nTrials) > 0
    nframesAllTrials = nTrials*(nOn+nOff);
    if nframesAllTrials > nFrames
        nTrials = floor(nFrames/(nOn+nOff));
    end 
end
tc = tc(1:(nOn+nOff)*nTrials,:);
tDirection = cell2mat(mworks.tGratingDirectionDeg);
tDirection = tDirection(1:nTrials);
tOrientation = tDirection;
tOrientation(tOrientation > 179) = tOrientation(tOrientation > 179) - 180;
[orientationInd, orientations] = findgroups(tOrientation);
nStim = length(orientations);
% theta = [orientations 180];

trialTC = reshape(tc,nOff+nOn,nTrials,nCells);
F = mean(trialTC(basewin,:,:),1);
dFF = bsxfun(@rdivide, bsxfun(@minus,trialTC,F),F);

tuningTC = nan(nOn+nOff,nCells,nStim);
avgResponseEaOri = nan(nCells,nStim);
semResponseEaOri = nan(nCells,nStim);
tuningResamp = nan(nCells,nStim,nBoot);
for istim = 1:nStim
    ind = find(orientationInd == istim);
    tuningTC(:,:,istim) = squeeze(mean(dFF(:,ind,:),2));
    avgResponseEaOri(:,istim) = squeeze(mean(mean(dFF(respwin,ind,:),1),2));
    semResponseEaOri(:,istim) = squeeze(ste(mean(dFF(respwin,ind,:),1),2));
    for iboot = 1:nBoot
        n = length(ind);
        randTrials = randsample(ind,n,1);
        tuningResamp(:,istim,iboot) = squeeze(mean(mean(...
            dFF(respwin,randTrials,:),1),2));
    end
end
% tuningResamp4Fit = cat(2,tuningResamp,tuningResamp(:,1,:));
% avgResp4Fit = cat(2,avgResponseEaOri,avgResponseEaOri(:,1));
% orientations = circshift(orientations,1);
[vonMisesFitAllCells,~,fitReliability,R_square] = vonmisesReliableFit(avgResponseEaOri,...
    tuningResamp,orientations,nBoot);

end