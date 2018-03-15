function maxDFF_ori = getAudDelayMaxDFF(params,mworks,data_reg,nFrPerRun)

frameRateHz = params.frameRate;
nBaselineFr = params.nBaselineMs./1000*frameRateHz;
[ypix,xpix,nfr] = size(data_reg);

motionCutoff = params.motionCutoff;
visRespS = 0.2;
nStimS = 0.5;
nStimFr = nStimS*frameRateHz;

visStimRespWin = nBaselineFr+params.nFramesVisDelay:(nBaselineFr+params.nFramesVisDelay+(0.2*frameRateHz));
%%
blk2 = cell2mat(mworks.tBlock2TrialNumber);
blk2BaseCon = mworks.block2BaseGratingContrast;
blk1BaseCon = mworks.baseGratingContrast;
if blk2BaseCon == 1
    pairedBlock = 1;
elseif blk1BaseCon == 1
    pairedBlock = 0;
elseif blk2BaseCon == blk1BaseCon
    error('Blocks 1 and 2 are the same')
end
isVisAudPaired = blk2 == pairedBlock;

firstStimOn = double(cell2mat(mworks.cFirstStim));
targetStimOn = double(cell2mat(mworks.cTargetOn));

tOrientation = double(celleqel2mat_padded(mworks.tGratingDirectionDeg));
orientations = unique(tOrientation);

tVisStimOn = targetStimOn;
tVisStimOn(isVisAudPaired) = firstStimOn(isVisAudPaired);

cmlvTrials = cumsum(mworks.trialsSinceReset);
cmlvFrames = cumsum(nFrPerRun);
for i = 1:length(nFrPerRun)
    if i == 1
        continue
    end
    tVisStimOn( (cmlvTrials(i-1)+1):cmlvTrials(i)) = ...
        tVisStimOn( (cmlvTrials(i-1)+1):cmlvTrials(i))+cmlvFrames(i-1);
end

if any((tVisStimOn + (4000./frameRateHz)) > nfr)
    ind = (tVisStimOn + (4000./frameRateHz)) > nfr;
    tVisStimOn = tVisStimOn(~ind);
    tOrientation = tOrientation(~ind);
end

%%
nTrials = length(tVisStimOn);
F_visAligned = nan(ypix,xpix,nBaselineFr+nStimFr,nTrials);
for i = 1:nTrials
    ind = (tVisStimOn(i) +1 - nBaselineFr):(tVisStimOn(i)+nStimFr);
    F_visAligned(:,:,:,i) = data_reg(:,:,ind);
end

F = mean(F_visAligned(:,:,1:nBaselineFr,:),3);
dFF = (F_visAligned - F)./F;
clear F_visAligned
%%
dFF_resp = squeeze(mean(dFF(:,:,visStimRespWin,:),3));
maxDFF = max(dFF_resp,[],3);
clear dFF

bwout = imCellEditInteractive(maxDFF);
mask = bwlabel(bwout);
tc = stackGetTimeCourses(data_reg,mask);
nROI = length(unique(mask))-1;

tc_visAligned = nan(nBaselineFr+nStimFr,nROI,nTrials);
for i = 1:nTrials
    ind = (tVisStimOn(i) +1 - nBaselineFr):(tVisStimOn(i)+nStimFr);
    tc_f = tc(ind,:);
    tc_f0 = mean(tc_f(1:nBaselineFr,:),1);
    tc_dff = (tc_f - tc_f0)./tc_f0;
    tc_visAligned(:,:,i) = tc_dff;
end

motionInd = max(diff(squeeze(mean(tc_visAligned,2)))) > motionCutoff;

dFF_resp = dFF_resp(:,:,~motionInd);
tOrientation = tOrientation(~motionInd);
%%
nOris = length(orientations);
maxDFF_ori = nan(ypix,xpix,nOris);
for i = 1:nOris
    ind = tOrientation == orientations(i);
    maxDFF_ori(:,:,i) = max(dFF_resp(:,:,ind),[],3);
end

end