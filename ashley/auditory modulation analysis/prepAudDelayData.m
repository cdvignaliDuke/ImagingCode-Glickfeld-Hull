function prepAudDelayData(ds)
rc = behavConstsAV;
eval(ds)
%%
mice = unique({expt.SubNum});
exptInd = zeros(1,length(mice));

%%
ms = struct;
for iexp = 1:size(expt,2)

subnum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

msInd = strcmp(mice,subnum);
exptInd(msInd) = exptInd(msInd)+1;

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
if ~exist(fullfile(fn,'data processing'),'dir')
    mkdir(fn,'data processing')
end
fnout = fullfile(fn,'data processing');
%%

load(fullfile(fnout,'timecourses_cells.mat'))
dataTC = data_tc_subnp;
[nfr,nCells] = size(dataTC);

nFrPerRun = zeros(1,expt(iexp).nrun);
for irun = 1:expt(iexp).nrun

    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
   
    input = loadMworksFile(subnum,expDate,expTime,rc.behavData);
    
    if strcmp(ds, 'audDelay_V1_EMX') & irun == 3 & strcmp(expDate,'180302')
        nFrPerRun(irun) = 31174;
    else
        load(fullfile(rc.ashleyData,mouse,'two-photon imaging',expDate,runFolder,fName))
        nFrPerRun(irun) = info.config.frames;
    end
    if irun == 1
        mworks = input;
    else
        mworks = [mworks input];
    end
    clear input
end
mworks = concatenateDataBlocks(mworks);
%%
frameRateHz = params.frameRate;
if frameRateHz ~= mworks.frameRateHz
    error(sprintf('Expt Frame Rate is %s',num2str(mworks.frameRateHz)));
end
stimOnTimeMs = params.stimOnTime;
if stimOnTimeMs ~= mworks.stimOnTimeMs
    error(sprintf('Vis stim time is %sms',num2str(mworks.stimOnTimeMs)))
end
nBaselineFr = params.nBaselineMs./1000*frameRateHz;
nStimFr_visAlign = params.nStimMs_visAlign./1000*frameRateHz;
nStimFr_audAlign = params.nStimMs_audAlign./1000*frameRateHz;
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

tPulseIntervalFr = targetStimOn - firstStimOn;
tPulseIntervalMs = round((tPulseIntervalFr./frameRateHz).*1000,3,'significant');

audVisDelayIntervalMs = tPulseIntervalMs;
audVisDelayIntervalMs(isVisAudPaired) = 0;

tOrientation = double(celleqel2mat_padded(mworks.tGratingDirectionDeg));

tAudStimOn = firstStimOn;
tVisStimOn = targetStimOn;
tVisStimOn(isVisAudPaired) = firstStimOn(isVisAudPaired);

%%
cmlvTrials = cumsum(mworks.trialsSinceReset);
cmlvFrames = cumsum(nFrPerRun);
for i = 1:expt(iexp).nrun
    if i == 1
        continue
    end
    tAudStimOn((cmlvTrials(i-1)+1):cmlvTrials(i)) = ...
        tAudStimOn((cmlvTrials(i-1)+1):cmlvTrials(i))+cmlvFrames(i-1);
    tVisStimOn( (cmlvTrials(i-1)+1):cmlvTrials(i)) = ...
        tVisStimOn( (cmlvTrials(i-1)+1):cmlvTrials(i))+cmlvFrames(i-1);
end

if any((tVisStimOn + nStimFr_audAlign) > nfr)
    ind = (tVisStimOn + nStimFr_audAlign) > nfr;
    tAudStimOn = tAudStimOn(~ind);
    tVisStimOn = tVisStimOn(~ind);
    tOrientation = tOrientation(~ind);
    audVisDelayIntervalMs = audVisDelayIntervalMs(~ind);
end
%% vis and aud aligned df/f
nTrials = length(tAudStimOn);
F_audAligned = nan(nBaselineFr+nStimFr_audAlign,nCells,nTrials);
for i = 1:nTrials
    ind = (tAudStimOn(i) +1 - nBaselineFr):(tAudStimOn(i)+nStimFr_audAlign);
    F_audAligned(:,:,i) = dataTC(ind,:);
end
F = mean(F_audAligned(1:nBaselineFr,:,:),1);
dFF_audAligned = (F_audAligned - F)./F;

F_visAligned = nan(nBaselineFr+nStimFr_visAlign,nCells,nTrials);
for i = 1:nTrials
    ind = (tVisStimOn(i) +1 - nBaselineFr):(tVisStimOn(i)+nStimFr_visAlign);
    F_visAligned(:,:,i) = dataTC(ind,:);
end
F = mean(F_visAligned(1:nBaselineFr,:,:),1);
dFF_visAligned = (F_visAligned - F)./F;

%% take out motion trials, cutoff trials

motionInd = max(diff(squeeze(mean(dFF_audAligned,2)))) > params.motionCutoff ...
    | max(diff(squeeze(mean(dFF_visAligned,2)))) > params.motionCutoff;


tOrientation = tOrientation(~motionInd);
audVisDelayIntervalMs = audVisDelayIntervalMs(~motionInd);

dFF_audAligned = dFF_audAligned(:,:,~motionInd);
dFF_visAligned = dFF_visAligned(:,:,~motionInd);
%% organize by expt

ms(msInd).expt(exptInd(msInd)).info.name = subnum;
ms(msInd).expt(exptInd(msInd)).info.date = expDate;
ms(msInd).expt(exptInd(msInd)).info.frameRateHz = frameRateHz;
ms(msInd).expt(exptInd(msInd)).info.stimOnTimeMs = stimOnTimeMs;

ms(msInd).expt(exptInd(msInd)).tOrientation = tOrientation;
ms(msInd).expt(exptInd(msInd)).tAudVisDelay = audVisDelayIntervalMs;
ms(msInd).expt(exptInd(msInd)).tVisAlignResp = dFF_visAligned;
ms(msInd).expt(exptInd(msInd)).tAudAlignResp = dFF_audAligned;

%%
save(fullfile(rc.ashleyAnalysis,'Expt Summaries',ds,'dataSummary'),'ms')
end
end
