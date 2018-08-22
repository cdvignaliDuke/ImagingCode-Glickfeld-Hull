clear all
close all
ds = 'adaptDREADDS_V1_SOM_control';
doRedOnly = 0;
motionThreshold = 0.05;
slct_expt = 1:4;
%%

rc = behavConstsAV;
eval(ds)

nBaselineFr = (params.nBaselineMs/1000)*params.frameRate;
nRespwinFr = 0.1*params.frameRate;
nPostStimFr = (params.nPostStimMs/1000)*params.frameRate;
nTrialFrames = nBaselineFr+nPostStimFr;
%%

for iexp = slct_expt
    mouse = expt(iexp).mouse;
    subnum = mouse;
    expDate = expt(iexp).date;
    folder = expt(iexp).adapt{1};
    t = expt(iexp).adapt{2};
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    disp([mouse ' ' expDate])
    
    %% adapt
    % load data
    load(fullfile(fn,'data processing','redTC.mat'));
    redTC = redTC{1};
    load(fullfile(fn,'data processing','dataTC.mat'));
    greenTC= dataTC{1};
    mw = loadMworksFile(subnum,expDate,t);

    secondStimDelay = cell2mat(mw.tStimOffTimeMs);
            
    trialStartInd = cell2mat(mw.cFirstStim);
    secondStimInd = cell2mat(mw.cTargetOn);
    [nFr, nRedCells] = size(redTC);
    [~, nGreenCells] = size(greenTC);
    nTrials = length(trialStartInd);
    
    greenDFF = cell(1,2);
    redDFF = cell(1,2);

    greenF = nan(nTrialFrames,nGreenCells,nTrials);
    redF = nan(nTrialFrames,nRedCells,nTrials);
    for itrial = 1:nTrials
        ind = (trialStartInd(itrial)-nBaselineFr):...
            (trialStartInd(itrial)+nPostStimFr-1);
        greenF(:,:,itrial) = greenTC(ind,:);
        redF(:,:,itrial) = redTC(ind,:);
    end
    
    greenF0 = mean(greenF(1:nBaselineFr,:,:),1);
    redF0 = mean(redF(1:nBaselineFr,:,:),1);
    
    greenDFF{1} = (greenF - greenF0)./greenF0;
    redDFF{1} = (redF - redF0)./redF0;

    greenF = nan(nTrialFrames,nGreenCells,nTrials);
    redF = nan(nTrialFrames,nRedCells,nTrials);
    for itrial = 1:nTrials
        ind = (secondStimInd(itrial)-nBaselineFr):...
            (secondStimInd(itrial)+nPostStimFr-1);
        greenF(:,:,itrial) = greenTC(ind,:);
        redF(:,:,itrial) = redTC(ind,:);
    end
    
    greenDFF{2} = (greenF - greenF0)./greenF0;
    redDFF{2} = (redF - redF0)./redF0;
    
    save(fullfile(fn,folder,'adaptResponse'),'greenDFF','redDFF','secondStimDelay')
end
    