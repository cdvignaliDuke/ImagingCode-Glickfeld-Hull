clear all
close all
ds = 'adaptDREADDS_V1_SOM_control';
doRedOnly = 0;
motionThreshold = 0.05;
%%

rc = behavConstsAV;
eval(ds)

saline = 1;
cno = 2;
isNotLabeled = 1;
isLabeled = 2;

nBaselineFr = (params.nBaselineMs/1000)*params.frameRate;
nRespwinFr = 0.1*params.frameRate;
nPostStimFr = (params.nPostStimMs/1000)*params.frameRate;
nTrialFrames = nBaselineFr+nPostStimFr;

respwin = (1+nBaselineFr+params.nFramesVisDelay_FS):...
    (nBaselineFr+params.nFramesVisDelay_FS+nRespwinFr);
basewin = (nBaselineFr - nRespwinFr+1):nBaselineFr;

%%
ms = struct;
exptDone = [];

for iexp = 1:size(expt,2)
    mouse = expt(iexp).mouse;
    subnum = mouse;
    expDate = expt(iexp).date;
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging');
    disp([mouse ' ' expDate])

    matchExpt = find(strcmp({expt.matchExptDate},expDate));
    label = expt(iexp).label;
    drug1 = expt(iexp).drug;
    drug2 = expt(matchExpt).drug;
    
    if any(exptDone == iexp)
        continue
    end
    if iexp == 1
        exptN = 1;
    else
        exptN = exptN + 1;
    end
    exptDone = cat(1,exptDone,iexp,matchExpt);
    
    ms(exptN).name = subnum;
    ms(exptN).dates = {expDate; expt(matchExpt).date};
    
    %% adapt
    folder1 = expt(iexp).adapt{1};
    folder2 = expt(matchExpt).adapt{1};
    % load data
    switch drug1
        case 'saline'
            load(fullfile(fn,expDate,folder1,'adaptResponse.mat')); 
            
            delays = unique(secondStimDelay);
            nDelay = length(delays);            
            delayGreenTC = cell(1,2);
            delayRedTC = cell(1,2);
            for istim = 1:2
                tc = greenDFF{istim};
                nCells = size(tc,2);
                delayGreenTC{istim} = nan(nTrialFrames,nCells,nDelay);
                for idelay = 1:nDelay
                    ind = secondStimDelay == delays(idelay);
                    delayGreenTC{istim}(:,:,idelay) = mean(tc(:,:,ind),3);
                end
                tc = redDFF{istim};
                nCells = size(tc,2);
                delayRedTC{istim} = nan(nTrialFrames,nCells,nDelay);
                for idelay = 1:nDelay
                    ind = secondStimDelay == delays(idelay);
                    delayRedTC{istim}(:,:,idelay) = mean(tc(:,:,ind),3);
                end           
            end            
            delayGreenResp = cellfun(@(x)... 
                squeeze(mean(x(respwin,:,:),1) - mean(x(basewin,:,:),1)),...
                delayGreenTC,'unif',0);
            delayRedResp = cellfun(@(x)... 
                squeeze(mean(x(respwin,:,:),1) - mean(x(basewin,:,:),1)),...
                delayRedTC,'unif',0);
            
            ms(exptN).drug(saline).name = 'saline';
            ms(exptN).drug(saline).label(isNotLabeled).name = 'unlabeled';
            ms(exptN).drug(saline).label(isNotLabeled).tc = delayGreenTC;
            ms(exptN).drug(saline).label(isNotLabeled).resp = delayGreenResp;
            ms(exptN).drug(saline).label(isLabeled).name = label;
            ms(exptN).drug(saline).label(isLabeled).tc = delayRedTC;
            ms(exptN).drug(saline).label(isLabeled).resp = delayRedResp;
            
            load(fullfile(fn,expt(matchExpt).date,folder2,'adaptResponse.mat'));            
            
            delays = unique(secondStimDelay);
            nDelay = length(delays);            
            delayGreenTC = cell(1,2);
            delayRedTC = cell(1,2);
            for istim = 1:2
                tc = greenDFF{istim};
                nCells = size(tc,2);
                delayGreenTC{istim} = nan(nTrialFrames,nCells,nDelay);
                for idelay = 1:nDelay
                    ind = secondStimDelay == delays(idelay);
                    delayGreenTC{istim}(:,:,idelay) = mean(tc(:,:,ind),3);
                end
                tc = redDFF{istim};
                nCells = size(tc,2);
                delayRedTC{istim} = nan(nTrialFrames,nCells,nDelay);
                for idelay = 1:nDelay
                    ind = secondStimDelay == delays(idelay);
                    delayRedTC{istim}(:,:,idelay) = mean(tc(:,:,ind),3);
                end           
            end            
            delayGreenResp = cellfun(@(x)... 
                squeeze(mean(x(respwin,:,:),1) - mean(x(basewin,:,:),1)),...
                delayGreenTC,'unif',0);
            delayRedResp = cellfun(@(x)... 
                squeeze(mean(x(respwin,:,:),1) - mean(x(basewin,:,:),1)),...
                delayRedTC,'unif',0);
            
            ms(exptN).drug(cno).name = 'cno';
            ms(exptN).drug(cno).label(isNotLabeled).name = 'unlabeled';
            ms(exptN).drug(cno).label(isNotLabeled).tc = delayGreenTC;
            ms(exptN).drug(cno).label(isNotLabeled).resp = delayGreenResp;
            ms(exptN).drug(cno).label(isLabeled).name = label;
            ms(exptN).drug(cno).label(isLabeled).tc = delayRedTC;
            ms(exptN).drug(cno).label(isLabeled).resp = delayRedResp;
        case 'cno'
            load(fullfile(fn,expt(matchExpt).date,folder2,'adaptResponse.mat')); 
            delays = unique(secondStimDelay);
            nDelay = length(delays);            
            delayGreenTC = cell(1,2);
            delayRedTC = cell(1,2);
            for istim = 1:2
                tc = greenDFF{istim};
                nCells = size(tc,2);
                delayGreenTC{istim} = nan(nTrialFrames,nCells,nDelay);
                for idelay = 1:nDelay
                    ind = secondStimDelay == delays(idelay);
                    delayGreenTC{istim}(:,:,idelay) = mean(tc(:,:,ind),3);
                end
                tc = redDFF{istim};
                nCells = size(tc,2);
                delayRedTC{istim} = nan(nTrialFrames,nCells,nDelay);
                for idelay = 1:nDelay
                    ind = secondStimDelay == delays(idelay);
                    delayRedTC{istim}(:,:,idelay) = mean(tc(:,:,ind),3);
                end           
            end            
            delayGreenResp = cellfun(@(x)... 
                squeeze(mean(x(respwin,:,:),1) - mean(x(basewin,:,:),1)),...
                delayGreenTC,'unif',0);
            delayRedResp = cellfun(@(x)... 
                squeeze(mean(x(respwin,:,:),1) - mean(x(basewin,:,:),1)),...
                delayRedTC,'unif',0);
            
            ms(exptN).drug(saline).name = 'saline';
            ms(exptN).drug(saline).label(isNotLabeled).name = 'unlabeled';
            ms(exptN).drug(saline).label(isNotLabeled).tc = delayGreenTC;
            ms(exptN).drug(saline).label(isNotLabeled).resp = delayGreenResp;
            ms(exptN).drug(saline).label(isLabeled).name = label;
            ms(exptN).drug(saline).label(isLabeled).tc = delayRedTC;
            ms(exptN).drug(saline).label(isLabeled).resp = delayRedResp;
            
            
            load(fullfile(fn,expDate,folder1,'adaptResponse.mat')); 
            delays = unique(secondStimDelay);
            nDelay = length(delays);            
            delayGreenTC = cell(1,2);
            delayRedTC = cell(1,2);
            for istim = 1:2
                tc = greenDFF{istim};
                nCells = size(tc,2);
                delayGreenTC{istim} = nan(nTrialFrames,nCells,nDelay);
                for idelay = 1:nDelay
                    ind = secondStimDelay == delays(idelay);
                    delayGreenTC{istim}(:,:,idelay) = mean(tc(:,:,ind),3);
                end
                tc = redDFF{istim};
                nCells = size(tc,2);
                delayRedTC{istim} = nan(nTrialFrames,nCells,nDelay);
                for idelay = 1:nDelay
                    ind = secondStimDelay == delays(idelay);
                    delayRedTC{istim}(:,:,idelay) = mean(tc(:,:,ind),3);
                end           
            end            
            delayGreenResp = cellfun(@(x)... 
                squeeze(mean(x(respwin,:,:),1) - mean(x(basewin,:,:),1)),...
                delayGreenTC,'unif',0);
            delayRedResp = cellfun(@(x)... 
                squeeze(mean(x(respwin,:,:),1) - mean(x(basewin,:,:),1)),...
                delayRedTC,'unif',0);
            
            ms(exptN).drug(cno).name = 'cno';
            ms(exptN).drug(cno).label(isNotLabeled).name = 'unlabeled';
            ms(exptN).drug(cno).label(isNotLabeled).tc = delayGreenTC;
            ms(exptN).drug(cno).label(isNotLabeled).resp = delayGreenResp;
            ms(exptN).drug(cno).label(isLabeled).name = label;
            ms(exptN).drug(cno).label(isLabeled).tc = delayRedTC;
            ms(exptN).drug(cno).label(isLabeled).resp = delayRedResp;
    end
    
    %% size tuning
    
    %% receptive field
    
    %% direction tuning
end

%% save 
% fnout = fullfile(rc.ashleyAnalysis,'Expt Summaries',['awData_audMod' ds]);
% save(fullfile(fnout,['awData_audMod' ds '_CaSummary']),'ms');
