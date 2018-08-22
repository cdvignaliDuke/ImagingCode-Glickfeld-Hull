clear all
close all
ds = 'FSAV_V1_decode';
cellsOrDendrites = 1;
%%
rc = behavConstsAV;
dataGroup = ds;
eval(dataGroup)

titleStr = ds;
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];

if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));    
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_trOutcomeStruct_dendrites' ds(5:end) '.mat']));    
end        
if cellsOrDendrites == 1
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_trOut_']); 
elseif cellsOrDendrites == 2
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr 'trOut_den_']); 
end

fnout = fullfile(rc.ashleyAnalysis,'FSAV Choice');
%%
visualTrials = 1;
auditoryTrials = 2;
pressAlign = 1;
faAlign = 2;
crAlign = 3;
targetAlign = 4;

visRT = [150 570];
audRT = [150 450];

nFrMonitorDelay = 3; %using cLeverDown
nFrMonitorDelay_target = 2; %using cTargetOn
optimalVisDelayInd = 10;
nFrVisDelay = -6:10; 
nRespwin = length(nFrVisDelay);
nFrRespwin = 3;
respCutoffDff = 0.01;
respAlpha = 0.01;
% frameRateHz = 30;


nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
cycTimeFr = mouse(1).expt(1).info.cycTimeFrames;
frameRate = 30;

respwin = cell(1,nRespwin);
respwin_target = cell(1,nRespwin);
for i = 1:nRespwin
    respwin{i} = (nBaselineFr+nFrVisDelay(i)+nFrMonitorDelay):...
        (nBaselineFr+nFrVisDelay(i)+nFrMonitorDelay+nFrRespwin - 1);
    respwin_target{i} = (nBaselineFr+nFrVisDelay(i)+nFrMonitorDelay_target):...
        (nBaselineFr+nFrVisDelay(i)+nFrMonitorDelay_target+nFrRespwin - 1);
end

basewin = (nBaselineFr-nFrRespwin+nFrMonitorDelay):...
    (nBaselineFr+nFrMonitorDelay-1);
basewin_target = (nBaselineFr-nFrRespwin+nFrMonitorDelay_target):...
    (nBaselineFr+nFrMonitorDelay_target-1);
%%
dcExpt = struct;

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if imouse == 1 && iexp == 1
            t = 1;
        else 
            t = t + 1;
        end

    ms = mouse(imouse).mouse_name;
    dt = mouse(imouse).expt(iexp).date;
    
    [~,fitPeak] = max(mouse(imouse).expt(iexp).oriTuning.oriFit,[],1);
    tuning = fitPeak-1;
    tuningReliability = mouse(imouse).expt(iexp).oriTuning.oriFitReliability;

    disp([ms '-' dt])
    dcExpt(t).mouse = ms;
    dcExpt(t).date = dt;
    dcExpt(t).visRTwin = visRT;
    dcExpt(t).audRTwin = audRT;

        for rw = 1:nRespwin

            dcExpt(t).respwin(rw).winMs = (([respwin{rw}(1) respwin{rw}(end)]...
                -(nBaselineFr))./frameRate).*1000;

            d = mouse(imouse).expt(iexp).av(visualTrials);

            tc = d.align(pressAlign).respTC;
            firstStimTrialInd = d.align(pressAlign).nCycles >= 2;
            firstStimResp = squeeze(mean(tc(respwin{rw},:,firstStimTrialInd),1));
            firstStimBL = squeeze(mean(tc(basewin,:,firstStimTrialInd),1));
            minRespCells = (mean(firstStimResp,2) - mean(firstStimBL,2)) > respCutoffDff;
            if rw == optimalVisDelayInd
                firstBaseResponsiveCells = ttest(firstStimResp,firstStimBL,'dim',2,...
                    'tail','right','alpha',respAlpha) & minRespCells;
            end

            tc = d.align(faAlign).respTC;
            rt = d.align(faAlign).reactTime;
            faTrialInd = rt > visRT(1) & rt < visRT(2);
            faResp = squeeze(mean(tc(respwin{rw},:,faTrialInd),1));
            faBL = squeeze(mean(tc(basewin,:,faTrialInd),1));
            faCycN = d.align(faAlign).nCycles(faTrialInd);
            faOut = repmat({'fa'},[1,sum(faTrialInd)]);
            faOri = zeros(1,sum(faTrialInd));
            faRT = rt(faTrialInd);

            tc = d.align(crAlign).respTC;
            crResp = squeeze(mean(tc(respwin{rw},:,:),1));
            crBL = squeeze(mean(tc(basewin,:,:),1));
            crCycN = d.align(crAlign).nCycles;
            crOut = repmat({'cr'},[1,length(crCycN)]);
            crOri = zeros(1,length(crCycN));
            crRT = nan(1,length(crCycN));

            tc = d.align(targetAlign).respTC;
            rt = d.align(targetAlign).reactTime;
            targetTrialInd = rt > visRT(1);
            targetResp = squeeze(mean(tc(respwin_target{rw},:,targetTrialInd),1));
            targetBL = squeeze(mean(tc(basewin_target,:,targetTrialInd),1));
            minRespCells = (mean(targetResp,2) - mean(targetBL,2)) > respCutoffDff;
            tOri = d.align(targetAlign).ori(targetTrialInd);
            oris = unique(tOri);
            targetRespOri = cell(1,length(oris));
            targetBLOri = cell(1,length(oris));
            for i = 1:length(oris)
                targetRespOri{i} = targetResp(:,tOri == oris(i));
                targetBLOri{i} = targetBL(:,tOri == oris(i));
            end
            if rw == optimalVisDelayInd
            targetAlpha = respAlpha/(length(oris)-1);
                eaTargetResponsiveCells = cellfun(@(x,y) ttest(x,y,'dim',2,'tail','right',...
                    'alpha',targetAlpha),targetRespOri,targetBLOri,'unif',0);
                allTargetResponsiveCells = ttest(targetResp,targetBL,'dim',2,'tail','right',...
                    'alpha',respAlpha);    
                targetResponsiveCells = (sum(cell2mat(eaTargetResponsiveCells),2) > 0 | ...
                    allTargetResponsiveCells) & minRespCells;
            end
            targetCycN = d.align(targetAlign).nCycles(targetTrialInd);
            targetOut = d.align(targetAlign).outcome(targetTrialInd);
            targetOut(strcmp(targetOut,'success')) = {'h'};
            targetOut(strcmp(targetOut,'ignore')) = {'m'};
            targetRT = rt(targetTrialInd);
            targetRT(strcmp(targetOut,'m')) = nan;

            dcExpt(t).respwin(rw).stimResp = cat(2,(faResp-faBL),(crResp-crBL),...
                (targetResp - targetBL));
            dcExpt(t).respwin(rw).trialOutcome = cat(2,faOut,crOut,targetOut);
            dcExpt(t).respwin(rw).trialOrientation = cat(2,faOri,crOri,tOri);
            if rw == optimalVisDelayInd
                dcExpt(t).info = 'cells x trials';
                dcExpt(t).oriPref = tuning;
                dcExpt(t).oriFitTheta90 = tuningReliability;
                dcExpt(t).signifResponsiveCells = firstBaseResponsiveCells | ...
                    targetResponsiveCells;
                dcExpt(t).visRT = cat(2,faRT,crRT,targetRT);
            end

            d = mouse(imouse).expt(iexp).av(auditoryTrials);

            tc = d.align(faAlign).respTC;
            rt = d.align(faAlign).reactTime;
            faTrialInd = rt > audRT(1) & rt < audRT(2);
            faResp = squeeze(mean(tc(respwin{rw},:,faTrialInd),1));
            faBL = squeeze(mean(tc(basewin,:,faTrialInd),1));
            faCycN = d.align(faAlign).nCycles(faTrialInd);
            faOut = repmat({'fa'},[1,sum(faTrialInd)]);
            faAmp = zeros(1,sum(faTrialInd));
            faRT = rt(faTrialInd);

            tc = d.align(crAlign).respTC;
            crResp = squeeze(mean(tc(respwin{rw},:,:),1));
            crBL = squeeze(mean(tc(basewin,:,:),1));
            crCycN = d.align(crAlign).nCycles;
            crOut = repmat({'cr'},[1,length(crCycN)]);
            crAmp = zeros(1,length(crCycN));
            crRT = nan(1,length(crCycN));

            tc = d.align(targetAlign).respTC;
            rt = d.align(targetAlign).reactTime;
            targetTrialInd = rt > audRT(1);
            targetResp = squeeze(mean(tc(respwin_target{rw},:,targetTrialInd),1));
            targetBL = squeeze(mean(tc(basewin_target,:,targetTrialInd),1));
            minRespCells = (mean(targetResp,2) - mean(targetBL,2)) > respCutoffDff;
            tAmp = d.align(targetAlign).amp(targetTrialInd);
            amps = unique(tAmp);
            targetRespAmp = cell(1,length(amps));
            targetBLAmp = cell(1,length(amps));
            for i = 1:length(amps)
                targetRespAmp{i} = targetResp(:,tAmp == amps(i));
                targetBLAmp{i} = targetBL(:,tAmp == amps(i));
            end
            if length(amps) > 1
                targetAlpha = respAlpha/(length(amps)-1);
            else
                targetAlpha = respAlpha;
            end
            if rw == optimalVisDelayInd
                eaTargetResponsiveCells = cellfun(@(x,y) ttest(x,y,'dim',2,'tail','right',...
                    'alpha',targetAlpha),targetRespAmp,targetBLAmp,'unif',0);
                allTargetResponsiveCells = ttest(targetResp,targetBL,'dim',2,'tail','right',...
                    'alpha',respAlpha);    
                targetResponsiveCells = (sum(cell2mat(eaTargetResponsiveCells),2) > 0 | ...
                    allTargetResponsiveCells) & minRespCells;
            end
            targetCycN = d.align(targetAlign).nCycles(targetTrialInd);
            targetOut = d.align(targetAlign).outcome(targetTrialInd);
            targetOut(strcmp(targetOut,'success')) = {'h'};
            targetOut(strcmp(targetOut,'ignore')) = {'m'};
            targetRT = rt(targetTrialInd);
            targetRT(strcmp(targetOut,'m')) = nan;

            dcExpt(t).respwin(rw).audStimResp = cat(2,(faResp-faBL),(crResp-crBL),(targetResp - targetBL));
            dcExpt(t).respwin(rw).audTrialOutcome = cat(2,faOut,crOut,targetOut);
            dcExpt(t).respwin(rw).audTrialAmp = cat(2,faAmp,crAmp,tAmp);
            if rw == optimalVisDelayInd
                dcExpt(t).audRT = cat(2,faRT,crRT,targetRT);
            end
        end    
    end
end
save(fullfile(fnout,'FSAV_respwinData'),'dcExpt')