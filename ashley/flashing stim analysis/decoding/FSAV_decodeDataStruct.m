clear all
close all
ds = 'FSAV_attentionV1';
cellsOrDendrites = 1;
bxMice = true;
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

fnout = fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,'FSAV Choice');
%%
visualTrials = 1;
auditoryTrials = 2;
pressAlign = 1;
faAlign = 2;
crAlign = 3;
targetAlign = 4;
catchAlign = 5;

if bxMice
    visRT = [250 570];
    audRT = [150 450];
else
    visRT = [0 570];
    audRT = [0 570];
end
nFrMonitorDelay = 3; %using cLeverDown
nFrMonitorDelay_target = 2; %using cTargetOn
nFrVisDelay = 3; 
nFrRespwin = 3;
respCutoffDff = 0.002;
respAlpha = 0.01;
% frameRateHz = 30;


nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
cycTimeFr = mouse(1).expt(1).info.cycTimeFrames;

respwin = (nBaselineFr+nFrVisDelay+nFrMonitorDelay):...
    (nBaselineFr+nFrVisDelay+nFrMonitorDelay+nFrRespwin - 1);
basewin = (nBaselineFr-nFrRespwin+nFrMonitorDelay+1):...
    (nBaselineFr+nFrMonitorDelay);

respwin_target = (nBaselineFr+nFrVisDelay+nFrMonitorDelay_target):...
    (nBaselineFr+nFrVisDelay+nFrMonitorDelay_target+nFrRespwin - 1);
basewin_target = (nBaselineFr-nFrRespwin+nFrMonitorDelay_target+1):...
    (nBaselineFr+nFrMonitorDelay_target);
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
        
        oriFit = mouse(imouse).expt(iexp).oriTuning.oriFit;
        oriFit = oriFit(1:end-1,:);
        osi_fit = getOSIeaCell(oriFit,1:180);

        disp([ms '-' dt])

        d = mouse(imouse).expt(iexp).av(visualTrials);
        da = mouse(imouse).expt(iexp).av(auditoryTrials);
        % first responders
        tc = d.align(pressAlign).respTC;
        firstStimTrialInd = d.align(pressAlign).nCycles >= 2;
        firstStimResp = squeeze(mean(tc(respwin,:,firstStimTrialInd),1));
        firstStimBL = squeeze(mean(tc(basewin,:,firstStimTrialInd),1));
        minRespCells = (mean(firstStimResp,2) - mean(firstStimBL,2)) > respCutoffDff;
        firstStimResp_av = cat(2,firstStimResp,...
            squeeze(mean(da.align(pressAlign).respTC(...
            respwin,:,da.align(pressAlign).nCycles >= 2),1)));
        firstStimBL_av = cat(2,firstStimBL,...
            squeeze(mean(da.align(pressAlign).respTC(basewin,:,...
            da.align(pressAlign).nCycles >= 2),1)));
%         firstBaseResponsiveCells = ttest(firstStimResp,firstStimBL,'dim',2,...
%             'tail','right','alpha',respAlpha) & minRespCells;
        firstBaseResponsiveCells = ttest(firstStimResp_av,firstStimBL_av,'dim',2,...
            'tail','right','alpha',respAlpha) & minRespCells;
        if bxMice
            % false alarms
            tc = d.align(faAlign).respTC;
            rt = d.align(faAlign).reactTime;
            faTrialInd = rt > visRT(1) & rt < visRT(2);
            faResp = squeeze(mean(tc(respwin,:,faTrialInd),1));
            faBL = squeeze(mean(tc(basewin,:,faTrialInd),1));
            faCycN = d.align(faAlign).nCycles(faTrialInd);
            faOut = repmat({'fa'},[1,sum(faTrialInd)]);
            faOri = zeros(1,sum(faTrialInd));
        else
            faResp = [];
            faBL = [];
            faCycN = [];
            faOut = [];
            faOri = [];
        end

        % correct rejects
        tc = d.align(crAlign).respTC;
        crResp = squeeze(mean(tc(respwin,:,:),1));
        crBL = squeeze(mean(tc(basewin,:,:),1));
        crCycN = d.align(crAlign).nCycles;
        crOut = repmat({'cr'},[1,length(crCycN)]);
        crOri = zeros(1,length(crCycN));
        
        % hits and misses
        tc = d.align(targetAlign).respTC;
        rt = d.align(targetAlign).reactTime;
        targetTrialInd = rt > visRT(1);
        targetResp = squeeze(mean(tc(respwin_target,:,targetTrialInd),1));
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
        targetAlpha = respAlpha/(length(oris)-1);
        eaTargetResponsiveCells = cellfun(@(x,y) ttest(x,y,'dim',2,'tail','right',...
            'alpha',targetAlpha),targetRespOri,targetBLOri,'unif',0);
        allTargetResponsiveCells = ttest(targetResp,targetBL,'dim',2,'tail','right',...
            'alpha',respAlpha);    
        targetResponsiveCells = (sum(cell2mat(eaTargetResponsiveCells),2) > 0 | ...
            allTargetResponsiveCells) & minRespCells;
        targetCycN = d.align(targetAlign).nCycles(targetTrialInd);
        targetOut = d.align(targetAlign).outcome(targetTrialInd);
        targetOut(strcmp(targetOut,'success')) = {'h'};
        targetOut(strcmp(targetOut,'ignore')) = {'m'};
        
        if size(d.align,2) == catchAlign
            % catch trials - hits & misses
            tc = d.align(catchAlign).respTC;
            rt = d.align(catchAlign).reactTime;
            catchOut = d.align(catchAlign).outcome;
            catchTrialInd = rt > audRT(1) & cellfun(@(x) ~isempty(x),catchOut);
            if sum(catchTrialInd) > 2
                catchResp = squeeze(mean(tc(respwin_target,:,catchTrialInd),1));
                catchBL = squeeze(mean(tc(basewin_target,:,catchTrialInd),1));
                catchOri = d.align(catchAlign).ori(catchTrialInd);
                catchOut = d.align(catchAlign).outcome(catchTrialInd);
                catchOut(strcmp(catchOut,'FA')) = {'h'};
                catchOut(strcmp(catchOut,'CR')) = {'m'};
                % catch trials - correct reject
                tc = d.align(catchAlign).crRespTC;
                crCatchResp = squeeze(mean(tc(respwin,:,:),1));
                crCatchBL = squeeze(mean(tc(basewin,:,:),1));
                crCatchCycN = d.align(catchAlign).crNCycles;
                crCatchOut = repmat({'cr'},[1,length(crCatchCycN)]);
                crCatchOri = zeros(1,length(crCatchCycN));
            else
                catchResp = [];
                catchBL = [];
                catchOri = [];
                catchOut = [];
                crCatchResp = [];
                crCatchBL = [];
                crCatchOut = [];
                crCatchOri = [];
            end
        end

        dcExpt(t).mouse = ms;
        dcExpt(t).date = dt;
        dcExpt(t).info = 'cells x trials';
        dcExpt(t).attention.task = mouse(imouse).expt(iexp).attnTask;
        dcExpt(t).attention.effect = mouse(imouse).expt(iexp).hasAttn;
%         dcExpt(t).attention.rewarded = mouse(imouse).expt(iexp).catchRew;
        dcExpt(t).oriPref = tuning;
        dcExpt(t).oriFitTheta90 = tuningReliability;
        dcExpt(t).osi_fit = osi_fit;
        dcExpt(t).stimResp = cat(2,(faResp-faBL),(crResp-crBL),(targetResp - targetBL));
        dcExpt(t).trialOutcome = cat(2,faOut,crOut,targetOut);
        dcExpt(t).trialOrientation = cat(2,faOri,crOri,tOri);
        dcExpt(t).signifResponsiveCells = firstBaseResponsiveCells | targetResponsiveCells;
        dcExpt(t).firstBaseResponsiveCells = firstBaseResponsiveCells;
        dcExpt(t).targetResponsiveCells = targetResponsiveCells;
        dcExpt(t).firstStimResp = mean(firstStimResp,2);
        

        d = mouse(imouse).expt(iexp).av(auditoryTrials);
        % auditory false alarms
        if bxMice
            tc = d.align(faAlign).respTC;
            rt = d.align(faAlign).reactTime;
            faTrialInd = rt > audRT(1) & rt < audRT(2);
            faResp = squeeze(mean(tc(respwin,:,faTrialInd),1));
            faBL = squeeze(mean(tc(basewin,:,faTrialInd),1));
            faCycN = d.align(faAlign).nCycles(faTrialInd);
            faOut = repmat({'fa'},[1,sum(faTrialInd)]);
            faAmp = zeros(1,sum(faTrialInd));
        else
            faResp = [];
            faBL = [];
            faCycN = [];
            faOut = [];
            faAmp = [];
        end
        % auditory correct rejects
        tc = d.align(crAlign).respTC;
        crResp = squeeze(mean(tc(respwin,:,:),1));
        crBL = squeeze(mean(tc(basewin,:,:),1));
        crCycN = d.align(crAlign).nCycles;
        crOut = repmat({'cr'},[1,length(crCycN)]);
        crAmp = zeros(1,length(crCycN));
        % auditory hits and misses
        tc = d.align(targetAlign).respTC;
        rt = d.align(targetAlign).reactTime;
        targetTrialInd = rt > audRT(1);
        targetResp = squeeze(mean(tc(respwin_target,:,targetTrialInd),1));
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
        eaTargetResponsiveCells = cellfun(@(x,y) ttest(x,y,'dim',2,'tail','right',...
            'alpha',targetAlpha),targetRespAmp,targetBLAmp,'unif',0);
        allTargetResponsiveCells = ttest(targetResp,targetBL,'dim',2,'tail','right',...
            'alpha',respAlpha);    
        targetResponsiveCells = (sum(cell2mat(eaTargetResponsiveCells),2) > 0 | ...
            allTargetResponsiveCells) & minRespCells;
        targetCycN = d.align(targetAlign).nCycles(targetTrialInd);
        targetOut = d.align(targetAlign).outcome(targetTrialInd);
        targetOut(strcmp(targetOut,'success')) = {'h'};
        targetOut(strcmp(targetOut,'ignore')) = {'m'};

        dcExpt(t).audStimResp = cat(2,(faResp-faBL),(crResp-crBL),(targetResp - targetBL));
        dcExpt(t).audTrialOutcome = cat(2,faOut,crOut,targetOut);
        dcExpt(t).audTrialAmp = cat(2,faAmp,crAmp,tAmp);
        
        
        if size(mouse(imouse).expt(iexp).av(visualTrials).align,2) == catchAlign
            dcExpt(t).catchResp = cat(2,(faResp-faBL),(crCatchResp-crCatchBL),(catchResp-catchBL));
            dcExpt(t).catchOutcome = cat(2,faOut,crCatchOut,catchOut);
            dcExpt(t).catchOrientation = cat(2,zeros(1,length(faAmp)),crCatchOri,catchOri);
        end
        
    end
end
if exist(fnout,'dir') == 1
    save(fullfile(fnout,'FSAV_decodeData'),'dcExpt')
else
    mkdir(fnout)
    save(fullfile(fnout,'FSAV_decodeData'),'dcExpt')
end