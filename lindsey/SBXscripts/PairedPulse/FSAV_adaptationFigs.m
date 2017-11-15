clear all
close all
ds = 'FSAV_V1_100ms';
cellsOrDendrites = 1;
%%
rc = behavConstsAV;
if strcmp(rc.name,'ashle') & isempty(ds)
    dataGroup = ['awFSAVdatasets' ds];
elseif strcmp(rc.name,'ashle') & strcmp(ds(1:3),'FSA')
    dataGroup = ds;
elseif strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
else
    dataGroup = [];
end
eval(dataGroup)
titleStr = ds;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
elseif strcmp(titleStr(1),'_')
    titleStr = titleStr(2:end);    
elseif strcmp(titleStr(1:4), 'FSAV')
    titleStr = titleStr(6:end);
end
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];

if isempty(ds)
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' ds '.mat']));
elseif strcmp(ds(1:3),'FSA')
    if cellsOrDendrites == 1
        load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary_cells' ds(5:end) '.mat']));    
    elseif cellsOrDendrites == 2
        load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary_dendrites' ds(5:end) '.mat']));    
    end        
else
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' ds '.mat']));
end
% load(fullfile(rc.caOutputDir,dataGroup,[titleStr '_' mouse_str '_modCells.mat']));
if cellsOrDendrites == 1
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_startAlign']); 
elseif cellsOrDendrites == 2
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_den_startAlign']); 
end
    

%%
outcomecolmat = {'k';'r'};
avcolmat = {'k','c'};
avmat = {'visual';'auditory'};
outcomemat = {'hit';'miss'};
pressAlign = 1;
targetAlign = 2;
visualTrials = 1;
auditoryTrials = 2;
prevTrialVisual = 0;
prevTrialAuditory = 1;
nav = 2;
hitTrials = 1;
missTrials = 2;
nout = 2;
cycLengthFr = mouse(1).expt(1).info.cyc_time;
frRateHz = expt(1).frame_rate;
onTimeFr = 0.1*frRateHz;
nMonitorDelayFr = 2;
nBaselineFr = mouse(1).expt(1).pre_event_frames;
timestamp_1cyc = (-nBaselineFr:cycLengthFr-1);
plotTimeRange_1cyc = [-5 chop(timestamp_1cyc(end),2)];
nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
exptDates = {expt.date};
exptMice = {expt.SubNum};
%% orientation tuning
tuningReliabilityThresh = 45;
tuningDownSampFactor = 10;
avgResponseEaOri = cell(1,nexp);
semResponseEaOri = cell(1,nexp);
fitReliability = cell(1,nexp);
vonMisesFitAllCells = cell(1,nexp);
oriTuningTCs = cell(1,nexp);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if imouse == 1 && iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        ms = mouse(imouse).expt(iexp).mouse_name;
        dt = mouse(imouse).expt(iexp).date;
        ind = strcmp(exptMice,ms) & strcmp(exptDates,dt);
        dataPath = fullfile(rc.ashleyAnalysis,expt(ind).mouse,...
            'two-photon imaging', dt, expt(ind).dirtuning);
        tuningData = load(fullfile(dataPath,'oriTuningAndFits.mat'));
        avgResponseEaOri{exptN} = tuningData.avgResponseEaOri;
        semResponseEaOri{exptN} = tuningData.semResponseEaOri;
        fitReliability{exptN} = tuningData.fitReliability;
        vonMisesFitAllCells{exptN} = squeeze(tuningData.vonMisesFitAllCells);
        load(fullfile(dataPath,'oriTuningTCs'))
        oriTuningTCs{exptN} = tuningTC;
    end
end
for iexp = 1:nexp
    noFit = sum(vonMisesFitAllCells{iexp}) == 0;
    thisExptFitReliability = fitReliability{iexp};
    thisExptFitReliability(noFit) = 90;
    fitReliability{iexp} = thisExptFitReliability;
end

isTuned = cellfun(@(x) x < tuningReliabilityThresh,fitReliability,'unif',0)';
%% responses
basewin = 32:34;
respwin = 36:38;
basewinTarget = 31:33;
respwinTarget = 35:37;
nCycles = 8;
cellGroupsAlpha = 0.05;
lateCycles = 5:nCycles;
lateWinFrames = cycLengthFr*length(lateCycles);
tcEaCellEaCycAllOutcome = cell(nCycles,nexp,nav);
trialTCEaCycEaCell = cell(nCycles,nexp,nav);
trialTCShuffleEaExpt = cell(nexp,nav);
nTrialsEaCycAllOutcome = zeros(nCycles,nexp,nav);
tcEaCellEaCycByOutcome = cell(nCycles,nexp,nav,nout);
nTrialsEaCycByOutcome = zeros(nCycles,nexp,nav,nout);
exptName = cell(1,nexp);
h = cell(1,nexp);
p = cell(nexp,nCycles);
responsiveCellsTtestEaCyc = cell(nexp,nCycles);
responsive2itiCellsTtestEaCyc = cell(nexp,nCycles);
suppressed2itiCellsTtestEaCyc = cell(nexp,nCycles);
suppressedCellsTtestEaCyc = cell(nexp,nCycles);
siEaCycEaExpt = cell(nCycles,nexp);
responsiveCellsV = cell(nexp,1);
responsiveCellsA = cell(nexp,1);
targetResponsiveCells = cell(nexp,1);
lateSIEaExpt = cell(nexp,1);
ttestAVEaExpt = cell(nCycles,nexp);
ttestLateAVEaExpt = cell(nexp,1);
respEarlyLateV = cell(2,nexp);
respEarlyLateA = cell(2,nexp);
longTrialTCAV = cell(1,nexp);
longTrialTCAV_sem = cell(1,nexp);
nTrialsEaExptEaCyc = zeros(nexp,nCycles);
nHits = zeros(nexp,nav,nCycles);
nMiss = zeros(nexp,nav,nCycles);
lateWinSIEaExpt = cell(nexp,1);
ttestLateWinAVEaExpt = cell(nexp,1);
cyc2toEndSIEaExpt = cell(nexp,1);
allCycSIEaExpt = cell(nexp,1);
ttestCyc2toEndAVEaExpt = cell(nexp,1);
ttestAllCycAVEaExpt = cell(nexp,1);
lateWinRespAVEaExpt = cell(nexp,nav);
lateResponsiveCellsAV = cell(nexp,1);
lateSuppressedCellsAV = cell(nexp,1);
lateResponsiveCellsA = cell(nexp,1);
lateSuppressedCellsA = cell(nexp,1);
lateTtestP = cell(nexp,3);
visTargetVsLastBaseAuroc = cell(nexp,1);
visTargetVsLastBaseUTest = cell(nexp,1);
visTargetTCsHM = cell(nexp,2);
visTargetTCsEaExpt = cell(nexp,1);
lateWinNsCorrAV = cell(nexp,nav);
visTargetVsFirstBaseAuroc{exptN} = cell(nexp,1);
visTargetVsFirstBaseUTest{exptN} = cell(nexp,1);
itiFEaExpt = cell(nexp,1);
lateWinFEaExpt = cell(nexp,1);
itiCorrEaExpt = cell(nexp,1);
firstStimSuppressedCellsV = cell(nexp,1);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end
        exptName{exptN} = [mouse(imouse).expt(iexp).mouse_name '-' ...
            mouse(imouse).expt(iexp).date];
%         dAV_cycle1 = cat(3,mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{1},...
%             mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{1});
        dV_cycle1 = mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{1};
        dA_cycle1 = mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{1};
        responsiveCellsV{exptN} = ttest(squeeze(mean(dV_cycle1(respwin,:,:),1)),...
            squeeze(mean(dV_cycle1(basewin,:,:),1)),'dim',2,'tail','right','alpha',cellGroupsAlpha)';
        firstStimSuppressedCellsV{exptN} = ttest(squeeze(mean(dV_cycle1(respwin,:,:),1)),...
            squeeze(mean(dV_cycle1(basewin,:,:),1)),'dim',2,'tail','left')';
        responsiveCellsA{exptN} = ttest(squeeze(mean(dA_cycle1(respwin,:,:),1)),...
            squeeze(mean(dA_cycle1(basewin,:,:),1)),'dim',2,'tail','right','alpha',cellGroupsAlpha)';
        targetRespTtest = cellfun(@(x) ttest(squeeze(mean(x(...
            respwinTarget,:,:),1)),squeeze(mean(x(basewinTarget,:,:),1)),'dim',2,'tail','right')...
            ,mouse(imouse).expt(iexp).align(targetAlign).av(visualTrials).outcome(...
            hitTrials).stimResp(2:end),'unif',0);
        nonEmptyTargets = cellfun(@(x) size(x,1),targetRespTtest) > 1;
        targetResponsiveCells{exptN} = (sum(cell2mat(targetRespTtest(nonEmptyTargets)),2) > 0)';
        visTargetTCsHM{exptN,hitTrials} = mean(mouse(imouse).expt(iexp)...
            .align(targetAlign).av(visualTrials).outcome(hitTrials).stimResp{3},3);
        visTargetTCsHM{exptN,missTrials} = mean(mouse(imouse).expt(iexp)...
            .align(targetAlign).av(visualTrials).outcome(missTrials).stimResp{3},3);
        
        responseAllTrials = cell(nCycles,nav);

        dAV_nCycles = cat(3,mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{nCycles},...
            mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{nCycles});
        longTrialTCAV{exptN} = mean(dAV_nCycles,3);
        longTrialTCAV_sem{exptN} = ste(dAV_nCycles,3);
        
        dV_nCycles = mouse(imouse).expt(iexp).align(pressAlign).av(...
            visualTrials).outcome(hitTrials).cmlvCycResp{nCycles};
        dA_nCycles = mouse(imouse).expt(iexp).align(pressAlign).av(...
            auditoryTrials).outcome(hitTrials).cmlvCycResp{nCycles};
        lateResponsiveCellsAV{exptN} = ttest(squeeze(mean(dV_nCycles(respwin(1):end,:,:),1)),...
            squeeze(mean(dV_nCycles(basewin,:,:),1)),'dim',2,'tail','right','alpha',cellGroupsAlpha)';
        lateSuppressedCellsAV{exptN} = ttest(squeeze(mean(dV_nCycles(respwin(1):end,:,:),1)),...
            squeeze(mean(dV_nCycles(basewin,:,:),1)),'dim',2,'tail','left','alpha',cellGroupsAlpha)';
        lateResponsiveCellsA{exptN} = ttest(squeeze(mean(dV_nCycles(respwin(1):end,:,:),1)),...
            squeeze(mean(dV_nCycles(basewin,:,:),1)),'dim',2,'tail','right','alpha',cellGroupsAlpha)';
        lateSuppressedCellsA{exptN} = ttest(squeeze(mean(dV_nCycles(respwin(1):end,:,:),1)),...
            squeeze(mean(dV_nCycles(basewin,:,:),1)),'dim',2,'tail','left','alpha',cellGroupsAlpha)';
        
        
        ntr = size(dAV_nCycles,3);
        randTrialsV = randsample(ntr,floor(ntr/2));
        randTrialsA = setdiff(1:ntr,randTrialsV);
        
        trialTCShuffleEaExpt{exptN,visualTrials} = mean(dAV_nCycles(:,:,randTrialsV),3);
        trialTCShuffleEaExpt{exptN,auditoryTrials} = mean(dAV_nCycles(:,:,randTrialsA),3);
        
        % target vs. last base roc
        visTargetTCs = mouse(imouse).expt(iexp).align(targetAlign).av(visualTrials).outcome(hitTrials).resp;
        visTargetResp =  squeeze(mean(visTargetTCs(respwinTarget,:,:),1)...
            - mean(visTargetTCs(basewinTarget,:,:),1));
        visLastBaseResp =  squeeze(mean(visTargetTCs(respwinTarget-cycLengthFr,:,:),1)...
            - mean(visTargetTCs(basewinTarget-cycLengthFr,:,:),1));
        
        visTargetTCsEaExpt(exptN) = {mean(visTargetTCs,3)};

        nCells = size(visTargetResp,1);
        

        r = zeros(nCells,1);
        ut = zeros(nCells,1);
        for icell = 1:nCells
                r(icell) = roc_gh(visLastBaseResp(icell,:),visTargetResp(icell,:));
                [~,ut(icell)] = ranksum(visLastBaseResp(icell,:),visTargetResp(icell,:));
        end

        visTargetVsLastBaseAuroc{exptN} = r;
        visTargetVsLastBaseUTest{exptN} = ut;
        
        % target vs. first base roc
        visTC = mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).resp;
        firstBaseResp = squeeze(mean(visTC(respwin,:,:),1));
        trNCycles = double(mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).tcyc);
        [ncells,ntrials] = size(firstBaseResp);
        targetResp = nan(ncells,ntrials);
        for itrial = 1:ntrials
            respwinThisTrial = respwin + ((trNCycles(itrial)+1)*cycLengthFr);
            basewinThisTrial = basewin + ((trNCycles(itrial)+1)*cycLengthFr);
            targetResp(:,itrial) = squeeze(mean(visTC(respwinThisTrial,:,itrial),1) - mean(visTC(basewinThisTrial,:,itrial),1));
        end
        
        r = zeros(ncells,1);
        ut = zeros(ncells,1);
        for icell = 1:ncells
                r(icell) = roc_gh(firstBaseResp(icell,:),targetResp(icell,:));
                [~,ut(icell)] = ranksum(firstBaseResp(icell,:),targetResp(icell,:));
        end
 
        visTargetVsFirstBaseAuroc{exptN} = r;
        visTargetVsFirstBaseUTest{exptN} = ut;
        
        f = mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).F;
        minLengthTrialsInd = trNCycles >=nCycles;
        fMinLengthTrials = f(:,:,minLengthTrialsInd);
        nCyclesTrFrames = (nCycles*cycLengthFr)+nBaselineFr;
        itiFEaExpt{exptN} = squeeze(mean(mean(fMinLengthTrials(1:nBaselineFr,:,:),3),1));
        lateWinFEaExpt{exptN} = squeeze(mean(mean(fMinLengthTrials(end-lateWinFrames:end,:,:),3),1));
        
        fAllTrials = permute(cat(3,mouse(imouse).expt(iexp).align(pressAlign).av(...
            visualTrials).outcome(hitTrials).F(1:nBaselineFr,:,:),...
            mouse(imouse).expt(iexp).align(pressAlign).av(...
            visualTrials).outcome(missTrials).F(1:nBaselineFr,:,:),...
            mouse(imouse).expt(iexp).align(pressAlign).av(...
            visualTrials).outcome(3).F(1:nBaselineFr,:,:)),[2,1,3]);
        fITI = reshape(fAllTrials,ncells,nBaselineFr*size(fAllTrials,3))';
        itiCorrEaExpt{exptN} = corr(fITI);
        
        
        
        
        for iav = 1:nav
            d = mouse(imouse).expt(iexp).align(pressAlign).av(iav);        
            for icyc = 1:nCycles
%                 tcThisCycAllOutcome = cat(3,d.outcome(hitTrials).cmlvCycResp{icyc},...
%                     d.outcome(missTrials).cmlvCycResp{icyc});
                tcThisCycAllOutcome = cat(3,d.outcome(hitTrials).cmlvCycResp{icyc});
                nHits(exptN,iav,icyc) = size(d.outcome(hitTrials).cmlvCycResp{icyc},3);
                nMiss(exptN,iav,icyc) = size(d.outcome(missTrials).cmlvCycResp{icyc},3);
%                 tcThisCycAllOutcome = d.outcome(hitTrials).cmlvCycResp{icyc};
                nTrialsEaExptEaCyc(exptN,icyc) = size(tcThisCycAllOutcome,3);
                extraFrames = cycLengthFr*(icyc-1);
                trialTCEaCycEaCell{icyc,exptN,iav} = mean(tcThisCycAllOutcome,3);
                for ioutcome = 1:nout
                    tcThisOutcome = d.outcome(ioutcome).cmlvCycResp{icyc};
                    nTrialsEaCycByOutcome(icyc,exptN,iav,ioutcome) = size(tcThisOutcome,3);
                    baselineByOutcome = mean(tcThisOutcome(basewin+extraFrames,:,:),1);
                    tcEaCellEaCycByOutcome{icyc,exptN,iav,ioutcome}...
                        = mean(bsxfun(@minus,tcThisOutcome(...
                        end-(cycLengthFr+nBaselineFr)+1:end,:,:),baselineByOutcome),3);
                end
                if iav == 1
                    if icyc == 1
                        base_resp_thisCyc = mean(tcThisCycAllOutcome(basewin,:,:),1);
                        base_resp_iti = base_resp_thisCyc;
                    else
                        base_resp_thisCyc = mean(tcThisCycAllOutcome(basewin+extraFrames,:,:),1);
                        base_resp_iti = mean(tcThisCycAllOutcome(basewin,:,:),1);
                    end
                    responsiveCellsTtestEaCyc{exptN,icyc} = zeros(1,size(tcThisCycAllOutcome,2));
                    responsive2itiCellsTtestEaCyc{exptN,icyc} = zeros(1,size(tcThisCycAllOutcome,2));
                    suppressed2itiCellsTtestEaCyc{exptN,icyc} = zeros(1,size(tcThisCycAllOutcome,2));
                    suppressedCellsTtestEaCyc{exptN,icyc} = zeros(1,size(tcThisCycAllOutcome,2));
                    hResp{exptN} = zeros(1,size(tcThisCycAllOutcome,2));
                    for iCell = 1:size(tcThisCycAllOutcome,2)
                        if icyc == 1
                            [hResp{exptN}(1,iCell), ~] = ttest(squeeze(mean(...
                                tcThisCycAllOutcome(respwin,iCell,:),1)),...
                                squeeze(base_resp_thisCyc(:,iCell,:)),'tail','right');
                            if hResp{exptN}(1,iCell) == 1
                                [max_val, max_ind{exptN}(1,iCell)] = ...
                                    max(diff(mean(...
                                    tcThisCycAllOutcome(26:end,iCell,:),3),1),[],1);
                                if max_ind{exptN}(1,iCell) < respwin(end)-25
                                    responsiveCellsTtestEaCyc{exptN,icyc}(1,iCell) = 1;
                                end
                            end
                        else
                            [responsiveCellsTtestEaCyc{exptN,icyc}(1,iCell), ~]...
                                = ttest2(squeeze(mean(...
                                tcThisCycAllOutcome(respwin,iCell,:),1)),...
                                squeeze(base_resp_thisCyc(:,iCell,:)),'tail','right');
                        end
                        [responsive2itiCellsTtestEaCyc{exptN,icyc}(1,iCell), ~]...
                            = ttest2(squeeze(mean(tcThisCycAllOutcome(respwin+extraFrames,iCell,:),1)),...
                            squeeze(base_resp_iti(:,iCell,:)),'tail','right');
                        [suppressed2itiCellsTtestEaCyc{exptN,icyc}(1,iCell), ~]...
                            = ttest2(squeeze(mean(tcThisCycAllOutcome(respwin+extraFrames,iCell,:),1)),...
                            squeeze(base_resp_iti(:,iCell,:)),'tail','left');
                        [suppressedCellsTtestEaCyc{exptN,icyc}(1,iCell), ~]...
                            = ttest2(squeeze(mean(tcThisCycAllOutcome(respwin+extraFrames,iCell,:),1)),...
                            squeeze(base_resp_thisCyc(:,iCell,:)),'tail','left');
                        
                    end
                end
%                 thisExptRespEaCyc = nan(nCycles,size(tcThisCycAllOutcome,2));
                baselineAllTrials = mean(tcThisCycAllOutcome(basewin+extraFrames,:,:),1);
                responseAllTrials{icyc,iav} = mean(tcThisCycAllOutcome(...
                    respwin+extraFrames,:,:),1) - baselineAllTrials;
%                 thisExptRespEaCyc(icyc,:) = squeeze(mean(responseAllTrials{icyc,iav} - ...
%                     baselineAllTrials,3));
                nTrialsEaCycAllOutcome(icyc,exptN,iav) = size(tcThisCycAllOutcome,3);
                tcEaCellEaCycAllOutcome{icyc,exptN,iav} = mean(bsxfun(@minus,tcThisCycAllOutcome(...
                    end-(cycLengthFr+nBaselineFr)+1:end,:,:),...
                    baselineAllTrials),3);
            end
        end
        longTrialTCNCyc = cell(1,nav);
        longTrialTCnCyc{visualTrials} = mouse(imouse).expt(iexp).align(...
            pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{nCycles};
        longTrialTCnCyc{auditoryTrials} = mouse(imouse).expt(iexp).align(...
            pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{nCycles};
        lateWinRespAV = cellfun(@(x) squeeze(mean(x(end-lateWinFrames:end,:,:),1))',...
            longTrialTCnCyc,'unif',0);
        lateWinRespAVEaExpt(exptN,:) = cellfun(@(x) mean(x,1),lateWinRespAV,'unif',0);
        lateWinSIEaExpt{exptN} = getSelectivityIndex(...
            lateWinRespAV{visualTrials},lateWinRespAV{auditoryTrials});
        lateWinNsCorrAV{exptN,visualTrials} = corrcoef(lateWinRespAV{visualTrials});
        lateWinNsCorrAV{exptN,auditoryTrials} = corrcoef(lateWinRespAV{auditoryTrials});
        
        responseAllTrials = cellfun(@(x)squeeze(x)',responseAllTrials,'unif',0);
        siEaCycEaExpt(:,exptN) = cellfun(@(x,y) getSelectivityIndex(x,y),...
            responseAllTrials(:,visualTrials),responseAllTrials(:,auditoryTrials),'unif',0);
        lateSIEaExpt(exptN) = {getSelectivityIndex(...
            cell2mat(responseAllTrials(lateCycles,visualTrials)),...
            cell2mat(responseAllTrials(lateCycles,auditoryTrials)))};
        cyc2toEndSIEaExpt(exptN) = {getSelectivityIndex(...
            cell2mat(responseAllTrials(2:end,visualTrials)),...
            cell2mat(responseAllTrials(2:end,auditoryTrials)))};
        allCycSIEaExpt(exptN) = {getSelectivityIndex(...
            cell2mat(responseAllTrials(:,visualTrials)),...
            cell2mat(responseAllTrials(:,auditoryTrials)))};
        ttestAVEaExpt(:,exptN) = cellfun(@(x,y) ttest2(x,y,'dim',1,'alpha',0.05/(nCycles-1)),...
            responseAllTrials(:,visualTrials),responseAllTrials(:,auditoryTrials),...
            'unif',0);
        lateRespV = cell2mat(responseAllTrials(lateCycles,visualTrials));
        lateRespA = cell2mat(responseAllTrials(lateCycles,auditoryTrials));
%         ttestLateAVEaExpt{exptN} = {ttest2(lateRespV,lateRespA,'dim',1,'alpha',0.05)};
        [ttestLateAVEaExpt{exptN},lateTtestP{exptN,1}] = ttest2(lateRespV,lateRespA,'dim',1,'alpha',0.05);
        [ttestCyc2toEndAVEaExpt{exptN},lateTtestP{exptN,2}] = ttest2(...
            cell2mat(responseAllTrials(2:end,visualTrials)),...
            cell2mat(responseAllTrials(2:end,auditoryTrials)),'dim',1,'alpha',0.05);
        [ttestAllCycAVEaExpt{exptN},lateTtestP{exptN,3}] = ttest2(...
            cell2mat(responseAllTrials(:,visualTrials)),...
            cell2mat(responseAllTrials(:,auditoryTrials)),'dim',1,'alpha',0.05);
        ttestLateWinAVEaExpt{exptN} = ttest2(...
            lateWinRespAV{visualTrials},lateWinRespAV{auditoryTrials},'dim',1,'alpha',0.05);
        respEarlyLateV(1,exptN) = {mean(responseAllTrials{1,visualTrials},1)};
        respEarlyLateA(1,exptN) = {mean(responseAllTrials{1,auditoryTrials},1)};
        respEarlyLateV(2,exptN) = {mean(cell2mat(responseAllTrials(lateCycles,visualTrials)),1)};
        respEarlyLateA(2,exptN) = {mean(cell2mat(responseAllTrials(lateCycles,auditoryTrials)),1)};
    end
end
%% expt info table
nCellEaExpt = cellfun(@length,responsiveCellsV);
nRespCellEaExpt = cellfun(@sum,responsiveCellsV);
nModRespCellEaExpt = cellfun(@(x,y) sum(x & y),ttestLateAVEaExpt,responsiveCellsV);
pctModCellEaExpt = chop((nModRespCellEaExpt./nRespCellEaExpt),2)*100;
nTrialsLong = sum(nHits(:,:,nCycles)+nMiss(:,:,nCycles),2);
exptTable = table(nCellEaExpt,nRespCellEaExpt,nModRespCellEaExpt,...
    pctModCellEaExpt,nTrialsLong,nHits(:,visualTrials,end),nMiss(:,visualTrials,end),...
    nHits(:,auditoryTrials,end),nMiss(:,auditoryTrials,end));
exptTable.Properties.VariableNames = {'nCellsEaExpt','nRespCellEaExpt',...
    'nModRespCellEaExpt','pctModRespCells','nTrialsLong','nHitsV','nMissV','nHitsA','nMissA'};
exptTable.Properties.RowNames = exptName';
exptTable.Properties.DimensionNames = {'expt', 'info'};
writetable(exptTable,[fnout '_exptInfoTable.xls'])
disp(exptTable)
%% previous trial sorted responses
prevVisTrialTCEaCycAV = cell(nCycles,nav,nexp);
prevAudTrialTCEaCycAV = cell(nCycles,nav,nexp);
prevVisTCEaCycAV = cell(nCycles,nav,nexp);
prevAudTCEaCycAV = cell(nCycles,nav,nexp);
respEarlyLatePrevVis = cell(2,nav,nexp);
respEarlyLatePrevAud = cell(2,nav,nexp);
ttestPrevAVAllCyc = cell(nexp,nav);
ttestPrevAVLateCyc = cell(nexp,nav);
ttestPrevAVFirstCyc = cell(nexp,nav);
siPrevAVAllCyc = cell(nexp,nav);
siPrevAVLateCyc = cell(nexp,nav);
siPrevAVFirstCyc = cell(nexp,nav);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if imouse == 1 && iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        d = mouse(imouse).expt(iexp).align(pressAlign);

        nExpCycles = length(d.av(1).outcome(1).cmlvCycResp);
        cmlvCycTCPrevVis = cell(nav,nExpCycles);
        cmlvCycTCPrevAud = cell(nav,nExpCycles);

        for iav = 1:nav
            dAV = d.av(iav).outcome(hitTrials);
            prevTrial = dAV.prevTrType;
            tCycles = dAV.tcyc;

            cmlvCycPrevTrial = cell(1,nExpCycles);
            for icyc = 1:nExpCycles
                thisCycCmlvCycInd = tCycles >= icyc;
                cmlvCycPrevTrial{icyc} = prevTrial(thisCycCmlvCycInd);
            end

            cmlvCycTCPrevVis(iav,:) = cellfun(@(x,y) x(:,:,y == prevTrialVisual),...
                dAV.cmlvCycResp,cmlvCycPrevTrial,'unif',0);

            cmlvCycTCPrevAud(iav,:) = cellfun(@(x,y) x(:,:,y == prevTrialAuditory),...
                dAV.cmlvCycResp,cmlvCycPrevTrial,'unif',0);
        end
        
        for icyc = 1:nCycles
            tcPrevVisThisCyc = cmlvCycTCPrevVis(:,icyc);
            tcPrevAudThisCyc = cmlvCycTCPrevAud(:,icyc);
            extraFrames = cycLengthFr*(icyc-1);
            
            prevVisTrialTCEaCycAV(icyc,:,exptN) = cellfun(...
                @(x) mean(x,3),tcPrevVisThisCyc,'unif',0);
            prevAudTrialTCEaCycAV(icyc,:,exptN) = cellfun(...
                @(x) mean(x,3),tcPrevAudThisCyc,'unif',0);
            
            for iav = 1:nav
                %prev vis
                tcThisCycAV = tcPrevVisThisCyc{iav,:};
                base_resp_thisCyc = mean(tcThisCycAV(basewin+extraFrames,:,:),1);
                prevVisTCEaCycAV{icyc,iav,exptN} = bsxfun(@minus,...
                    tcThisCycAV(end-(cycLengthFr+nBaselineFr)+1:end,:,:),...
                    base_resp_thisCyc);
                %prev aud
                tcThisCycAV = tcPrevAudThisCyc{iav,:};
                base_resp_thisCyc = mean(tcThisCycAV(basewin+extraFrames,:,:),1);
                prevAudTCEaCycAV{icyc,iav,exptN} = bsxfun(@minus,...
                    tcThisCycAV(end-(cycLengthFr+nBaselineFr)+1:end,:,:),...
                    base_resp_thisCyc);
            end    
        end
        respThisExptPrevVis = cellfun(@(x) squeeze(mean(x(respwin,:,:),1)),...
            prevVisTCEaCycAV(:,:,exptN),'unif',0);
        oneTrialCyc = cellfun(@(x) size(x,1),respThisExptPrevVis) == 1;
        if any(oneTrialCyc(:))
            respThisExptPrevVis(oneTrialCyc) = cellfun(...
                @(x) x',respThisExptPrevVis(oneTrialCyc),'unif',0);
        end
        respThisExptPrevAud = cellfun(@(x) squeeze(mean(x(respwin,:,:),1)),...
            prevAudTCEaCycAV(:,:,exptN),'unif',0);        
        oneTrialCyc = cellfun(@(x) size(x,1),respThisExptPrevAud) == 1;
        if any(oneTrialCyc(:))
            respThisExptPrevAud(oneTrialCyc) = cellfun(...
                @(x) x',respThisExptPrevAud(oneTrialCyc),'unif',0);
        end
                
        respAllTrialsPrevVis = cell(1,nav);
        respAllTrialsPrevAud = cell(1,nav);
        respLateTrialsPrevVis = cell(1,nav);
        respLateTrialsPrevAud = cell(1,nav);
        respFirstTrialsPrevVis = cell(1,nav);
        respFirstTrialsPrevAud = cell(1,nav);
        for iav  = 1:nav
            respEarlyLatePrevVis(2,iav,exptN) = {mean(cell2mat(...
                respThisExptPrevVis(3:end,iav)'),2)};
            respEarlyLatePrevAud(2,iav,exptN) = {mean(cell2mat(...
                respThisExptPrevAud(3:end,iav)'),2)};
            
            respAllTrialsPrevVis(iav) = {cell2mat(respThisExptPrevVis(:,iav)')};
            respAllTrialsPrevAud(iav) = {cell2mat(respThisExptPrevAud(:,iav)')};
            
            respLateTrialsPrevVis(iav) = {cell2mat(respThisExptPrevVis(3:end,iav)')};
            respLateTrialsPrevAud(iav) = {cell2mat(respThisExptPrevAud(3:end,iav)')};
            
            respFirstTrialsPrevVis(iav) = {cell2mat(respThisExptPrevVis(1,iav)')};
            respFirstTrialsPrevAud(iav) = {cell2mat(respThisExptPrevAud(1,iav)')};
        end
        
        respEarlyLatePrevVis(1,:,exptN) = cellfun(...
            @(x) mean(x,2),respThisExptPrevVis(1,:),'unif',0);
        respEarlyLatePrevAud(1,:,exptN) = cellfun(...
            @(x) mean(x,2),respThisExptPrevAud(1,:),'unif',0);
        
        ttestPrevAVAllCyc(exptN,visualTrials) = {ttest2(...
            respAllTrialsPrevVis{visualTrials},respAllTrialsPrevAud{visualTrials},...
            'dim',2)};
        ttestPrevAVAllCyc(exptN,auditoryTrials) = {ttest2(...
            respAllTrialsPrevVis{auditoryTrials},respAllTrialsPrevAud{auditoryTrials},...
            'dim',2)};
        siPrevAVAllCyc(exptN,visualTrials) = {getSelectivityIndex(...
            respAllTrialsPrevVis{visualTrials}',respAllTrialsPrevAud{visualTrials}')};
        siPrevAVAllCyc(exptN,auditoryTrials) = {getSelectivityIndex(...
            respAllTrialsPrevVis{auditoryTrials}',respAllTrialsPrevAud{auditoryTrials}')};
        
        ttestPrevAVLateCyc(exptN,visualTrials) = {ttest2(...
            respLateTrialsPrevVis{visualTrials},respLateTrialsPrevAud{visualTrials},...
            'dim',2)};
        ttestPrevAVLateCyc(exptN,auditoryTrials) = {ttest2(...
            respLateTrialsPrevVis{auditoryTrials},respLateTrialsPrevAud{auditoryTrials},...
            'dim',2)};
        siPrevAVLateCyc(exptN,visualTrials) = {getSelectivityIndex(...
            respLateTrialsPrevVis{visualTrials}',respLateTrialsPrevAud{visualTrials}')};
        siPrevAVLateCyc(exptN,auditoryTrials) = {getSelectivityIndex(...
            respLateTrialsPrevVis{auditoryTrials}',respLateTrialsPrevAud{auditoryTrials}')};
        
        ttestPrevAVFirstCyc(exptN,visualTrials) = {ttest2(...
            respFirstTrialsPrevVis{visualTrials},respFirstTrialsPrevAud{visualTrials},...
            'dim',2)};
        ttestPrevAVFirstCyc(exptN,auditoryTrials) = {ttest2(...
            respFirstTrialsPrevVis{auditoryTrials},respFirstTrialsPrevAud{auditoryTrials},...
            'dim',2)};
        siPrevAVFirstCyc(exptN,visualTrials) = {getSelectivityIndex(...
            respFirstTrialsPrevVis{visualTrials}',respFirstTrialsPrevAud{visualTrials}')};
        siPrevAVFirstCyc(exptN,auditoryTrials) = {getSelectivityIndex(...
            respFirstTrialsPrevVis{auditoryTrials}',respFirstTrialsPrevAud{auditoryTrials}')};
        
    end
end
%% figures
figure; 
suptitle('response to each vis stim, each experiment (n cells)')
for i = 1:nexp
    subplot(4,5,i);
    ind = find(responsiveCellsTtestEaCyc{i,1});
    nCells = length(ind);
    for icyc = 1:5
        plot(mean(tcEaCellEaCycAllOutcome{icyc,i,visualTrials}(25:end,ind),2)); 
        hold on; 
        title(num2str(length(ind)))
    end 
end
responseCutoffMet = cellfun(@(x) mean(x(respwin,:),1) > 0.002, ...
    tcEaCellEaCycAllOutcome(1,:,visualTrials),'unif',0)';
[tcEaCycResponsiveCells,avgResponseEaCycAV,semResponseEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    responsiveCellsV,visualTrials,respwin,1);
nCellsResp1Stim = size(tcEaCycResponsiveCells{1,1},2);

% lateResponseCells = combineIndAcrossCellArray(...
%     responsive2itiCellsTtestEaCyc(:,lateCycles));
% lateResponseCellInd = cellfun(@(x,y) x == 1 & y == 0,lateResponseCells,...
%     responsiveCellsV(:,1),'unif',0);
% [tcEaCycLateRespCells,avgLateRespEaCycAV,semLateRespEaCycAV] = ...
%     getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
%     lateResponseCellInd,visualTrials,respwin,0);
% nCellsResp6Stim = size(tcEaCycLateRespCells{1,1},2);

% lateSuppCells = combineIndAcrossCellArray(...
%     suppressed2itiCellsTtestEaCyc(:,lateCycles));
% lateSuppCellInd = cellfun(@(x,y) x == 1 & y == 0,lateSuppCells,...
%     responsiveCellsAV(:,1),'unif',0);
% [tcEaCycLateSuppCells,avgLateSuppEaCycAV,semLateSuppEaCycAV] = ...
%     getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
%     lateSuppCellInd,visualTrials,respwin,0);
% nCellsSupp6Stim = size(tcEaCycLateSuppCells{1,1},2);

lateResponseCellInd = cellfun(@(x,y) x == 1 & y == 0,lateResponsiveCellsAV,...
    responsiveCellsV(:,1),'unif',0);
[tcEaCycLateRespCells,avgLateRespEaCycAV,semLateRespEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    lateResponseCellInd,visualTrials,respwin,0);
nCellsResp6Stim = size(tcEaCycLateRespCells{1,1},2);

lateSuppCellInd = cellfun(@(x,y) x == 1 & y == 0,lateSuppressedCellsAV,...
    responsiveCellsV(:,1),'unif',0);
[tcEaCycLateSuppCells,avgLateSuppEaCycAV,semLateSuppEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    lateSuppCellInd,visualTrials,respwin,0);
nCellsSupp6Stim = size(tcEaCycLateSuppCells{1,1},2);
[tcEaCycTargetCells,avgTargetEaCycAV,semTargetEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    targetResponsiveCells,visualTrials,respwin,0);
nCellsTargetStim = size(tcEaCycTargetCells{1,1},2);

figure
suptitle('responsive cells, adaptation analysis window, visual vs auditory trials')
subplot(2,2,1)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycResponsiveCells{icyc,visualTrials},1))-33).*(1000/frRateHz),mean(tcEaCycResponsiveCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(['n = ' num2str(nCellsResp1Stim)])
xlabel('Time (ms)')
ylabel('dF/F')
vline((respwin([1 end])-33).*(1000/frRateHz))
hold on
vline((basewin([1 end])-33).*(1000/frRateHz),'k:')
figAxForm([])

avcolmat = {[0 0 0],[0 1 1]};
subplot(2,2,3)
for iav = 1:nav
    h = errorbar(1:nCycles, avgResponseEaCycAV(:,iav), semResponseEaCycAV(:,iav),'o');
    h.Color = avcolmat{iav};
    h.MarkerFaceColor = [1 1 1];
    hold on
end
figYAxis([],'Normalized dF/F',[0 0.018])
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
legend({'visual','auditory'})
figAxForm([])

normTCEaCycResponsiveCells = getNormTCByCondition(tcEaCycResponsiveCells,1,...
    visualTrials,respwin);

responsiveCellsPctChangeEaCyc = getPctChangeFromCyc1AV(...
    avgResponseEaCycAV,visualTrials);
lateRespCellsPctChangeEaCyc = getPctChangeFromCyc1AV(...
    avgLateRespEaCycAV,visualTrials);
lateSuppCellsPctChangeEaCyc = getPctChangeFromCyc1AV(...
    1+avgLateSuppEaCycAV,visualTrials);
targetCellsPctChangeEaCyc = getPctChangeFromCyc1AV(...
    1+avgTargetEaCycAV,visualTrials);

[avgNormTCEaCycAV,avgNormResponseEaCycAV,semNormResponseEaCycAV] = ...
    getNormCycRespAV(normTCEaCycResponsiveCells,visualTrials,respwin);


subplot(2,2,2)
for icyc = 1:nCycles
    plot(((26:size(avgNormTCEaCycAV,1))-33).*(1000/frRateHz),avgNormTCEaCycAV(26:end,icyc,visualTrials))
    hold on
end
title(['n = ' num2str(nCellsResp1Stim)])
xlabel('Time (ms)')
ylabel('Normalized dF/F')
vline((respwin([1 end])-33).*(1000/frRateHz))
hold on
vline((basewin([1 end])-33).*(1000/frRateHz),'k:')
figAxForm([])

subplot(2,2,4)
for iav = 1:nav
    h = errorbar(1:nCycles, avgNormResponseEaCycAV(:,iav), semNormResponseEaCycAV(:,iav),'o');
    h.Color = avcolmat{iav};
    h.MarkerFaceColor = [1 1 1];
    hold on
end
figYAxis([],'Normalized dF/F',[0 1.1])
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
legend({'visual','auditory'})
figAxForm([])


print([fnout 'adaptationAVResponsiveCells'],'-dpdf','-fillpage')
%% cell groups AV
setFigParams4Print('portrait')
responseLim = [-0.01 0.02];
normRespLim = [0 1.1];
cellType = {'respond to 1st stim';'responsive late in trial';...
    'suppressed late in trial';'respond to target'};
figure;

subplot(3,3,1)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycLateRespCells{icyc,visualTrials},1))-33).*(1000/frRateHz),...
        mean(tcEaCycLateRespCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(sprintf('%s (%s)',cellType{2}, num2str(nCellsResp6Stim)))

subplot(3,3,2)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycLateSuppCells{icyc,visualTrials},1))-33).*(1000/frRateHz),...
        mean(tcEaCycLateSuppCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(sprintf('%s (%s)',cellType{3}, num2str(nCellsSupp6Stim)))

subplot(3,3,3)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycTargetCells{icyc,visualTrials},1))-33).*(1000/frRateHz),...
        mean(tcEaCycTargetCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(sprintf('%s (%s)',cellType{4}, num2str(nCellsTargetStim)))

subplot(3,3,4)
for iav = 1:nav
    errorbar(1:nCycles, avgLateRespEaCycAV(:,iav), semLateRespEaCycAV(:,iav));
    hold on
end
subplot(3,3,5)
for iav = 1:nav
    errorbar(1:nCycles, avgLateSuppEaCycAV(:,iav), semLateSuppEaCycAV(:,iav));
    hold on
end
subplot(3,3,6)
for iav = 1:nav
    errorbar(1:nCycles, avgTargetEaCycAV(:,iav), semTargetEaCycAV(:,iav));
    hold on
end

for itype = 1:3
    subplot(3,3,itype)
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',responseLim)
    figAxForm([])
    
    subplot(3,3,itype+3)
    figXAxis([],'Stimulus number',[0 nCycles+1])
    figYAxis([],'dF/F',responseLim)
    figAxForm([])
    legend({'visual','auditory'})
    title(cellType{itype+1})
end

subplot(3,3,7)
normRespLim = [0 1.1];
for iav = 1:nav
    plot(1:nCycles, lateRespCellsPctChangeEaCyc(:,iav));
    hold on
end
figXAxis([],'Stimulus number',[0 nCycles+1])
figYAxis([],'% change',normRespLim)
figAxForm([])
legend({'visual','auditory'})
title(cellType{2})
subplot(3,3,8)
normRespLim = [0.99 1.01];
for iav = 1:nav
    plot(1:nCycles, lateSuppCellsPctChangeEaCyc(:,iav));
    hold on
end
figXAxis([],'Stimulus number',[0 nCycles+1])
figYAxis([],'% change',normRespLim)
figAxForm([])
legend({'visual','auditory'})
title(cellType{3})
print([fnout 'adaptationAVLateCellGroups'],'-dpdf','-fillpage')

%% cell group time-courses
avcolmat = {'k','c'};

allRespCellsEaExpt = cellfun(@(x,y,z) (x & y) | z,...
    responsiveCellsV,responseCutoffMet,lateResponseCellInd,'unif',0);

taskCellsEaExp = cellfun(@(x,y,z,a) (x & y) | z | a,...
    responsiveCellsV,responseCutoffMet,...
    lateResponseCellInd, lateSuppCellInd,'unif',0);
trialTCEaCycTaskCells = getCycResponse4CellIndAV(trialTCEaCycEaCell,...
    taskCellsEaExp,visualTrials,respwin,1);
trialTCEaCycResponsiveCells = getCycResponse4CellIndAV(trialTCEaCycEaCell,...
    responsiveCellsV,visualTrials,respwin,1);
trialTCEaCycLateRespCells = getCycResponse4CellIndAV(trialTCEaCycEaCell,...
    lateResponseCellInd,visualTrials,respwin,0);
trialTCEaCycLateSuppCells = getCycResponse4CellIndAV(trialTCEaCycEaCell,...
    lateSuppCellInd,visualTrials,respwin,0);
trialTCEaCycAllRespCells = getCycResponse4CellIndAV(trialTCEaCycEaCell,...
    allRespCellsEaExpt,visualTrials,respwin,0);

responseLim = [-0.04 0.065];
setFigParams4Print('portrait')
figure
subplot 321
for iav = 1:nav
    y = mean(trialTCEaCycResponsiveCells{6,iav},2);
    yerr = ste(trialTCEaCycResponsiveCells{6,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{1})
figAxForm([])
subplot 322
for iav = 1:nav
    nc = size(trialTCEaCycResponsiveCells{nCycles,iav},2);
    y = mean(trialTCEaCycResponsiveCells{nCycles,iav},2);
    yerr = ste(trialTCEaCycResponsiveCells{nCycles,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(sprintf('%s (%s)',cellType{1},num2str(nc)))
figAxForm([])

subplot 323
for iav = 1:nav
    y = mean(trialTCEaCycLateRespCells{6,iav},2);
    yerr = ste(trialTCEaCycLateRespCells{6,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{2})
figAxForm([])
subplot 324
for iav = 1:nav
    nc = size(trialTCEaCycLateRespCells{nCycles,iav},2);
    y = mean(trialTCEaCycLateRespCells{nCycles,iav},2);
    yerr = ste(trialTCEaCycLateRespCells{nCycles,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(sprintf('%s (%s)',cellType{2},num2str(nc)))
figAxForm([])

subplot 325
for iav = 1:nav
    y = mean(trialTCEaCycLateSuppCells{6,iav},2);
    yerr = ste(trialTCEaCycLateSuppCells{6,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{3})
figAxForm([])
subplot 326
for iav = 1:nav
    nc = size(trialTCEaCycLateSuppCells{nCycles,iav},2);
    y = mean(trialTCEaCycLateSuppCells{nCycles,iav},2);
    yerr = ste(trialTCEaCycLateSuppCells{nCycles,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(sprintf('%s (%s)',cellType{3},num2str(nc)))
figAxForm([])

for iplot = 1:6
    subplot(3,2,iplot)
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim)    
end
print([fnout 'tcAVLateCellGroups'],'-dpdf','-fillpage')

%%
tc = squeeze(trialTCEaCycEaCell(nCycles,:,:));

tcTaskCells = cell(1,nav);
tcTaskCellsShuffle = cell(1,nav);
for iav = 1:nav
    tcTaskCells(1,iav) = {cell2mat(cellfun(@(x,y) x(:,y),tc(:,iav),taskCellsEaExp,'unif',0)')};
    tcTaskCellsShuffle(1,iav) = {cell2mat(cellfun(@(x,y) x(:,y),trialTCShuffleEaExpt(:,iav),taskCellsEaExp,'unif',0)')};
end

n = size(tcTaskCells{1,1},2);

responseLim = [-0.01 0.02];
subLim =  [-0.005 0.01];
setFigParams4Print('portrait')
figure
subplot 311
for iav = 1:nav
    y = mean(tcTaskCells{iav},2);
    yerr = ste(tcTaskCells{iav},2);    
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav});
    hold on
end
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',responseLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title(sprintf('all anti-responsive cells (%s)',num2str(n)))

tcTaskCellsSub = tcTaskCells{visualTrials} - tcTaskCells{auditoryTrials};
tcTaskCellsShuffSub = tcTaskCellsShuffle{visualTrials} - tcTaskCells{auditoryTrials};

subplot 312
y = mean(tcTaskCellsSub,2);
yerr  = ste(tcTaskCellsSub,2);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcTaskCellsShuffSub,2);
yerr  = ste(tcTaskCellsShuffSub,2);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'b');
hold on
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',subLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
hold on
hline(0,'k--');
title('subtraction (black) and shufffled subtraction (blue)')

subplot 313
responseLim = [0 0.005];
for iav = 1:nav
    y = std(tcTaskCells{iav},[],2);   
    tt = ((26:length(y))-33).*(1000/frRateHz);
    plot(tt,y(26:length(y)),avcolmat{iav});
    hold on
end
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'standard deviation dF/F',responseLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title(sprintf('all anti-responsive cells (%s)',num2str(n)))

print([fnout '_tcAVallAntiCells'],'-dpdf','-fillpage')



%% baseline

[avgBaselineEaCycAV,semBaselineEaCycAV] = ...
    getCycBaseline4CellIndAV(trialTCEaCycEaCell,...
    responsiveCellsV,visualTrials,basewin,respwin,cycLengthFr,1);

[avgBaselineLateRespEaCycAV,semBaselineLateRespEaCycAV] = ...
    getCycBaseline4CellIndAV(trialTCEaCycEaCell,...
    lateResponseCellInd,visualTrials,basewin,respwin,cycLengthFr,0);

[avgBaselineLateSuppEaCycAV,semBaselineLateSuppEaCycAV] = ...
    getCycBaseline4CellIndAV(trialTCEaCycEaCell,...
    lateSuppCellInd,visualTrials,basewin,respwin,cycLengthFr,0);

[avgBaselineTargetCellsEaCycAV,semBaselineTargetCellsEaCycAV] = ...
    getCycBaseline4CellIndAV(trialTCEaCycEaCell,...
    targetResponsiveCells,visualTrials,basewin,respwin,cycLengthFr,1);

setFigParams4Print('landscape')
responseLim = [-0.025 0.04];
normRespLim = [0 1.1];
figure
suptitle('baseline each cycle')
for iav = 1:nav
    subplot(1,4,1)
    errorbar(1:nCycles, avgBaselineEaCycAV(:,iav), ...
        semBaselineEaCycAV(:,iav));
    hold on
    subplot(1,4,2)
    errorbar(1:nCycles, avgBaselineLateRespEaCycAV(:,iav), ...
        semBaselineLateRespEaCycAV(:,iav));
    hold on
    subplot(1,4,3)
    errorbar(1:nCycles, avgBaselineLateSuppEaCycAV(:,iav), ...
        semBaselineLateSuppEaCycAV(:,iav));
    hold on
    subplot(1,4,4)
    errorbar(1:nCycles, avgBaselineTargetCellsEaCycAV(:,iav), ...
        semBaselineTargetCellsEaCycAV(:,iav));
    hold on
end

for itype = 1:4
    subplot(1,4,itype)
    figXAxis([],'Stimulus number',[0 nCycles+1])
    figYAxis([],'dF/F',responseLim)
    figAxForm([])
    legend({'visual','auditory'},'location','northeastoutside')
    title(cellType{itype})
end
print([fnout 'baseAVCellGroups'],'-dpdf','-fillpage')

%% tuning
binEdgesOri = 0:1:181;
[~, tuningPeak] = cellfun(@(x) max(x),vonMisesFitAllCells,'unif',0);
figure;
ind = cellfun(@(x,y) x & logical(y),isTuned,responsiveCellsV,'unif',0);
tp = cell2mat(cellfun(@(x,y) x(y),tuningPeak',ind,'unif',0)'); 
tp_rad = deg2rad(tp) - (pi/2);
subplot(2,2,1)
histogram(tp,binEdgesOri)
title({sprintf('%s (%s/%s)',cellType{1},...
    num2str(sum(cellfun(@sum,ind))),num2str(nCellsResp1Stim))...
    ['avg= ' num2str(mean(tp)) ' (deg)']})
    
ind = cellfun(@(x,y) x & logical(y),isTuned,lateResponseCellInd,'unif',0);
tp = cell2mat(cellfun(@(x,y) x(y),tuningPeak',ind,'unif',0)');  
tp_rad = deg2rad(tp) - (pi/2);  
subplot(2,2,2)
histogram(tp,binEdgesOri)
title({sprintf('%s (%s/%s)',cellType{2},...
    num2str(sum(cellfun(@sum,ind))),num2str(nCellsResp6Stim))...
    ['avg= ' num2str(mean(tp)) ' (deg)']})

ind = cellfun(@(x,y) x & logical(y),isTuned,lateSuppCellInd,'unif',0);
tp = cell2mat(cellfun(@(x,y) x(y),tuningPeak',ind,'unif',0)');   
tp_rad = deg2rad(tp) - (pi/2); 
subplot(2,2,3)
histogram(tp,binEdgesOri)
title({sprintf('%s (%s/%s)',cellType{3},...
    num2str(sum(cellfun(@sum,ind))),num2str(nCellsSupp6Stim)),...
    ['avg= ' num2str(mean(tp)) ' (deg)']})

ind = cellfun(@(x,y) x & logical(y),isTuned,targetResponsiveCells,'unif',0);
tp = cell2mat(cellfun(@(x,y) x(y),tuningPeak',ind,'unif',0)');   
tp_rad = deg2rad(tp) - (pi/2); 
subplot(2,2,4)
histogram(tp,binEdgesOri)
title({sprintf('%s (%s/%s)',cellType{4},...
    num2str(sum(cellfun(@sum,ind))),num2str(sum(cellfun(@sum,targetResponsiveCells))))...
    ['avg= ' num2str(mean(tp)) ' (deg)']})
for i = 1:4
    subplot(2,2,i)
    figXAxis([],'tuning peak',[0 181])
    figYAxis([],'n cells',[])
end
%% second stim
secondStimRespInd = responsiveCellsTtestEaCyc(:,2);%cellfun(@(x,y) x == 0 & y == 1,responsiveCellsTtestEaCyc(:,1),responsiveCellsTtestEaCyc(:,2),'unif',0);

[tcEaCycSecondRespCells,avgSecondRespEaCycAV,semSecondRespEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    secondStimRespInd,visualTrials,respwin,1);
nCellsResp2Stim = size(tcEaCycSecondRespCells{1,1},2);

normTCEaCycSecondRespCells = getNormTCByCondition(tcEaCycSecondRespCells,1,...
    visualTrials,respwin);
[avgNormSecondRespTCEaCycAV,avgNormSecondRespEaCycAV,semNormSecondRespEaCycAV] = ...
    getNormCycRespAV(normTCEaCycSecondRespCells,visualTrials,respwin);

figure
suptitle('2nd stim responsive cells')
subplot(2,2,1)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycSecondRespCells{icyc,visualTrials},1))-33).*(1000/frRateHz),...
        mean(tcEaCycSecondRespCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(['n = ' num2str(nCellsResp2Stim)])
xlabel('Time (ms)')
ylabel('dF/F')
vline((respwin([1 end])-33).*(1000/frRateHz))

subplot(2,2,3)
for iav = 1:nav
    errorbar(1:nCycles, avgSecondRespEaCycAV(:,iav), semSecondRespEaCycAV(:,iav));
    hold on
end
xlim([0 nCycles+1])
xlabel('Stimulus number')
ylabel('dF/F')
legend({'visual','auditory'})

subplot(2,2,2)
for icyc = 1:nCycles
    plot(((26:size(avgNormSecondRespTCEaCycAV,1))-33).*(1000/frRateHz),...
        avgNormSecondRespTCEaCycAV(26:end,icyc,visualTrials))
    hold on
end
title(['n = ' num2str(nCellsResp2Stim)])
xlabel('Time (ms)')
ylabel('Normalized dF/F')
vline((respwin([1 end])-33).*(1000/frRateHz))

subplot(2,2,4)
for iav = 1:nav
    errorbar(1:nCycles, avgNormSecondRespEaCycAV(:,iav), ...
        semNormSecondRespEaCycAV(:,iav));
    hold on
end
ylim([0.8 1.1])
xlim([0 nCycles+1])
xlabel('Stimulus number')
ylabel('Normalized dF/F')
legend({'visual','auditory'})
print([fnout 'AVrespSecondStim'],'-dpdf','-fillpage')
%% cell groups
BLRespCells = logical(cell2mat(responsiveCellsV')) & cell2mat(responseCutoffMet');
TarRespCells = logical(cell2mat(targetResponsiveCells'));
lateBLRespCells = cell2mat(lateResponseCellInd');
lateBLSuppCells = cell2mat(lateSuppCellInd');
anticipationCells = BLRespCells | lateBLRespCells | lateBLSuppCells;

%% cell group scatters
lateWinRespV = cell2mat(lateWinRespAVEaExpt(:,visualTrials)');
lateWinRespA = cell2mat(lateWinRespAVEaExpt(:,auditoryTrials)');

respLim = [-0.17 0.3];

figure;
subplot 221
ind = lateBLRespCells | BLRespCells;
scatter(lateWinRespV(ind),lateWinRespA(ind),'ko');
hold on
plot(respLim,respLim,'k--')
hold on
[~,p] = ttest(lateWinRespV(ind),lateWinRespA(ind));
figXAxis([],'visual (dF/F)',respLim)
figYAxis([],'auditory (dF/F)',respLim)
figAxForm([])
title(sprintf('all responsive cells (%s), p = %s',num2str(sum(ind)),...
    num2str(p)))

subplot 222
scatter(lateWinRespV(BLRespCells),lateWinRespA(BLRespCells),'ko');
hold on
plot(respLim,respLim,'k--')
hold on
[~,p] = ttest(lateWinRespV(BLRespCells),lateWinRespA(BLRespCells));
figXAxis([],'visual (dF/F)',respLim)
figYAxis([],'auditory (dF/F)',respLim)
figAxForm([])
title(sprintf('first stim responsive (%s), p = %s',num2str(sum(BLRespCells)),...
    num2str(p)))

subplot 223
scatter(lateWinRespV(lateBLRespCells),lateWinRespA(lateBLRespCells),'ko');
hold on
plot(respLim,respLim,'k--')
hold on
[~,p] = ttest(lateWinRespV(lateBLRespCells),lateWinRespA(lateBLRespCells));
figXAxis([],'visual (dF/F)',respLim)
figYAxis([],'auditory (dF/F)',respLim)
figAxForm([])
title(sprintf('late responsive (%s), p = %s',num2str(sum(lateBLRespCells)),...
    num2str(p)))

subplot 224
scatter(lateWinRespV(lateBLSuppCells),lateWinRespA(lateBLSuppCells),'ko');
hold on
plot(respLim,respLim,'k--')
hold on
[~,p] = ttest(lateWinRespV(lateBLSuppCells),lateWinRespA(lateBLSuppCells));
figXAxis([],'visual (dF/F)',respLim)
figYAxis([],'auditory (dF/F)',respLim)
figAxForm([])
title(sprintf('late suppressed (%s), p = %s',num2str(sum(lateBLSuppCells)),...
    num2str(p)))
print([fnout '_scatterCellGroups'],'-dpdf','-fillpage')

%% selectivity 
% respCells = BLRespCells | lateBLRespCells |lateBLSuppCells;
% respCells = BLRespCells | lateBLRespCells;
respCells = BLRespCells;

earlyModAV = cell2mat(combineIndAcrossCellArray(ttestAVEaExpt(1:2,:)')');
lateModAV = cell2mat(ttestLateAVEaExpt');

lateSI = cell2mat(lateSIEaExpt');

binEdgesSI = -10:0.5:10;

siLim = [-10 10];
figure;
h = histogram(lateSI(respCells),binEdgesSI);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [1 1 1];
h = vline(0,'k-');
h.LineWidth = 2;
h = vline(mean(lateSI(respCells)),'r-');
h.LineWidth = 2;
hold on
h = histogram(lateSI(lateModAV & respCells),binEdgesSI);
h.FaceColor = 'k';
h.EdgeColor = [1 1 1];
figXAxis([],'late trial selectivity index (using cycles 5-8)',siLim)
figYAxis([],'n ncells',[])
figAxForm([])
title({sprintf('significantly modulated late (%s/%s)'...
    ,num2str(sum(lateModAV & respCells)),num2str(sum(respCells))),...
    'red: avg SI all cels, pink: avg SI signif cells'})
h = vline(mean(lateSI(respCells & lateModAV)),'m-');
h.LineWidth = 2;
print([fnout 'SIofModCells'],'-dpdf','-fillpage')

% response to each cyc, each modulated group
respV = cell2mat(cellfun(@(x) mean(x(respwin,:),1),...
    tcEaCellEaCycAllOutcome(:,:,visualTrials),'unif',0))';
respA = cell2mat(cellfun(@(x) mean(x(respwin,:),1),...
    tcEaCellEaCycAllOutcome(:,:,auditoryTrials),'unif',0))';

respLim = [0 0.015];
figure
suptitle('among responsive cells...')
subplot 121
ind = respCells & lateModAV & lateSI > 0;
y = mean(respV(ind,:),1);
yerr = ste(respV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateModAV & lateSI > 0;
y = mean(respA(ind,:),1);
yerr = ste(respA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells, V > A (%s)', num2str(sum(ind))))
legend({'visual','auditory'},'location','southwest')

subplot 122
ind = respCells & lateModAV & lateSI < 0;
y = mean(respV(ind,:),1);
yerr = ste(respV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateModAV & lateSI < 0;
y = mean(respA(ind,:),1);
yerr = ste(respA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells, V < A (%s)', num2str(sum(ind))))
print([fnout 'cycRespModCells'],'-dpdf','-fillpage')

%% adaptatoin across +/- SI groups

respLim = [0 0.035];
figure
suptitle('among responsive cells...')
subplot 221
ind = BLRespCells & lateSI > 0;
y = mean(respV(ind,:),1);
yerr = ste(respV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateSI > 0;
y = mean(respA(ind,:),1);
yerr = ste(respA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells (cyc 5-8), V > A (%s)', num2str(sum(ind))))
legend({'visual','auditory'},'location','northeast')

subplot 222
ind = BLRespCells & lateSI < 0;
y = mean(respV(ind,:),1);
yerr = ste(respV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateSI < 0;
y = mean(respA(ind,:),1);
yerr = ste(respA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells (cyc 5-8), V < A (%s)', num2str(sum(ind))))
legend({'visual','auditory'},'location','northeast')

normRespV = respV./respV(:,1);
normRespA = respA./respA(:,1);

respLim = [-0.1 1.1];
subplot 223
ind = BLRespCells & lateSI > 0;
y = mean(normRespV(ind,:),1);
yerr = ste(normRespV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateSI > 0;
y = mean(normRespA(ind,:),1);
yerr = ste(normRespA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'Normalized dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells (cyc 5-8), V > A (%s)', num2str(sum(ind))))
legend({'visual','auditory'},'location','northeast')

subplot 224
ind = BLRespCells & lateSI < 0;
y = mean(normRespV(ind,:),1);
yerr = ste(normRespV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateSI < 0;
y = mean(normRespA(ind,:),1);
yerr = ste(normRespA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'Normalized dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells (cyc 5-8), V < A (%s)', num2str(sum(ind))))
legend({'visual','auditory'},'location','northeast')

print([fnout 'cycRespSIGroupsRespCells'],'-dpdf','-fillpage')


earlyLateV = cell2mat(respEarlyLateV)';
earlyLateA = cell2mat(respEarlyLateA)';

earlyLateVNorm = earlyLateV./earlyLateV(:,1);
earlyLateANorm = earlyLateA./earlyLateA(:,1);

figure;
suptitle('1st stim responsive cells')
subplot 121
ind = BLRespCells & lateSI > 0;
h = bar(1,mean(earlyLateVNorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{1};
hold on
errorbar(1,mean(earlyLateVNorm(ind,2)),ste(earlyLateVNorm(ind,2),1),'k');
h = bar(2,mean(earlyLateANorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{2};
hold on
errorbar(2,mean(earlyLateANorm(ind,2)),ste(earlyLateANorm(ind,2),1),'c');
figXAxis([],'trial type',[0 3],1:2,{'visual';'auditory'})
figYAxis([],'% adapt',respLim)
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title(sprintf('V > A (%s)',num2str(sum(ind))))
subplot 122
ind = BLRespCells & lateSI < 0;
h = bar(1,mean(earlyLateVNorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{1};
hold on
errorbar(1,mean(earlyLateVNorm(ind,2)),ste(earlyLateVNorm(ind,2),1),'k');
h = bar(2,mean(earlyLateANorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{2};
hold on
errorbar(2,mean(earlyLateANorm(ind,2)),ste(earlyLateANorm(ind,2),1),'c');
figXAxis([],'trial type',[0 3],1:2,{'visual';'auditory'})
figYAxis([],'% adapt',respLim)
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title(sprintf('V < A modulated cells (%s)',num2str(sum(ind))))
print([fnout 'SIGroupsRespCellsNormBar'],'-dpdf','-fillpage')
%% SI calc across alternative n cycles
respLim = [0 0.035];
figure
suptitle('among modulated responsive cells...')
subplot 231
ind = respCells & lateModAV & lateSI > 0;
y = mean(respV(ind,:),1);
yerr = ste(respV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateModAV & lateSI > 0;
y = mean(respA(ind,:),1);
yerr = ste(respA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells (cyc 5-8), V > A (%s)', num2str(sum(ind))))
legend({'visual','auditory'},'location','southwest')

cyc2toEndModAV = cell2mat(ttestCyc2toEndAVEaExpt');
cyc2toEndSI = cell2mat(cyc2toEndSIEaExpt');

subplot 232
ind = BLRespCells & cyc2toEndModAV & cyc2toEndSI > 0;
y = mean(respV(ind,:),1);
yerr = ste(respV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = BLRespCells & cyc2toEndModAV & cyc2toEndSI > 0;
y = mean(respA(ind,:),1);
yerr = ste(respA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells (cyc 2-8), V > A (%s)', num2str(sum(ind))))

allCycModAV = cell2mat(ttestAllCycAVEaExpt');
allCycSI = cell2mat(allCycSIEaExpt');

subplot 233
ind = BLRespCells & allCycModAV & allCycSI > 0;
y = mean(respV(ind,:),1);
yerr = ste(respV(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = BLRespCells & allCycModAV & allCycSI > 0;
y = mean(respA(ind,:),1);
yerr = ste(respA(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells (all cyc), V > A (%s)', num2str(sum(ind))))


subplot 234
histogram(cell2mat(lateTtestP(:,1)'),20)
figXAxis([],'p-value',[-0.1 1.1])
figYAxis([],'n cells',[0 110])
figAxForm([])
title('t-test stimuli 5-8')

subplot 235
histogram(cell2mat(lateTtestP(:,2)'),20)
figXAxis([],'p-value',[-0.1 1.1])
figYAxis([],'n cells',[0 110])
figAxForm([])
title('t-test stimuli 2-8')

subplot 236
histogram(cell2mat(lateTtestP(:,3)'),20)
figXAxis([],'p-value',[-0.1 1.1])
figYAxis([],'n cells',[0 110])
figAxForm([])
title('t-test all stimuli')

print([fnout 'cycRespModCells_multiModWin'],'-dpdf','-fillpage')

figure
subplot 121
ind = BLRespCells & cyc2toEndModAV & cyc2toEndSI > 0;
h = bar(1,mean(earlyLateVNorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{1};
hold on
errorbar(1,mean(earlyLateVNorm(ind,2)),ste(earlyLateVNorm(ind,2),1),'k');
h = bar(2,mean(earlyLateANorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{2};
hold on
errorbar(2,mean(earlyLateANorm(ind,2)),ste(earlyLateANorm(ind,2),1),'c');
figXAxis([],'trial type',[0 3],1:2,{'visual';'auditory'})
figYAxis([],'% adapt',respLim)
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title(sprintf('V > A (%s)',num2str(sum(ind))))
print([fnout 'adaptRespModCells_multiModWin'],'-dpdf','-fillpage')
%% selectivity of suppressed cells (selectivity calc by late resp window for cycles 4-7)
lateRespSI = cell2mat(lateWinSIEaExpt');
lateRespModAV = cell2mat(ttestLateWinAVEaExpt');
figure;
h = histogram(lateRespSI(lateBLSuppCells),binEdgesSI);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [1 1 1];
h = vline(0,'k-');
h.LineWidth = 2;
h = vline(mean(lateRespSI(lateBLSuppCells)),'r-');
h.LineWidth = 2;
hold on
h = histogram(lateRespSI(lateRespModAV & lateBLSuppCells),binEdgesSI);
h.FaceColor = 'k';
h.EdgeColor = [1 1 1];
figXAxis([],'late resp selectivity index (using cycles 5-8)',siLim)
figYAxis([],'n ncells',[])
figAxForm([])
title({sprintf('significantly modulated late (%s/%s)'...
    ,num2str(sum(lateRespModAV & lateBLSuppCells)),num2str(sum(lateBLSuppCells))),...
    'red: avg SI all cels, pink: avg SI signif cells'})
h = vline(mean(lateRespSI(lateBLSuppCells & lateRespModAV)),'m-');
h.LineWidth = 2;
print([fnout 'SISuppModCellsLateWin'],'-dpdf','-fillpage')
%% selectivity across late window for responsive and suppressed cells
lateWinSI = cell2mat(lateWinSIEaExpt');
lateWinModAV = cell2mat(ttestLateWinAVEaExpt');

setFigParams4Print('portrait')
binEdgesSI = -10:0.5:10;
siLim = [-10 10];
figure
subplot 211
ind = BLRespCells | lateBLRespCells;
h=histogram(lateWinSI(ind),binEdgesSI);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [1 1 1];
h = vline(0,'k-');
h.LineWidth = 2;
h = vline(mean(lateWinSI(ind)),'r-');
h.LineWidth = 2;
hold on
h = histogram(lateWinSI(lateWinModAV & ind),binEdgesSI);
h.FaceColor = 'k';
h.EdgeColor = [1 1 1];
h = vline(mean(lateWinSI(ind & lateWinModAV)),'m-');
h.LineWidth = 2;
figXAxis([],'late trial selectivity index (using cycles 5-8)',siLim)
figYAxis([],'n ncells',[])
figAxForm([])
title({sprintf('responsive cells (%s/%s)'...
    ,num2str(sum(lateWinModAV & ind)),num2str(sum(ind)));...
    'red: avg SI all cels, pink: avg SI signif cells'})

subplot 212
ind = lateBLSuppCells;
h=histogram(lateWinSI(ind),binEdgesSI);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [1 1 1];
h = vline(0,'k-');
h.LineWidth = 2;
h = vline(mean(lateWinSI(ind)),'r-');
h.LineWidth = 2;
hold on
h = histogram(lateWinSI(lateWinModAV & ind),binEdgesSI);
h.FaceColor = 'k';
h.EdgeColor = [1 1 1];
h = vline(mean(lateWinSI(ind & lateWinModAV)),'m-');
h.LineWidth = 2;
figXAxis([],'late trial selectivity index (using cycles 5-8)',siLim)
figYAxis([],'n ncells',[])
figAxForm([])
title(sprintf('suppressed cells (%s/%s)'...
    ,num2str(sum(lateWinModAV & ind)),num2str(sum(ind))))
print([fnout '_SICellGroupsLateWin'],'-dpdf','-fillpage')
%% tc subtractions and shuffled subtractions
subLim = [-0.011 0.011];
tcAllCells = cell(1,nav);
tcShuffleAllCells = cell(1,nav);
for iav = 1:nav
    tcAllCells{iav} = cell2mat(tc(:,iav)');
    tcShuffleAllCells{iav} = cell2mat(trialTCShuffleEaExpt(:,iav)');
end
tcSub = tcAllCells{visualTrials} - tcAllCells{auditoryTrials};
tcShuffleSub = tcShuffleAllCells{visualTrials} - tcShuffleAllCells{auditoryTrials};

setFigParams4Print('portrait')
figure;
suptitle('subtraction (black) and shufffled subtraction (blue)')
subplot 311
y = mean(tcSub(:,BLRespCells),2);
yerr = ste(tcSub(:,BLRespCells),2);    
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcShuffleSub(:,BLRespCells),2);
yerr = ste(tcShuffleSub(:,BLRespCells),2);    
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'b');
hold on
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',subLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
hold on
hline(0,'k--');
title('respond to 1st')

subplot 312
y = mean(tcSub(:,lateBLRespCells),2);
yerr = ste(tcSub(:,lateBLRespCells),2);    
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcShuffleSub(:,lateBLRespCells),2);
yerr = ste(tcShuffleSub(:,lateBLRespCells),2);    
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'b');
hold on
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',subLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
hold on
hline(0,'k--');
title('respond late')

subplot 313
y = mean(tcSub(:,lateBLSuppCells),2);
yerr = ste(tcSub(:,lateBLSuppCells),2);    
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcShuffleSub(:,lateBLSuppCells),2);
yerr = ste(tcShuffleSub(:,lateBLSuppCells),2);    
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'b');
hold on
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',subLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
hold on
hline(0,'k--');
title('suppressed late')

print([fnout '_tcSubtractionCellGroups'],'-dpdf','-fillpage')

%% 
responseLim = [-0.035 0.05];
figure
suptitle('visual (black), auditory (cyan)')

subplot 221
ind = BLRespCells | lateBLRespCells;
y = mean(tcAllCells{visualTrials}(:,ind),2);
yerr = ste(tcAllCells{visualTrials}(:,ind),2);
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcAllCells{auditoryTrials}(:,ind),2);
yerr = ste(tcAllCells{auditoryTrials}(:,ind),2);
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'c');
hold on
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',responseLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
hold on
hline(0,'k--');
title('all responsive cells')

subplot 222
y = mean(tcAllCells{visualTrials}(:,BLRespCells),2);
yerr = ste(tcAllCells{visualTrials}(:,BLRespCells),2);
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcAllCells{auditoryTrials}(:,BLRespCells),2);
yerr = ste(tcAllCells{auditoryTrials}(:,BLRespCells),2);
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'c');
hold on
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',responseLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
hold on
hline(0,'k--');
title('respond to 1st')

subplot 223
y = mean(tcAllCells{visualTrials}(:,lateBLRespCells),2);
yerr = ste(tcAllCells{visualTrials}(:,lateBLRespCells),2);
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcAllCells{auditoryTrials}(:,lateBLRespCells),2);
yerr = ste(tcAllCells{auditoryTrials}(:,lateBLRespCells),2);
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'c');
hold on
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',responseLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
hold on
hline(0,'k--');
title('respond late')


subplot 224
y = mean(tcAllCells{visualTrials}(:,lateBLSuppCells),2);
yerr = ste(tcAllCells{visualTrials}(:,lateBLSuppCells),2);
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcAllCells{auditoryTrials}(:,lateBLSuppCells),2);
yerr = ste(tcAllCells{auditoryTrials}(:,lateBLSuppCells),2);
tt = ((26:length(y))-33).*(1000/frRateHz);
shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'c');
hold on
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',responseLim)  
f = gca;
f.Box = 'off';
f.TickDir = 'out';
hold on
hline(0,'k--');
title('suppressed late')

print([fnout '_tcAVCellGroups'],'-dpdf','-fillpage')
%% tuning of mod cells
tuningPeakAllCells = cell2mat(tuningPeak);
isTunedAllCells = cell2mat(isTuned');

binEdgesOri = 0:22.5:180;
nBinnedRespCells = histcounts(tuningPeakAllCells(respCells & isTunedAllCells)-1,binEdgesOri);
[b, k, R, u, ~, R_squareRespCells] = miaovonmisesfit_ori(deg2rad(binEdgesOri(1:end-1)),nBinnedRespCells);
popFitRespCells = b+R.*exp(k.*(cos(2.*(deg2rad(binEdgesOri(2:end))-u))-1));

nBinnedModCells = histcounts(...
    tuningPeakAllCells(respCells & lateModAV & lateSI > 0 & isTunedAllCells)-1,binEdgesOri);
[b, k, R, u, ~, R_squareModCells] = miaovonmisesfit_ori(deg2rad(binEdgesOri(1:end-1)),nBinnedModCells);
popFitModCells = b+R.*exp(k.*(cos(2.*(deg2rad(binEdgesOri(2:end))-u))-1));

figure;
h = bar(binEdgesOri(1:end-1),nBinnedRespCells);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [1 1 1];
h.BarWidth = 1;
hold on
h = plot(binEdgesOri(1:end-1),popFitRespCells,'r-');
h.LineWidth = 2;
hold on
h = bar(binEdgesOri(1:end-1),nBinnedModCells);
h.FaceColor = [0 0 0];
h.EdgeColor = [1 1 1];
h.BarWidth = 1;
hold on
h = plot(binEdgesOri(1:end-1),popFitModCells,'m-');
h.LineWidth = 2;
title({'tuning of all anticipation responsive/suppressed cells (gray) and those that are modulated (black)';...
    sprintf('von mises fit to population, red: responsive cells (R2 %s),pink: mod cells (%s)',...
    num2str(R_squareRespCells),num2str(R_squareModCells))})
print([fnout 'tuningRespAndModCells'],'-dpdf','-fillpage')

%%
tuningMod90AllCells = tuningPeakAllCells;
tuningMod90AllCells(tuningMod90AllCells > 90) = 180 - tuningMod90AllCells(tuningMod90AllCells > 90);

n = cell(1,4);
figure
subplot 121
cdfplot(tuningPeakAllCells(respCells & isTunedAllCells))
n{1} = num2str(sum(respCells & isTunedAllCells));
hold on
cdfplot(tuningPeakAllCells(respCells & isTunedAllCells & lateSI > 0 & lateModAV))
n{2} = num2str(sum(respCells & isTunedAllCells & lateSI > 0 & lateModAV));
hold on
cdfplot(tuningPeakAllCells(respCells & isTunedAllCells & lateSI < 0 & lateModAV))
n{3} = num2str(sum(respCells & isTunedAllCells & lateSI < 0 & lateModAV));
hold on
cdfplot(tuningPeakAllCells(respCells & isTunedAllCells & ~lateModAV))
n{4} = num2str(sum(respCells & isTunedAllCells & ~lateModAV));
[ks180Result,ks180pval] = kstest2(...
    tuningPeakAllCells(respCells & isTunedAllCells & lateSI > 0 & lateModAV),...
    tuningPeakAllCells(respCells & isTunedAllCells & ~lateModAV));
title({'distribution of orientation tuning';['compare vis selective to NS,ks test p = ' num2str(ks180pval)]})
legend(cellfun(@(x,y) sprintf(x,y),...
    {'tuned responsive cells (%s)','visual selective (%s)','auditory selective (%s)','non-selective (%s)'}...
    ,n,'unif',0),'location','northwest')
figXAxis([],'pref orientation',[0 180],[0:45:180],[0:45:180])
figYAxis([],'fraction of cells',[])
figAxForm([])

subplot 122
cdfplot(tuningMod90AllCells(respCells & isTunedAllCells))
hold on
cdfplot(tuningMod90AllCells(respCells & isTunedAllCells & lateSI > 0 & lateModAV))
hold on
cdfplot(tuningMod90AllCells(respCells & isTunedAllCells & lateSI < 0 & lateModAV))
hold on
cdfplot(tuningMod90AllCells(respCells & isTunedAllCells & ~lateModAV))
[ks90Result,ks90pval] = kstest2(...
    tuningMod90AllCells(respCells & isTunedAllCells & lateSI > 0 & lateModAV),...
    tuningMod90AllCells(respCells & isTunedAllCells & ~lateModAV));
title({'folded tuning curve';['compare vis selective to NS,ks test p = ' num2str(ks90pval)]})
figXAxis([],'orientation',[0 90],[0:45:90],[0:45:90])
figYAxis([],'fraction of cells',[])
figAxForm([])
print([fnout 'tuningDistribution'],'-dpdf','-fillpage')


%% tuning of responsive cells
[x,y,oriBinIndex] = histcounts(tuningPeakAllCells,binEdgesOri);
nBins = length(binEdgesOri)-1;
nOri = nBins/2;
oriInd = cell(1,nOri);
for iori = 1:nOri
    if iori == 1
        oriInd{iori} = oriBinIndex == 1 | oriBinIndex == nBins;
    else
        oriInd{iori} = oriBinIndex == iori+(iori-2) | oriBinIndex == iori+(iori-1);
    end
end

orientations = [0 45 90 135];

respLim = [-0.005 0.035];
figure
suptitle({'responsive & tuned cells'; 'tuning bins selected by grouping 22.5 deg around tested orientations'})
[nrows,ncols] = optimizeSubplotDim(nOri);
for iori = 1:nOri
    subplot(nrows,ncols,iori)
    ind = respCells & oriInd{iori} & isTunedAllCells;
    y = mean(respV(ind,:),1);
    yerr = ste(respV(ind,:),1);
    h = errorbar(1:nCycles,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(respA(ind,:),1);
    yerr = ste(respA(ind,:),1);
    h = errorbar(1:nCycles,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];   
    figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',respLim)
    figAxForm([])
    title(sprintf('%s deg tuned(%s)',num2str(orientations(iori)), num2str(sum(ind)))) 
end
print([fnout 'cycRespTunedCells'],'-dpdf','-fillpage')


% respLim = [-0.005 0.025];
figure
suptitle('visual selective cells')
[nrows,ncols] = optimizeSubplotDim(nOri);
for iori = 1:nOri
    subplot(nrows,ncols,iori)
    ind = respCells & oriInd{iori} & lateSI > 0 & lateModAV & isTunedAllCells;
    y = mean(respV(ind,:),1);
    yerr = ste(respV(ind,:),1);
    h = errorbar(1:nCycles,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(respA(ind,:),1);
    yerr = ste(respA(ind,:),1);
    h = errorbar(1:nCycles,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];   
    figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',respLim)
    figAxForm([])
    title(sprintf('%s deg tuned(%s)',num2str(orientations(iori)), num2str(sum(ind)))) 
end
print([fnout 'cycRespModAndTunedCells'],'-dpdf','-fillpage')

%% normalized responses of mod cells
respVNorm = respV./respV(:,1);
respANorm = respA./respA(:,1);

respLim = [-0.1 1.5];
figure
suptitle('visual selective cells, normalized to 1st stim')
[nrows,ncols] = optimizeSubplotDim(nOri);
for iori = 1:nOri
    subplot(nrows,ncols,iori)
    ind = BLRespCells & oriInd{iori} & lateSI > 0 & lateModAV & isTunedAllCells;
    y = mean(mean(respVNorm(ind,1),2),1);
    yerr = ste(mean(respVNorm(ind,1:2),2),1);
    h = errorbar(1,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(mean(respANorm(ind,1:2),2),1);
    yerr = ste(mean(respANorm(ind,1:2),2),1);
    h = errorbar(1,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(mean(respVNorm(ind,3:8),2),1);
    yerr = ste(mean(respVNorm(ind,3:8),2),1);
    h = errorbar(2,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(mean(respANorm(ind,3:8),2),1);
    yerr = ste(mean(respANorm(ind,3:6),2),1);
    h = errorbar(2,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];   
    figXAxis([],'stimulus number',[0 3],[1 2],{'1-2','3-8'})
    figYAxis([],'normalized dF/F',respLim)
    figAxForm([])
    title(sprintf('%s deg tuned(%s)',num2str(orientations(iori)), num2str(sum(ind)))) 
end
print([fnout 'earlyLateRespModAndTunedCellsNormalized'],'-dpdf','-fillpage')
%%
figure
suptitle('responsive cells, normalized to 1st stim')
[nrows,ncols] = optimizeSubplotDim(nOri);
for iori = 1:nOri
    subplot(nrows,ncols,iori)
    ind = BLRespCells & oriInd{iori} & isTunedAllCells;
    y = mean(mean(respVNorm(ind,1:2),2),1);
    yerr = ste(mean(respVNorm(ind,1:2),2),1);
    h = errorbar(1,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(mean(respANorm(ind,1:2),2),1);
    yerr = ste(mean(respANorm(ind,1:2),2),1);
    h = errorbar(1,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(mean(respVNorm(ind,3:8),2),1);
    yerr = ste(mean(respVNorm(ind,3:8),2),1);
    h = errorbar(2,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(mean(respANorm(ind,3:8),2),1);
    yerr = ste(mean(respANorm(ind,3:8),2),1);
    h = errorbar(2,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];   
    figXAxis([],'stimulus number',[0 3],[1 2],{'1','3-8'})
    figYAxis([],'normalized dF/F',respLim)
    figAxForm([])
    title(sprintf('%s deg tuned(%s)',num2str(orientations(iori)), num2str(sum(ind)))) 
end
print([fnout 'earlyLateRespCellsNormalized'],'-dpdf','-fillpage')

respLim = [-0.5 2];
figure
suptitle('visual selective cells, normalized to 1st stim')
[nrows,ncols] = optimizeSubplotDim(nOri);
for iori = 1:nOri
    subplot(nrows,ncols,iori)
    ind = BLRespCells & oriInd{iori} & lateSI > 0 & lateModAV & isTunedAllCells;
    y = mean(respVNorm(ind,:),1);
    yerr = ste(respVNorm(ind,:),1);
    h = errorbar(1:nCycles,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(respANorm(ind,:),1);
    yerr = ste(respANorm(ind,:),1);
    h = errorbar(1:nCycles,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];   
    figXAxis([],'stimulus number',[0 nCycles+1])
    figYAxis([],'normalized dF/F',respLim)
    figAxForm([])
    title(sprintf('%s deg tuned(%s)',num2str(orientations(iori)), num2str(sum(ind)))) 
end
print([fnout 'cycRespModAndTunedCellsNormalized'],'-dpdf','-fillpage')

figure
suptitle({'responsive & tuned cells'; 'tuning bins selected by grouping 22.5 deg around tested orientations'})
[nrows,ncols] = optimizeSubplotDim(nOri);
for iori = 1:nOri
    subplot(nrows,ncols,iori)
    ind = BLRespCells & oriInd{iori} & isTunedAllCells;
    y = mean(respVNorm(ind,:),1);
    yerr = ste(respVNorm(ind,:),1);
    h = errorbar(1:nCycles,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(respANorm(ind,:),1);
    yerr = ste(respANorm(ind,:),1);
    h = errorbar(1:nCycles,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];   
    figXAxis([],'stimulus number',[0 nCycles+1])
    figYAxis([],'normalized dF/F',respLim)
    figAxForm([])
    title(sprintf('%s deg tuned(%s)',num2str(orientations(iori)), num2str(sum(ind)))) 
end
print([fnout 'cycRespTunedCellsNormalized'],'-dpdf','-fillpage')

figure
suptitle('among anticipation responsive or suppressed cells (normalized)...')
subplot 121
ind = respCells & lateModAV & lateSI > 0;
y = mean(respVNorm(ind,:),1);
yerr = ste(respVNorm(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateModAV & lateSI > 0;
y = mean(respANorm(ind,:),1);
yerr = ste(respANorm(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'stimulus number',[0 nCycles+1])
figYAxis([],'normalized dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells, V > A (%s)', num2str(sum(ind))))
legend({'visual','auditory'},'location','southwest')

subplot 122
ind = respCells & lateModAV & lateSI < 0;
y = mean(respVNorm(ind,:),1);
yerr = ste(respVNorm(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'ko');
h.MarkerFaceColor = [1 1 1];
hold on
ind = respCells & lateModAV & lateSI < 0;
y = mean(respANorm(ind,:),1);
yerr = ste(respANorm(ind,:),1);
h = errorbar(1:nCycles,y,yerr,'co');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'stimulus number',[0 nCycles+1])
figYAxis([],'normalized dF/F',respLim)
figAxForm([])
title(sprintf('late modulated cells, V < A (%s)', num2str(sum(ind))))
print([fnout 'cycRespModCellsNormalized'],'-dpdf','-fillpage')

%% compare early and late
earlyLateV = cell2mat(respEarlyLateV)';
earlyLateA = cell2mat(respEarlyLateA)';

earlyLateVNorm = earlyLateV./earlyLateV(:,1);
earlyLateANorm = earlyLateA./earlyLateA(:,1);

respLim = [0 1.1];
figure;
subplot 121
ind = BLRespCells;
h = errorbar(mean(earlyLateVNorm(ind,:)),ste(earlyLateVNorm(ind,:),1),'ko-');
h.MarkerFaceColor = [1 1 1];
hold on
h = errorbar(mean(earlyLateANorm(ind,:)),ste(earlyLateANorm(ind,:),1),'co-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'stimulus number',[0 3],1:2,{'1','3-8'})
figYAxis([],'normalized dF/F',respLim)
figAxForm([])
title(['responsive cells, n = ' num2str(sum(ind))])
subplot 122
ind = BLRespCells & lateSI > 0 & lateModAV;
h = errorbar(mean(earlyLateVNorm(ind,:)),ste(earlyLateVNorm(ind,:),1),'ko-');
h.MarkerFaceColor = [1 1 1];
hold on
h = errorbar(mean(earlyLateANorm(ind,:)),ste(earlyLateANorm(ind,:),1),'co-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'stimulus number',[0 3],1:2,{'1','3-8'})
figYAxis([],'normalized dF/F',respLim)
figAxForm([])
title(['late modulated cells, n = ' num2str(sum(ind))])
print([fnout 'earlyLateRespNormalized'],'-dpdf','-fillpage')

figure;
subplot 121
ind = BLRespCells;
h = bar(1,mean(earlyLateVNorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{1};
hold on
errorbar(1,mean(earlyLateVNorm(ind,2)),ste(earlyLateVNorm(ind,2),1),'k');
h = bar(2,mean(earlyLateANorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{2};
hold on
errorbar(2,mean(earlyLateANorm(ind,2)),ste(earlyLateANorm(ind,2),1),'c');
figXAxis([],'trial type',[0 3],1:2,{'visual';'auditory'})
figYAxis([],'% adapt',respLim)
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title(sprintf('all resp cells (%s)',num2str(sum(ind))))
subplot 122
ind = BLRespCells & lateSI > 0 & lateModAV;
h = bar(1,mean(earlyLateVNorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{1};
hold on
errorbar(1,mean(earlyLateVNorm(ind,2)),ste(earlyLateVNorm(ind,2),1),'k');
h = bar(2,mean(earlyLateANorm(ind,2)));
h.BarWidth = 1;
h.EdgeColor = [1 1 1];
h.FaceColor = avcolmat{2};
hold on
errorbar(2,mean(earlyLateANorm(ind,2)),ste(earlyLateANorm(ind,2),1),'c');
figXAxis([],'trial type',[0 3],1:2,{'visual';'auditory'})
figYAxis([],'% adapt',respLim)
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title(sprintf('modulated cells (%s)',num2str(sum(ind))))
print([fnout 'earlyLateRespNormalizedBar'],'-dpdf','-fillpage')



modRespCells = BLRespCells & lateSI > 0 & lateModAV & isTunedAllCells;

pctTunedCellsMod = nan(1,nOri+2);
pctTunedCellsMod(1) = sum(modRespCells) / sum(BLRespCells & isTunedAllCells);
pctTunedCellsMod(end) = sum(BLRespCells & lateSI > 0 & lateModAV) /sum(BLRespCells)

meanPctOfFirstStim = nan(nav,nOri+1);
semPctOfFirstStim = nan(nav,nOri+1);

meanPctOfFirstStim(visualTrials,1) = mean(earlyLateVNorm(BLRespCells & isTunedAllCells,2),1);
meanPctOfFirstStim(auditoryTrials,1) = mean(earlyLateANorm(BLRespCells & isTunedAllCells,2),1);
semPctOfFirstStim(visualTrials,1) = ste(earlyLateVNorm(BLRespCells & isTunedAllCells,2),1);
semPctOfFirstStim(auditoryTrials,1) = ste(earlyLateANorm(BLRespCells & isTunedAllCells,2),1);

meanModPctOfFirstStim = nan(nav,nOri+1);
semModPctOfFirstStim = nan(nav,nOri+1);

meanModPctOfFirstStim(visualTrials,1) = mean(earlyLateVNorm(modRespCells,2),1);
meanModPctOfFirstStim(auditoryTrials,1) = mean(earlyLateANorm(modRespCells,2),1);
semModPctOfFirstStim(visualTrials,1) = ste(earlyLateVNorm(modRespCells,2),1);
semModPctOfFirstStim(auditoryTrials,1) = ste(earlyLateANorm(modRespCells,2),1);

cellGroupsName = {'all resp tuned';'0';'45';'90';'135'};

SIcdfFig = figure;
hold on
cdfplot(lateSI(BLRespCells))

tuningCDFRespCells = cell(nOri,1);
tuningCDFSuppCells = cell(nOri,1);

normRespOriFig = figure;
suptitle('responsive cells, normalized to 1st stim')
[nrows,ncols] = optimizeSubplotDim(nOri);
for iori = 1:nOri
    
    
    figure(normRespOriFig)
    subplot(nrows,ncols,iori)
    ind = BLRespCells & oriInd{iori} & isTunedAllCells;
    
    meanPctOfFirstStim(visualTrials,iori+1) = mean(earlyLateVNorm(ind,2),1);
    meanPctOfFirstStim(auditoryTrials,iori+1) = mean(earlyLateANorm(ind,2),1);
    semPctOfFirstStim(visualTrials,iori+1) = ste(earlyLateVNorm(ind,2),1);
    semPctOfFirstStim(auditoryTrials,iori+1) = ste(earlyLateANorm(ind,2),1);
    
    ind2 = modRespCells & oriInd{iori};
    pctTunedCellsMod(iori+1) = sum(ind2) / sum(BLRespCells & isTunedAllCells & oriInd{iori});
    
    meanModPctOfFirstStim(visualTrials,iori+1) = mean(earlyLateVNorm(ind2,2),1);
    meanModPctOfFirstStim(auditoryTrials,iori+1) = mean(earlyLateANorm(ind2,2),1);
    semModPctOfFirstStim(visualTrials,iori+1) = ste(earlyLateVNorm(ind2,2),1);
    semModPctOfFirstStim(auditoryTrials,iori+1) = ste(earlyLateANorm(ind2,2),1);
    
    y = mean(earlyLateVNorm(ind,1),1);
    yerr = ste(earlyLateVNorm(ind,1),1);
    h = errorbar(1,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(earlyLateANorm(ind,1),1);
    yerr = ste(earlyLateANorm(ind,1),1);
    h = errorbar(1,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(earlyLateVNorm(ind,2),1);
    yerr = ste(earlyLateVNorm(ind,2),1);
    h = errorbar(2,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    y = mean(earlyLateANorm(ind,2),1);
    yerr = ste(earlyLateANorm(ind,2),1);
    h = errorbar(2,y,yerr,'co');
    h.MarkerFaceColor = [1 1 1];
    hold on   
    figXAxis([],'stimulus number',[0 3],[1 2],{'1','3-8'})
    figYAxis([],'normalized dF/F',respLim)
    figAxForm([])
    title(sprintf('%s deg tuned(%s)',num2str(orientations(iori)), num2str(sum(ind)))) 
    
    figure(SIcdfFig)
    hold on
    cdfplot(lateWinSI(ind))
    
    tuningCDFRespCells{iori} = lateWinSI(ind);
    tuningCDFSuppCells{iori} = lateWinSI(lateBLSuppCells & oriInd{iori} & isTunedAllCells);
end
figure(normRespOriFig)
print([fnout 'earlyLateRespCellsNormalized'],'-dpdf','-fillpage')

siLim = [-4 5];
figure;
subplot 221
for iori = 1:nOri
    cdfplot(tuningCDFRespCells{iori})
    hold on
end
legend({'0';'45';'90';'135'})
title('resp')
figXAxis([],'selectivity',siLim)
figYAxis([],'fraction of cells',[0 1])
figAxForm([])
subplot 222
cdfplot(tuningCDFRespCells{2})
hold on
cdfplot(cell2mat(tuningCDFRespCells([1,3,4])'))
[~,p] = kstest2(tuningCDFRespCells{2},cell2mat(tuningCDFRespCells([1,3,4])'));
legend({'45';'other'})
title(['resp; p = ' num2str(p)])
figXAxis([],'selectivity',siLim)
figYAxis([],'fraction of cells',[0 1])
figAxForm([])

subplot 223
for iori = 1:nOri
    cdfplot(tuningCDFSuppCells{iori})
    hold on
end
legend({'0';'45';'90';'135'})
title('supp')
figXAxis([],'selectivity',siLim)
figYAxis([],'fraction of cells',[0 1])
figAxForm([])
subplot 224
cdfplot(tuningCDFSuppCells{2})
hold on
cdfplot(cell2mat(tuningCDFSuppCells([1,3,4])'))
[~,p] = kstest2(tuningCDFSuppCells{2},cell2mat(tuningCDFSuppCells([1,3,4])'));
legend({'45';'other'})
title(['supp; p = ' num2str(p)])
figXAxis([],'selectivity',siLim)
figYAxis([],'fraction of cells',[0 1])
figAxForm([])
print([fnout 'cdfSIOriRespSupp'],'-dpdf','-fillpage')


figure(SIcdfFig)
figXAxis([],'selectivity index',[])
figYAxis([],'fraction of cells in group',[0 1])
hold on
vline(0,'k--')
title('distribution of SI for tuning groups')
legend(cellGroupsName)
print([fnout 'cdfSIOriCellGroups'],'-dpdf','-fillpage')

figure
avcolmat = {'k','c'};
h = bar(1:nOri+1,meanPctOfFirstStim','grouped');
xoffset = [-0.1 0.1];
for iav = 1:nav
    h(iav).BarWidth = 1;
    h(iav).EdgeColor = [1 1 1];
    h(iav).FaceColor = avcolmat{iav};
    hold on
    
    errpoint = errorbar((1:nOri+1)+xoffset(iav),meanPctOfFirstStim(iav,:),semPctOfFirstStim(iav,:),avcolmat{iav});
    errpoint.LineStyle = 'none';
end
figXAxis([],'cell group',[0 nOri+2],[1:nOri+1],cellGroupsName)
figYAxis([],'% adapted',[-0.5 1]);
figAxForm([])
title('percent adaptation from first stim to late stim within cell groups')
legend({'visual';'auditory'},'location','northeastoutside')
print([fnout 'pctAdaptationCellGroupsBar'],'-dpdf','-fillpage')

figure
avcolmat = {'k','c'};
h = bar(1:nOri+1,meanModPctOfFirstStim','grouped');
xoffset = [-0.1 0.1];
for iav = 1:nav
    h(iav).BarWidth = 1;
    h(iav).EdgeColor = [1 1 1];
    h(iav).FaceColor = avcolmat{iav};
    hold on
    
    errpoint = errorbar((1:nOri+1)+xoffset(iav),meanModPctOfFirstStim(iav,:),semModPctOfFirstStim(iav,:),avcolmat{iav});
    errpoint.LineStyle = 'none';
end
figXAxis([],'cell group',[0 nOri+2],[1:nOri+1],cellGroupsName)
figYAxis([],'% adapted',[0 1]);
figAxForm([])
title('percent adaptation from first stim to late stim within cell groups')
legend({'visual';'auditory'},'location','northeastoutside')
print([fnout 'pctModAdaptationCellGroupsBar'],'-dpdf','-fillpage')


figure
for iav = 1:nav    
    hold on
    errpoint = errorbar((1:nOri+1),meanPctOfFirstStim(iav,:),semPctOfFirstStim(iav,:),[avcolmat{iav} 'o']);
    errpoint.LineStyle = 'none';
    errpoint.MarkerFaceColor = [1 1 1];
end
figXAxis([],'cell group',[0 nOri+2],[1:nOri+1],cellGroupsName)
figYAxis([],'% adapted',[0 1]);
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title('percent adaptation from first stim to late stim within cell groups')
legend({'visual';'auditory'},'location','northeastoutside')
print([fnout 'pctAdaptationCellGroups'],'-dpdf','-fillpage')

figure
for iav = 1:nav    
    hold on
    errpoint = errorbar((1:nOri+1),meanModPctOfFirstStim(iav,:),semModPctOfFirstStim(iav,:),[avcolmat{iav} 'o']);
    errpoint.LineStyle = 'none';
    errpoint.MarkerFaceColor = [1 1 1];
end
figXAxis([],'cell group',[0 nOri+2],[1:nOri+1],cellGroupsName)
figYAxis([],'% adapted',[0 1]);
f = gca;
f.Box = 'off';
f.TickDir = 'out';
title('percent adaptation from first stim to late stim within cell groups')
legend({'visual';'auditory'},'location','northeastoutside')
print([fnout 'pctAdaptationModCellGroups'],'-dpdf','-fillpage')

%pct of tuned modulated cells
figure;
subplot 121
h = barh(1:nOri+1,fliplr(pctTunedCellsMod(1:5)*100));
h.FaceColor = [1 1 1];
figXAxis([],'% modulated cells in group',[0 25]);
figYAxis([],'cell group',[],1:nOri+1,fliplr(cellGroupsName))
figAxForm([])
title(['% = ' num2str(round(fliplr(pctTunedCellsMod(1:5)*100),1))])
subplot 122
h = bar(1:nOri+1,pctTunedCellsMod(1:5)*100);
h.FaceColor = [1 1 1];
figYAxis([],'% modulated cells in group',[0 25]);
figXAxis([],'cell group',[],1:nOri+1,cellGroupsName)
figAxForm([])
print([fnout 'pctModulatedCellsInGroup'],'-dpdf','-fillpage')

% %% previous trial type
% AVmodCells = lateSI > 0 & lateModAV;
% 
% avgRespEaCycPrevVisAV = cellfun(@(x) mean(x(BLRespCells)),respPrevVisEaCycAV);
% avgRespEaCycPrevAudAV = cellfun(@(x) mean(x(BLRespCells)),respPrevAudEaCycAV);
% semRespEaCycPrevVisAV = cellfun(@(x) ste(x(BLRespCells),2),respPrevVisEaCycAV);
% semRespEaCycPrevAudAV = cellfun(@(x) ste(x(BLRespCells),2),respPrevAudEaCycAV);
% 
% avTitleName = {'visual trials';'auditory trials'};
% legLabels = {'prev visual';'prev auditory'};
% prevAVColors = cat(3,[0 0 0; 0.5 0.5 0.5],[0 0.5 0.5; 0 1 1]);
% 
% respLim = [0 0.015];
% figure
% for iav = 1:nav
%     subplot(1,2,iav)
%     h = errorbar(1:nCycles,avgRespEaCycPrevVisAV(:,iav),semRespEaCycPrevVisAV(:,iav),'o');
%     h.Color = prevAVColors(visualTrials,:,iav);
%     h.MarkerFaceColor = [1 1 1];
%     hold on
%     h = errorbar(1:nCycles,avgRespEaCycPrevAudAV(:,iav),semRespEaCycPrevAudAV(:,iav),'o');
%     h.Color = prevAVColors(auditoryTrials,:,iav);
%     h.MarkerFaceColor = [1 1 1];
%     title(avTitleName(iav))
%     legend(legLabels)
%     figXAxis([],'stimulus number',[0 nCycles+1])
%     figYAxis([],'dF/F',respLim)
%     figAxForm([])
% end
% print([fnout 'prevTrialRespEaCyc'],'-dpdf')
% 
% avgRespPrevVisAV = cellfun(@(x) mean(x(BLRespCells)),respPrevVisAllCells);
% avgRespPrevAudAV = cellfun(@(x) mean(x(BLRespCells)),respPrevAudAllCells);
% semRespPrevVisAV = cellfun(@(x) ste(x(BLRespCells),2),respPrevVisAllCells);
% semRespPrevAudAV = cellfun(@(x) ste(x(BLRespCells),2),respPrevAudAllCells);
% 
% figure
% for iav = 1:nav
%     subplot(1,2,iav)
%     h = errorbar(1:2,avgRespPrevVisAV(:,iav),semRespPrevVisAV(:,iav),'o');
%     h.Color = prevAVColors(visualTrials,:,iav);
%     hold on
%     h = errorbar(1:2,avgRespPrevAudAV(:,iav),semRespPrevAudAV(:,iav),'o');
%     h.Color = prevAVColors(auditoryTrials,:,iav);
%     title(avTitleName(iav))
%     legend(legLabels,'location','southwest')
%     figXAxis([],'stimulus number',[0 3],1:2,{'1','3-8'})
%     figYAxis([],'dF/F',respLim)
%     figAxForm([])
% end
% 
% 
% avgRespEaCycPrevVisAV_AVmodcells = cellfun(@(x) mean(x(BLRespCells & AVmodCells)),respPrevVisEaCycAV);
% avgRespEaCycPrevAudAV_AVmodcells = cellfun(@(x) mean(x(BLRespCells & AVmodCells)),respPrevAudEaCycAV);
% semRespEaCycPrevVisAV_AVmodcells = cellfun(@(x) ste(x(BLRespCells & AVmodCells),2),respPrevVisEaCycAV);
% semRespEaCycPrevAudAV_AVmodcells = cellfun(@(x) ste(x(BLRespCells & AVmodCells),2),respPrevAudEaCycAV);
% 
% figure
% for iav = 1:nav
%     subplot(1,2,iav)
%     h = errorbar(1:nCycles,avgRespEaCycPrevVisAV_AVmodcells(:,iav),semRespEaCycPrevVisAV_AVmodcells(:,iav),'o');
%     h.Color = prevAVColors(visualTrials,:,iav);
%     h.MarkerFaceColor = [1 1 1];
%     hold on
%     h = errorbar(1:nCycles,avgRespEaCycPrevAudAV_AVmodcells(:,iav),semRespEaCycPrevAudAV_AVmodcells(:,iav),'o');
%     h.Color = prevAVColors(auditoryTrials,:,iav);
%     h.MarkerFaceColor = [1 1 1];
%     title(avTitleName(iav))
%     legend(legLabels)
%     figXAxis([],'stimulus number',[0 nCycles+1])
%     figYAxis([],'dF/F',respLim)
%     figAxForm([])
% end
% print([fnout 'prevTrialRespEaCycModCells'],'-dpdf')
% 
% avgRespPrevVisAV_AVmodcells = cellfun(@(x) mean(x(BLRespCells & AVmodCells)),respPrevVisAllCells);
% avgRespPrevAudAV_AVmodcells = cellfun(@(x) mean(x(BLRespCells & AVmodCells)),respPrevAudAllCells);
% semRespPrevVisAV_AVmodcells = cellfun(@(x) ste(x(BLRespCells & AVmodCells),2),respPrevVisAllCells);
% semRespPrevAudAV_AVmodcells = cellfun(@(x) ste(x(BLRespCells & AVmodCells),2),respPrevAudAllCells);
% 
% figure
% for iav = 1:nav
%     subplot(1,2,iav)
%     h = errorbar(1:2,avgRespPrevVisAV_AVmodcells(:,iav),semRespPrevVisAV_AVmodcells(:,iav),'o');
%     h.Color = prevAVColors(visualTrials,:,iav);
%     hold on
%     h = errorbar(1:2,avgRespPrevAudAV_AVmodcells(:,iav),semRespPrevAudAV_AVmodcells(:,iav),'o');
%     h.Color = prevAVColors(auditoryTrials,:,iav);
%     title(avTitleName(iav))
%     legend(legLabels,'location','southwest')
%     figXAxis([],'stimulus number',[0 3],1:2,{'1','3-8'})
%     figYAxis([],'dF/F',respLim)
%     figAxForm([])
% end
% 
% 
% 
% prevSILateCycVisTrials = cell2mat(siPrevAVLateCyc(:,visualTrials)');
% prevSILateCycAudTrials = cell2mat(siPrevAVLateCyc(:,auditoryTrials)');
% 
% prevModLateCycVisTrials = cell2mat(ttestPrevAVLateCyc(:,visualTrials))';
% prevModLateCycAudTrials = cell2mat(ttestPrevAVLateCyc(:,auditoryTrials))';
% 
% SILim = [-8 8];
% nCellsLim = [0 35];
% binEdgesSI = -8:0.5:8;
% 
% figure;
% suptitle({sprintf('previous trial significantly modulated cells (%s responsive cells)',...
%     num2str(sum(BLRespCells))),'cycles 3-8 used','red: avg SI all cels, pink: avg SI signif cells'})
% subplot 121
% h = histogram(prevSILateCycVisTrials(BLRespCells),binEdgesSI);
% h.FaceColor = [0.5 0.5 0.5];
% h.EdgeColor = [1 1 1];
% hold on
% h = histogram(prevSILateCycVisTrials(prevModLateCycVisTrials & BLRespCells),binEdgesSI);
% h.FaceColor = 'k';
% h.EdgeColor = [1 1 1];
% figXAxis([],{'prev trial selectivity index';'prev diff <-- --> prev same'},siLim)
% figYAxis([],'n ncells',nCellsLim)
% figAxForm([])
% h = vline(0,'k-');
% h.LineWidth = 2;
% h = vline(mean(prevSILateCycVisTrials(respCells)),'r-');
% h.LineWidth = 2;
% h = vline(mean(prevSILateCycVisTrials(prevModLateCycVisTrials & BLRespCells)),'m-');
% h.LineWidth = 2;
% title(sprintf('visual trials (%s modulated)',num2str(sum(...
%     prevModLateCycVisTrials & BLRespCells))))
% 
% subplot 122
% h = histogram(prevSILateCycAudTrials(BLRespCells),binEdgesSI);
% h.FaceColor = [0.5 0.5 0.5];
% h.EdgeColor = [1 1 1];
% hold on
% h = histogram(prevSILateCycAudTrials(prevModLateCycAudTrials & BLRespCells),binEdgesSI);
% h.FaceColor = 'k';
% h.EdgeColor = [1 1 1];
% figXAxis([],{'prev trial selectivity index';'prev same <-- --> prev diff'},siLim)
% figYAxis([],'n ncells',nCellsLim)
% figAxForm([])
% h = vline(0,'k-');
% h.LineWidth = 2;
% h = vline(mean(prevSILateCycAudTrials(respCells)),'r-');
% h.LineWidth = 2;
% h = vline(mean(prevSILateCycAudTrials(prevModLateCycAudTrials & BLRespCells)),'m-');
% h.LineWidth = 2;
% title(sprintf('auditory trials (%s modulated)',num2str(sum(...
%     prevModLateCycAudTrials & BLRespCells))))
% print([fnout 'prevTrialSIHist'],'-dpdf')
% 
% %% previous trial TC
% taskCells = BLRespCells | lateBLRespCells | lateBLSuppCells;
% trialTCVisPrevTrialEaExpt = squeeze(prevVisTrialTCEaCycAV(nCycles,:,:));
% trialTCAudPrevTrialEaExpt = squeeze(prevAudTrialTCEaCycAV(nCycles,:,:));
% 
% visPrevTrialTC = cell(1,2);
% visPrevTrialTC_sem = cell(1,2);
% audPrevTrialTC = cell(1,2);
% audPrevTrialTC_sem = cell(1,2);
% 
% tc = cell2mat(trialTCVisPrevTrialEaExpt(1,:));
% visPrevTrialTC{1} = mean(tc(:,taskCells),2);
% visPrevTrialTC_sem{1} = ste(tc(:,taskCells),2);
% 
% tc = cell2mat(trialTCAudPrevTrialEaExpt(1,:));
% visPrevTrialTC{2} = mean(tc(:,taskCells),2);
% visPrevTrialTC_sem{2} = ste(tc(:,taskCells),2);
% 
% tc = cell2mat(trialTCVisPrevTrialEaExpt(2,:));
% audPrevTrialTC{1} = mean(tc(:,taskCells),2);
% audPrevTrialTC_sem{1} = ste(tc(:,taskCells),2);
% 
% tc = cell2mat(trialTCAudPrevTrialEaExpt(2,:));
% audPrevTrialTC{2} = mean(tc(:,taskCells),2);
% audPrevTrialTC_sem{2} = ste(tc(:,taskCells),2);
% 
% 
% responseLim = [-0.004 0.015];
% prevtrcol = {[0 0 0],[.5 .5 .5];[0 1 1],[0 .5 .5]};
% prevtrpatch = {[0.85 0.85 0.85],[.95 .95 .95];[.85 1 1],[0.5 .95 .95]};
% setFigParams4Print('portrait')
% figure
% suptitle('all anticipation cells')
% for iav = 1:nav
%     subplot 211
%     y = visPrevTrialTC{iav};
%     yerr = visPrevTrialTC_sem{iav};
%     tt = ((26:length(y))-33).*(1000/frRateHz);
%     h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
%     h.mainLine.Color = prevtrcol{1,iav};
%     h.patch.FaceColor = prevtrpatch{1,iav};
%     title('visual')
%     figXAxis([],'time (ms)',[tt(1) tt(end)])
%     figYAxis([],'dF/F',responseLim) 
%     f = gca;
%     f.Box = 'off';
%     f.TickDir = 'out';
%     hold on
%     subplot 212
%     y = audPrevTrialTC{iav};
%     yerr = audPrevTrialTC_sem{iav};
%     tt = ((26:length(y))-33).*(1000/frRateHz);
%     h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'c');
%     title('auditory')
%     h.mainLine.Color = prevtrcol{2,iav};
%     h.patch.FaceColor = prevtrpatch{1,iav};
%     figXAxis([],'time (ms)',[tt(1) tt(end)])
%     f = gca;
%     f.Box = 'off';
%     f.TickDir = 'out';
%     figYAxis([],'dF/F',responseLim) 
%     hold on
% end
% 
% print([fnout '_tcPrevTrial_allAnticipationCells'],'-dpdf','-fillpage')

%% 1st stim responsive cells
taskCells = BLRespCells;
trialTCVisPrevTrialEaExpt = squeeze(prevVisTrialTCEaCycAV(nCycles,:,:));
trialTCAudPrevTrialEaExpt = squeeze(prevAudTrialTCEaCycAV(nCycles,:,:));

visPrevTrialTC = cell(1,2);
visPrevTrialTC_sem = cell(1,2);
audPrevTrialTC = cell(1,2);
audPrevTrialTC_sem = cell(1,2);

tc = cell2mat(trialTCVisPrevTrialEaExpt(1,:));
visPrevTrialTC{1} = mean(tc(:,taskCells),2);
visPrevTrialTC_sem{1} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCAudPrevTrialEaExpt(1,:));
visPrevTrialTC{2} = mean(tc(:,taskCells),2);
visPrevTrialTC_sem{2} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCVisPrevTrialEaExpt(2,:));
audPrevTrialTC{1} = mean(tc(:,taskCells),2);
audPrevTrialTC_sem{1} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCAudPrevTrialEaExpt(2,:));
audPrevTrialTC{2} = mean(tc(:,taskCells),2);
audPrevTrialTC_sem{2} = ste(tc(:,taskCells),2);


responseLim = [-0.004 0.03];
prevtrcol = {[0 0 0],[.5 .5 .5];[0 1 1],[0 .5 .5]};
prevtrpatch = {[0.85 0.85 0.85],[.95 .95 .95];[.85 1 1],[0.5 .95 .95]};
setFigParams4Print('portrait')
figure
suptitle('1st stim responsive cells')
for iav = 1:nav
    subplot 211
    y = visPrevTrialTC{iav};
    yerr = visPrevTrialTC_sem{iav};
    tt = ((26:length(y))-33).*(1000/frRateHz);
    h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
    h.mainLine.Color = prevtrcol{1,iav};
    h.patch.FaceColor = prevtrpatch{1,iav};
    title('visual')
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim) 
    f = gca;
    f.Box = 'off';
    f.TickDir = 'out';
    hold on
    subplot 212
    y = audPrevTrialTC{iav};
    yerr = audPrevTrialTC_sem{iav};
    tt = ((26:length(y))-33).*(1000/frRateHz);
    h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'c');
    title('auditory')
    h.mainLine.Color = prevtrcol{2,iav};
    h.patch.FaceColor = prevtrpatch{1,iav};
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    f = gca;
    f.Box = 'off';
    f.TickDir = 'out';
    figYAxis([],'dF/F',responseLim) 
    hold on
end

print([fnout '_tcPrevTrial_firstStimResp'],'-dpdf','-fillpage')

%% late resp
taskCells = lateBLRespCells;
trialTCVisPrevTrialEaExpt = squeeze(prevVisTrialTCEaCycAV(nCycles,:,:));
trialTCAudPrevTrialEaExpt = squeeze(prevAudTrialTCEaCycAV(nCycles,:,:));

visPrevTrialTC = cell(1,2);
visPrevTrialTC_sem = cell(1,2);
audPrevTrialTC = cell(1,2);
audPrevTrialTC_sem = cell(1,2);

tc = cell2mat(trialTCVisPrevTrialEaExpt(1,:));
visPrevTrialTC{1} = mean(tc(:,taskCells),2);
visPrevTrialTC_sem{1} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCAudPrevTrialEaExpt(1,:));
visPrevTrialTC{2} = mean(tc(:,taskCells),2);
visPrevTrialTC_sem{2} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCVisPrevTrialEaExpt(2,:));
audPrevTrialTC{1} = mean(tc(:,taskCells),2);
audPrevTrialTC_sem{1} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCAudPrevTrialEaExpt(2,:));
audPrevTrialTC{2} = mean(tc(:,taskCells),2);
audPrevTrialTC_sem{2} = ste(tc(:,taskCells),2);


responseLim = [-0.01 0.03];
prevtrcol = {[0 0 0],[.5 .5 .5];[0 1 1],[0 .5 .5]};
prevtrpatch = {[0.85 0.85 0.85],[.95 .95 .95];[.85 1 1],[0.5 .95 .95]};
setFigParams4Print('portrait')
figure
suptitle('late responsive cells')
for iav = 1:nav
    subplot 211
    y = visPrevTrialTC{iav};
    yerr = visPrevTrialTC_sem{iav};
    tt = ((26:length(y))-33).*(1000/frRateHz);
    h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
    h.mainLine.Color = prevtrcol{1,iav};
    h.patch.FaceColor = prevtrpatch{1,iav};
    title('visual')
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim) 
    f = gca;
    f.Box = 'off';
    f.TickDir = 'out';
    hold on
    subplot 212
    y = audPrevTrialTC{iav};
    yerr = audPrevTrialTC_sem{iav};
    tt = ((26:length(y))-33).*(1000/frRateHz);
    h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'c');
    title('auditory')
    h.mainLine.Color = prevtrcol{2,iav};
    h.patch.FaceColor = prevtrpatch{1,iav};
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    f = gca;
    f.Box = 'off';
    f.TickDir = 'out';
    figYAxis([],'dF/F',responseLim) 
    hold on
end

print([fnout '_tcPrevTrial_lateResp'],'-dpdf','-fillpage')
%% late supp
taskCells = lateBLSuppCells;
trialTCVisPrevTrialEaExpt = squeeze(prevVisTrialTCEaCycAV(nCycles,:,:));
trialTCAudPrevTrialEaExpt = squeeze(prevAudTrialTCEaCycAV(nCycles,:,:));

visPrevTrialTC = cell(1,2);
visPrevTrialTC_sem = cell(1,2);
audPrevTrialTC = cell(1,2);
audPrevTrialTC_sem = cell(1,2);

tc = cell2mat(trialTCVisPrevTrialEaExpt(1,:));
visPrevTrialTC{1} = mean(tc(:,taskCells),2);
visPrevTrialTC_sem{1} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCAudPrevTrialEaExpt(1,:));
visPrevTrialTC{2} = mean(tc(:,taskCells),2);
visPrevTrialTC_sem{2} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCVisPrevTrialEaExpt(2,:));
audPrevTrialTC{1} = mean(tc(:,taskCells),2);
audPrevTrialTC_sem{1} = ste(tc(:,taskCells),2);

tc = cell2mat(trialTCAudPrevTrialEaExpt(2,:));
audPrevTrialTC{2} = mean(tc(:,taskCells),2);
audPrevTrialTC_sem{2} = ste(tc(:,taskCells),2);


responseLim = [-0.04 0.01];
prevtrcol = {[0 0 0],[.5 .5 .5];[0 1 1],[0 .5 .5]};
prevtrpatch = {[0.85 0.85 0.85],[.95 .95 .95];[.85 1 1],[0.5 .95 .95]};
setFigParams4Print('portrait')
figure
suptitle('late supp cells')
for iav = 1:nav
    subplot 211
    y = visPrevTrialTC{iav};
    yerr = visPrevTrialTC_sem{iav};
    tt = ((26:length(y))-33).*(1000/frRateHz);
    h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'k');
    h.mainLine.Color = prevtrcol{1,iav};
    h.patch.FaceColor = prevtrpatch{1,iav};
    title('visual')
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim) 
    f = gca;
    f.Box = 'off';
    f.TickDir = 'out';
    hold on
    subplot 212
    y = audPrevTrialTC{iav};
    yerr = audPrevTrialTC_sem{iav};
    tt = ((26:length(y))-33).*(1000/frRateHz);
    h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),'c');
    title('auditory')
    h.mainLine.Color = prevtrcol{2,iav};
    h.patch.FaceColor = prevtrpatch{1,iav};
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    f = gca;
    f.Box = 'off';
    f.TickDir = 'out';
    figYAxis([],'dF/F',responseLim) 
    hold on
end

print([fnout '_tcPrevTrial_lateSupp'],'-dpdf','-fillpage')


%%
respLim = [-0.005 0.02];
tt = ((26:size(avgNormTCEaCycAV,1))-33).*(1000/frRateHz);
figure; 
colmat = ['k','c'];
% norm_tcs_cum = cell(1,nav);
% temp_tcs_cum = cell(1,nav);
for icyc = 1:nCycles
    subplot(2,4,icyc)
    for iav = 1:nav
        
        shadedErrorBar(tt, mean(tcEaCycResponsiveCells{icyc,iav}(26:end,:),2), ste(tcEaCycResponsiveCells{icyc,iav}(26:end,:),2),colmat(:,iav));
        hold on
%         if icyc > 2 & icyc < 7
%             norm_tcs_cum{iav} = cat(3, norm_tcs_cum{iav}, normTCEaCycResponsiveCells{icyc,iav});
%             temp_tcs_cum{iav} = cat(3, temp_tcs_cum{iav}, tcEaCycResponsiveCells{icyc,iav}(:,big_ind));
%         end
    end
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'normalized dF/F',respLim)
    figAxForm([])
    title(['Stimulus #' num2str(icyc)])
end
print([fnout 'eaCycTCResponsiveCells'],'-dpdf','-fillpage')

%% heatmaps of time-courses all anticipation cells
anticipationCells = BLRespCells | lateBLRespCells | lateBLSuppCells;
nAntiCells = sum(anticipationCells);

visTCAllCells = cell2mat(trialTCEaCycEaCell(nCycles,:,visualTrials));
audTCAllCells = cell2mat(trialTCEaCycEaCell(nCycles,:,auditoryTrials));
AVTCAllCells = cell2mat(longTrialTCAV);

visTCAntiCells = visTCAllCells(:,anticipationCells);
audTCAntiCells = audTCAllCells(:,anticipationCells);
AVTCAntiCells = AVTCAllCells(:,anticipationCells);

% [~,sortCellsInd] = sort(mean(cat(2,...
%     respV(anticipationCells,6:nCycles),respA(anticipationCells,6:nCycles)),2));
[~,sortAntiCellsInd] = sort(mean(AVTCAntiCells(end-cycLengthFr*3+1:end,:),1));

[~,sortAllCellsInd] = sort(mean(AVTCAllCells(end-cycLengthFr*3+1:end,:),1));

visTCHeatmap = flipud(visTCAntiCells(:,sortAntiCellsInd)');
audTCHeatmap = flipud(audTCAntiCells(:,sortAntiCellsInd)');
AVTCHeatmap = flipud(AVTCAllCells(:,sortAllCellsInd)');
tt = ((26:size(AVTCHeatmap,2))-33).*(1000/frRateHz);
ttAxisTick = find(ismember(floor(tt),0:500:3000));
ttAxisLabel = floor(tt(ttAxisTick));

VSubAHeatmap = visTCHeatmap - audTCHeatmap;
[~, VSubASortInd] = sort(mean(VSubAHeatmap(:,end-cycLengthFr*3+1:end),2));

respLim = [-0.1 0.1];
figure
colormap(brewermap([],'*RdBu'));
suptitle(sprintf('anticipation responsive cells (%s)',num2str(nAntiCells)));

subplot 221
h = imagesc(visTCHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Visual Trials');

subplot 222
h = imagesc(audTCHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Auditory Trials');

subplot 223
h = imagesc(AVTCHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('All Trials');

subplot 224
% h = imagesc(flipud(VSubAHeatmap(:,26:end)));
h = imagesc(flipud(VSubAHeatmap(VSubASortInd,26:end)));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cells resorted by subtraction magnitude',[])
figAxForm([])
colorbar
caxis(respLim)
title('V - A');

print([fnout 'heatmapsAnticipation'],'-dpdf','-fillpage')

%% example neuron responses
visStimOnTimeMs = ((0:cycLengthFr:cycLengthFr*(nCycles-1))/frRateHz)*1000;
visStimOffTimeMs = visStimOnTimeMs+100;
% exExpt = 10;

exRespCell = 738; %randsample(find(BLRespCells),1);
exLateRespCell = 319; %randsample(find(lateBLRespCells & isTunedAllCells),1);
exLateSuppCell = 1157; %randsample(find(lateBLSuppCells & isTunedAllCells),1);
% exNRCell = randsample(find(~(responsiveCellsAV{exExpt} | ...
%     lateResponseCellInd{exExpt} | lateSuppCellInd{exExpt})),1);

exCellTC = AVTCAllCells(26:end,:);
AVTCAllCellsSem = cell2mat(longTrialTCAV_sem);
exCellTCsem = AVTCAllCellsSem(26:end,:);

respLim = [-0.1 0.25];
setFigParams4Print('portrait')
figure
subplot 311
shadedErrorBar(tt,exCellTC(:,exRespCell),exCellTCsem(:,exRespCell),'k');
hold on
hline(0,'k--')
figYAxis([],'dF/F',respLim)
title(sprintf('respond to 1st stim, #%s',num2str(exRespCell)))
subplot 312
shadedErrorBar(tt,exCellTC(:,exLateRespCell),exCellTCsem(:,exLateRespCell),'k');
hold on
hline(0,'k--')
figYAxis([],'dF/F',respLim)
title(sprintf('respond late, #%s',num2str(exLateRespCell)))
subplot 313
shadedErrorBar(tt,exCellTC(:,exLateSuppCell),exCellTCsem(:,exLateSuppCell),'k');
hold on
hline(0,'k--')
figYAxis([],'dF/F',respLim)
title(sprintf('suppressed late, #%s',num2str(exLateSuppCell)))
% subplot 224
% shadedErrorBar(tt,exCellTC(:,exNRCell),exCellTCsem(:,exNRCell),'k');
% title(num2str(exNRCell))
for i = 1:3
    subplot(3,1,i)
    figXAxis([],'time (ms)',[tt(1) tt(end)],ttAxisLabel,ttAxisLabel)
    f = gca;
    f.Box = 'off';
    f.TickDir = 'out';
    hold on
    vline(visStimOnTimeMs,'k--')
    vline(visStimOffTimeMs,'k--')
end
print([fnout 'exampleCellsTC'],'-dpdf','-fillpage')

% direction tuning
avgOriRespAllCells = cell2mat(avgResponseEaOri');
semOriRespAllCells = cell2mat(semResponseEaOri');
fitOriAllCells = cell2mat(vonMisesFitAllCells);

tunLim = {[-0.05 0.5],[-0.05 0.1],[-0.05 0.1]};
exCells = [exRespCell exLateRespCell exLateSuppCell];
figure
for i = 1:3
    subplot(1,3,i)
    y = avgOriRespAllCells(exCells(i),:);
    yerr = semOriRespAllCells(exCells(i),:);
    yfit = fitOriAllCells(:,exCells(i));
    
    h = errorbar(orientations,y,yerr,'ko');
    h.MarkerFaceColor = [1 1 1];
    hold on
    plot(0:180,yfit,'k-')
    
    figXAxis([],'orientation',[-5 181],orientations,orientations)
    figYAxis([],'dF/F',tunLim{i})
    figAxForm([])
end
print([fnout 'exampleCellsTuning'],'-dpdf','-fillpage')

AVTCAntiCells = AVTCAllCells(:,anticipationCells);
AVTCHeatmapAntiCells = flipud(AVTCAntiCells(:,sortAntiCellsInd)');

exCellsAntiCellInd = find(ismember(find(anticipationCells),exCells));

% marked heatmap
respLim = [-0.1 0.1];
figure
colormap(brewermap([],'*RdBu'))
% sortedExCellInd = find(ismember(flip(sortAllCellsInd),exCells));
sortedExCellInd = find(ismember(flip(sortAntiCellsInd),exCellsAntiCellInd));
% img = AVTCHeatmap(:,26:end);
% negInd = ones(size(img));
% negInd(img < 0) = -1;
% imagesc(log(abs(img)).*negInd);
imagesc(AVTCHeatmapAntiCells(:,26:end))
hold on
hline(sortedExCellInd,'k-')
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('All Trials')
print([fnout '_heatmapExampleCellsMarked'],'-dpdf','-fillpage')
%% tuning all cells
oriTCAllCells = cell(1,4);
for iori = 1:4
    oriTCAllCells{iori} = cell2mat(cellfun(@(x) x(:,:,iori),oriTuningTCs,'unif',0));
end
oriTCRespCells = cellfun(@(x) x(:,anticipationCells),oriTCAllCells,'unif',0);
[~,sortBLRespCellsInd] = sort(mean(...
    AVTCAllCells(end-cycLengthFr*3+1:end,anticipationCells),1));
oriTCRespCellsSort = cellfun(@(x) x(:,fliplr(sortBLRespCellsInd))',oriTCRespCells,'unif',0);

ttori = ([-3 0 3 6])/3;
respLim = [-1 1];

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));
suptitle('drifting grating responses, responsive cells')
for iori = 1:4
    subplot(1,4,iori)
    imagesc(oriTCRespCellsSort{iori}(:,9:18))
    figXAxis([],'time (s)',[],[1:3:9],ttori)
    figYAxis([],'cell #',[])
    f = gca;
    f.Box = 'off';
    f.TickDir = 'out';
    colorbar
    caxis(respLim)
end
print([fnout 'oriTuningHeatmap'],'-dpdf','-fillpage')

avgOriRespAllCells = cell2mat(avgResponseEaOri');
semOriRespAllCells = cell2mat(semResponseEaOri');
fitOriAllCells = cell2mat(vonMisesFitAllCells);

exCell_0 = 593;
exCell_90 = 387;

respLim = [-0.02 0.1];
figure
subplot 121
h = errorbar(orientations,avgOriRespAllCells(exCell_0,:),semOriRespAllCells(exCell_0,:),'ko');
h.MarkerFaceColor = [1 1 1];
hold on
plot(0:180,fitOriAllCells(:,exCell_0))
figXAxis([],'orientation',[-5 181],orientations,orientations)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(num2str(exCell_0))
subplot 122
h = errorbar(orientations,avgOriRespAllCells(exCell_90,:),...
    semOriRespAllCells(exCell_90,:),'ko');
h.MarkerFaceColor = [1 1 1];
hold on
plot(0:180,fitOriAllCells(:,exCell_90))
figXAxis([],'orientation',[-5 181],orientations,orientations)
figYAxis([],'dF/F',respLim)
figAxForm([])
title(num2str(exCell_90))
print([fnout 'oriTuningExCells0_90'],'-dpdf','-fillpage')

%% suppressed cell tuning
antiTuning = avgOriRespAllCells(anticipationCells & isTunedAllCells,:);
[antiTuningPeakResp, antiTuningPeakInd] = max(antiTuning,[],2);
nAntiCell = length(antiTuningPeakResp);

suppTuning = avgOriRespAllCells(lateBLSuppCells & isTunedAllCells,:);
[suppTuningPeakResp, suppTuningPeakInd] = max(suppTuning,[],2);
nSuppCell = length(suppTuningPeakResp);

orthogonalLookup = [3 4 1 2];
antiOffPeakResp = nan(nAntiCell,1);
for icell = 1:nAntiCell
    orthInd = orthogonalLookup(antiTuningPeakInd(icell));
    antiOffPeakResp(icell) = antiTuning(icell,orthInd);
end
    
suppOffPeakResp = nan(nSuppCell,1);
for icell = 1:nSuppCell
    orthInd = orthogonalLookup(suppTuningPeakInd(icell));
    suppOffPeakResp(icell) = suppTuning(icell,orthInd);
end

respLim = [-0.15 0.4];
figure
subplot 121
scatter(antiTuningPeakResp,antiOffPeakResp,'ko')
hold on
scatter(suppTuningPeakResp,suppOffPeakResp,'ro')
plot(respLim,respLim,'k--')
figXAxis([],'pref ori resp (dF/F)',respLim)
figYAxis([],'other ori resp (dF/F)',respLim)
figAxForm([])
hold on
vline(0,'b')
hline(0,'b')
title('red: suppressed cells; blk: all anticipation cells')
subplot 122
h = scatter(mean(antiTuningPeakResp),mean(antiOffPeakResp),'ko');
h.MarkerFaceColor = 'k';
hold on
h = scatter(mean(suppTuningPeakResp),mean(suppOffPeakResp),'ro');
h.MarkerFaceColor = 'r';
plot(respLim,respLim,'k--')
figXAxis([],'pref ori resp (dF/F)',respLim)
figYAxis([],'other ori resp (dF/F)',respLim)
figAxForm([])
hold on
vline(0,'b')
hline(0,'b')
title('mean of cell groups')
print([fnout 'scatterSuppCellsOriResp'],'-dpdf','-fillpage')


binEdgesOri = 0:22.5:180;
nBinnedSuppCells = histcounts(...
    tuningPeakAllCells(lateBLSuppCells & isTunedAllCells)-1,binEdgesOri);
nBinnedAnticipationCells = histcounts(...
    tuningPeakAllCells(anticipationCells & isTunedAllCells)-1,binEdgesOri);

figure;
h = bar(binEdgesOri(1:end-1),nBinnedAnticipationCells);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [1 1 1];
h.BarWidth = 1;
hold on
h = bar(binEdgesOri(1:end-1),nBinnedSuppCells);
h.FaceColor = [0.25 0.25 0.25];
h.EdgeColor = [1 1 1];
h.BarWidth = 1;
figXAxis([],'ori pref (deg)',[])
figYAxis([],'n cells',[])
figAxForm([])
title(sprintf('tuning of all suppressed cells (%s/%s)',...
    num2str(sum(lateBLSuppCells & isTunedAllCells)),num2str(sum(lateBLSuppCells))))
print([fnout 'histTuningSuppCells'],'-dpdf','-fillpage')

%% polor plot of selectivity and tuning
figure;
subplot 121
p = polarplot(deg2rad(tuningPeakAllCells((BLRespCells | lateBLRespCells) & isTunedAllCells)), lateWinSI((BLRespCells | lateBLRespCells)  & isTunedAllCells),'ok');
hold on
p = polarplot(deg2rad(tuningPeakAllCells((BLRespCells | lateBLRespCells)   & lateWinModAV & isTunedAllCells)), lateWinSI((BLRespCells | lateBLRespCells)  & lateWinModAV & isTunedAllCells),'or');
title({'Selectivity index by tuning for all tuned cells- selective cells are red';'first stim responsive'})

subplot 122
p = polarplot(deg2rad(tuningPeakAllCells(lateBLSuppCells & isTunedAllCells)), lateWinSI(lateBLSuppCells & isTunedAllCells),'ok');
hold on
p = polarplot(deg2rad(tuningPeakAllCells(lateBLSuppCells  & lateWinModAV & isTunedAllCells)), lateWinSI(lateBLSuppCells & lateWinModAV & isTunedAllCells),'or');
title({'Selectivity index by tuning for all tuned cells- selective cells are red';'suppressed'})
print([fnout 'selectivityTuningPolarPlot'],'-dpdf','-fillpage')

%% target hit vs miss
targetTCsHits = cell2mat(visTargetTCsHM(:,hitTrials)');
targetTCsHitsBL = targetTCsHits - mean(targetTCsHits(basewinTarget,:));
targetTCsMiss = cell2mat(visTargetTCsHM(:,missTrials)');
targetTCsMissBL = targetTCsMiss - mean(targetTCsMiss(basewinTarget,:));
targetCells = cell2mat(targetResponsiveCells');

targetResponsiveTCsHM = cat(1,{targetTCsHitsBL(:,targetCells)'}...
    ,{targetTCsMissBL(:,targetCells)'});
targetHitResp = mean(targetResponsiveTCsHM{hitTrials}(respwinTarget,:),1)...
    - mean(targetResponsiveTCsHM{1}(basewinTarget,:),1);
targetMissResp = mean(targetResponsiveTCsHM{missTrials}(respwinTarget,:),1)...
    - mean(targetResponsiveTCsHM{1}(basewinTarget,:),1);

[~,sortTargetRespInd] = sort(mean(targetResponsiveTCsHM{missTrials}(:,respwin),2));
sortTargetRespInd = flip(sortTargetRespInd);

timeAfterTargetMs = 400;
timeAfterTargetFr = timeAfterTargetMs*frRateHz/1000;
ttLim = [1 ((frRateHz*(timeAfterTargetMs./1000))+(frRateHz-25))];
ttTarget = ((26:size(targetTCsHits,1))-32).*(1000/frRateHz);
ttTargetAxisTick = find(ismember(round(ttTarget),[-100 0:100:timeAfterTargetMs]));
ttTargetAxisLabel = floor(ttTarget(ttTargetAxisTick));

respLim = [-0.2 0.2];
figure
colormap(brewermap([],'*RdBu'))
suptitle('target responsive cells')
subplot 221
imagesc(targetResponsiveTCsHM{hitTrials}...
    (sortTargetRespInd,26:timeAfterTargetFr+frRateHz))
figXAxis([],'time(ms)',[],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Hit Trials');

subplot 222
imagesc(targetResponsiveTCsHM{missTrials}...
    (sortTargetRespInd,26:timeAfterTargetFr+frRateHz))
figXAxis([],'time(ms)',[],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Miss Trials');

subplot 223
y = mean(targetResponsiveTCsHM{hitTrials}(:,26:timeAfterTargetFr+frRateHz),1);
yerr = ste(targetResponsiveTCsHM{hitTrials}(:,26:timeAfterTargetFr+frRateHz),1);
shadedErrorBar([],y,yerr,'k');
hold on
y = nanmean(targetResponsiveTCsHM{missTrials}(:,26:timeAfterTargetFr+frRateHz),1);
yerr = ste(targetResponsiveTCsHM{missTrials}(:,26:timeAfterTargetFr+frRateHz),1);
shadedErrorBar([],y,yerr,'r');
vline(find(26:timeAfterTargetFr+frRateHz == respwin(1)),'r:')
vline(find(26:timeAfterTargetFr+frRateHz == respwin(end)),'r:')
figXAxis([],'time(ms)',[],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',[])
figAxForm([])
title('target aligned')

subplot 224
for iav = 1:nav
    shadedErrorBar([], mean(tcEaCycResponsiveCells{1,iav}(26:end,:),2),...
        ste(tcEaCycResponsiveCells{1,iav}(26:end,:),2),colmat(:,iav));
    hold on
end
vline(find(26:timeAfterTargetFr+frRateHz == respwinTarget(1)),'r:')
vline(find(26:timeAfterTargetFr+frRateHz == respwinTarget(end)),'r:')
figXAxis([],'time(ms)',[],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',[])
figAxForm([])
title('press aligned, 1st stim responsive cells')

print([fnout '_cmpHitMissTC_target2'],'-dpdf','-fillpage')

%% auROC analysis
auROCTargetVsLastBase = cell2mat(visTargetVsLastBaseAuroc);
uTestTargetVsLastBase = logical(cell2mat(visTargetVsLastBaseUTest));

cellType = {'target resp';'bl resp';'supp'};
cellTypeGroups{1} = targetCells & ~BLRespCells & ~lateBLRespCells & ~lateBLSuppCells;
cellTypeGroups{2} = (BLRespCells | lateBLRespCells) & ~targetCells;
cellTypeGroups{3} = lateBLSuppCells & ~targetCells;

pctLim = [0 0.5];
nCellsCellType = cell(3,1);
auROCBinEdges = 0:0.05:1;
setFigParams4Print('portrait')
figure
subplot 211
ind = cellTypeGroups{1};
nCellsCellType{1} = num2str(sum(ind));
histogram(auROCTargetVsLastBase(ind),auROCBinEdges,...
    'Normalization','probability')
hold on
ind = cellTypeGroups{2} ;
nCellsCellType{2} = num2str(sum(ind));
histogram(auROCTargetVsLastBase(ind),auROCBinEdges,...
    'Normalization','probability')
hold on
ind = cellTypeGroups{3};
nCellsCellType{3} = num2str(sum(ind));
histogram(auROCTargetVsLastBase(ind),auROCBinEdges,...
    'Normalization','probability')
figXAxis([],'auROC',[0 1])
figYAxis([],'probability',pctLim)
vline(0.5,'k--')
legend(cellfun(@(x,y) sprintf('%s, n = %s',x,y),...
    cellType, nCellsCellType,'unif',0),'Location','Northwest')
title({'auROC calculated between target resp (all vis targets) and last base stim resp';...
    'cells are non-overlapping populations'})

subplot 212
ind = (targetCells & ~BLRespCells & ~lateBLRespCells & ~ lateBLSuppCells) & uTestTargetVsLastBase';
nCellsCellType{1} = num2str(sum(ind));
histogram(auROCTargetVsLastBase(ind),auROCBinEdges,...
    'Normalization','probability')
hold on
ind = (BLRespCells | lateBLRespCells) & uTestTargetVsLastBase' & ~targetCells;
nCellsCellType{2} = num2str(sum(ind));
histogram(auROCTargetVsLastBase(ind),auROCBinEdges,...
    'Normalization','probability')
hold on
ind = lateBLSuppCells & ~targetCells & uTestTargetVsLastBase';
nCellsCellType{3} = num2str(sum(ind));
histogram(auROCTargetVsLastBase(ind),auROCBinEdges,...
    'Normalization','probability')
figXAxis([],'auROC',[0 1])
figYAxis([],'probability',pctLim)
vline(0.5,'k--')
title('auROC significantly different than 0.5')
legend(nCellsCellType,'Location','Northwest')
print([fnout '_auROCHistCellGroups'],'-dpdf','-fillpage')

% auROC groups of anticipation & target cells
auROCGroups = {'target > base', 'target < base', 'auROC = 0.5'};
auROCColors = {'r','b','k'};
cellTypeGroupColors = {'b','r','y'};
cellTypeGroupColorsVals = {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]};
auROCGroupInd{1} = uTestTargetVsLastBase' & auROCTargetVsLastBase' > 0.5 ...
    & (anticipationCells | targetCells);
auROCGroupInd{2} = uTestTargetVsLastBase' & auROCTargetVsLastBase' < 0.5 ...
    & (anticipationCells | targetCells);
auROCGroupInd{3} = ~uTestTargetVsLastBase' ...
    & (anticipationCells | targetCells);

% plot target resp of auROC groups
visTargetTCs = cell2mat(visTargetTCsEaExpt');
visTargetBL = visTargetTCs - mean(visTargetTCs(basewinTarget,:),1);
ttTarget = ((1:frRateHz*2)-32)*(1000/frRateHz);
ttTargetAxisTick = -1000:250:1000;
ttTargetAxisLabel = ttTargetAxisTick;
ttRespwinTarget = (respwinTarget([1 end])-32)*(1000/frRateHz);
ttRespwinBL = (respwinTarget([1 end])-32-cycLengthFr)*(1000/frRateHz);

setFigParams4Print('landscape')
figure

subplot 121
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,cellTypeGroups{igroup}),2);
    yerr = ste(visTargetBL(1:frRateHz*2,cellTypeGroups{igroup}),2);
    h = shadedErrorBar(ttTarget,y,yerr,cellTypeGroupColors{igroup});
    h.mainLine.Color = cellTypeGroupColorsVals{igroup};
    h.mainLine.LineWidth = 2;
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',[])
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,cellType,'location','northwest')
title('anticipation/target modulation groups')

subplot 122
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup}),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup}),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',[])
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('auROC groups')
print([fnout '_targetTCGrouped'],'-dpdf','-fillpage')

respLim = [-0.06 0.1];
figure
suptitle('cells group by anticipation groups uses in previous figs')
subplot 131
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & (BLRespCells | lateBLRespCells)),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & (BLRespCells | lateBLRespCells)),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('responsive cells, auROC groups')

subplot 132
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & lateBLSuppCells),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & lateBLSuppCells),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('suppressed cells, auROC groups')

subplot 133
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & targetCells),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & targetCells),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('target cells, auROC groups')
print([fnout '_targetTCAnti&TargetGrouped'],'-dpdf','-fillpage')

%% auROC groups all anticipation cells
figure
legData = [];
auROCGroupsN = cell(1,3);
subplot 121
for igroup = 1:3
    ind = auROCGroupInd{igroup} & anticipationCells;
    y = mean(visTargetBL(1:frRateHz*2,ind),2);
    yerr = ste(visTargetBL(1:frRateHz*2,ind),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
    auROCGroupsN{igroup} = sprintf('%s, n = %s',auROCGroups{igroup},num2str(sum(ind)));
end
hline(0,'k--')
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',[-0.04 0.07])
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroupsN,'location','northwest')
title('anticipation cells, auROC groups')
colorbar

tarLim = [-0.1 0.15];
lastBLAligned = visTargetBL(1:frRateHz*2,:) - ...
    mean(visTargetBL(basewinTarget-cycLengthFr,:),1);
subplot 122
colormap(brewermap([],'*RdBu'))
t = mean(visTargetBL(respwinTarget,anticipationCells),1);
b = mean(lastBLAligned(respwinTarget-cycLengthFr,anticipationCells),1);
h = scatter(b,t,100,auROCTargetVsLastBase(anticipationCells),'.');
hold on
plot(tarLim,tarLim,'k--')
colorbar
caxis([0 1])
figXAxis([],'last base stim',tarLim)
figYAxis([],'target',tarLim)
figAxForm([])
print([fnout '_targetTC&ScatAllAntiCells'],'-dpdf','-fillpage')


%%
respLim = [-0.06 0.1];
figure
suptitle('non-overlapping groups (same as histogram)')
subplot 131
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & cellTypeGroups{2}),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & cellTypeGroups{2}),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('responsive cells, auROC groups')

subplot 132
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & cellTypeGroups{3}),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & cellTypeGroups{3}),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('suppressed cells, auROC groups')

subplot 133
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & cellTypeGroups{1}),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & cellTypeGroups{1}),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('target cells, auROC groups')
print([fnout '_targetTCAnti&TargetGroupedNoOverlap'],'-dpdf','-fillpage')

%% auROC groups, anticipation time-course
visTCnCycles = cell2mat(trialTCEaCycEaCell(nCycles,:,visualTrials));
audTCnCycles = cell2mat(trialTCEaCycEaCell(nCycles,:,auditoryTrials));


tt = ((26:size(visTCnCycles,1))-33).*(1000/frRateHz);
ttAxisTick = -100:500:3000;
ttAxisLabel = ttTargetAxisTick;


respLim = [-0.15 0.25];
setFigParams4Print('landscape')
aurocCellGroupsTCFig = figure;
aurocCellGroupsScatFig = figure;
suptitle('all anticipation cells')
for igroup = 1:3
    figure(aurocCellGroupsTCFig)
    responseLim = [-0.02 0.1];
    subplot(2,3,igroup)
    ind = (BLRespCells | lateBLRespCells) & auROCGroupInd{igroup};
    y = mean(visTCnCycles(26:end,ind),2);
    yerr = ste(visTCnCycles(26:end,ind),2);
    shadedErrorBar(tt,y,yerr,'k');
    hold on
    y = mean(audTCnCycles(26:end,ind),2);
    yerr = ste(audTCnCycles(26:end,ind),2);
    shadedErrorBar(tt,y,yerr,'c');
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim) 
    title(sprintf('%s , n = %s',auROCGroups{igroup},num2str(sum(ind))))
    
    figure(aurocCellGroupsScatFig)
    subplot(2,3,igroup)
    [~,p] = ttest(lateWinRespV(ind),lateWinRespA(ind));
    scatter(lateWinRespV(ind),lateWinRespA(ind),'ko');
    hold on
    plot(respLim,respLim,'k--')
    figXAxis([],'visual (dF/F)',respLim)
    figYAxis([],'auditory (dF/F)',respLim)
    figAxForm([])
    title(sprintf('%s (%s), p = %s',auROCGroups{igroup},num2str(sum(ind)),...
        num2str(p)))
    
    figure(aurocCellGroupsTCFig)
    responseLim = [-0.1 0.02];
    subplot(2,3,igroup+3)
    ind = lateBLSuppCells & auROCGroupInd{igroup};
    y = mean(visTCnCycles(26:end,ind),2);
    yerr = ste(visTCnCycles(26:end,ind),2);
    shadedErrorBar(tt,y,yerr,'k');
    hold on
    y = mean(audTCnCycles(26:end,ind),2);
    yerr = ste(audTCnCycles(26:end,ind),2);
    shadedErrorBar(tt,y,yerr,'c');
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim) 
    title(sprintf('%s , n = %s',auROCGroups{igroup},num2str(sum(ind))))
    if igroup == 3
    print([fnout '_antiTCaurocGroups'],'-dpdf','-fillpage')
    end   
    
    
    figure(aurocCellGroupsScatFig)
    subplot(2,3,igroup+3)
    [~,p] = ttest(lateWinRespV(ind),lateWinRespA(ind));
    scatter(lateWinRespV(ind),lateWinRespA(ind),'ko');
    hold on
    plot(respLim,respLim,'k--')
    figXAxis([],'visual (dF/F)',respLim)
    figYAxis([],'auditory (dF/F)',respLim)
    figAxForm([])
    title(sprintf('%s (%s), p = %s',auROCGroups{igroup},num2str(sum(ind)),...
        num2str(p)))
    if igroup == 3
    print([fnout '_antiScataurocGroups'],'-dpdf','-fillpage')
    end
    
end
% print([fnout '_antiTCAnti&TargetGroupedNoOverlap'],'-dpdf','-fillpage')

%% auROC histograms
figure
subplot 211
histogram(auROCTargetVsLastBase(BLRespCells | lateBLRespCells),20)
hold on
vline(mean(auROCTargetVsLastBase(BLRespCells | lateBLRespCells)),'r-')
vline(0.5,'k--')
figXAxis([],'auROC',[0 1])
title('responsive cells')

subplot 212
histogram(auROCTargetVsLastBase(lateBLSuppCells),20)
hold on
vline(mean(auROCTargetVsLastBase(lateBLSuppCells)),'r-')
vline(0.5,'k--')
figXAxis([],'auROC',[0 1])
title('suppressed cells')

%% auROC target vs first base stim
auROCTargetVsFirstBase = cell2mat(visTargetVsFirstBaseAuroc');
uTestTargetVsFirstBase = logical(cell2mat(visTargetVsFirstBaseUTest'));

pctLim = [0 0.5];
nCellsCellType = cell(3,1);
auROCBinEdges = 0:0.05:1;
setFigParams4Print('portrait')
figure
subplot 211
ind = cellTypeGroups{1};
nCellsCellType{1} = num2str(sum(ind));
histogram(auROCTargetVsFirstBase(ind),auROCBinEdges,...
    'Normalization','probability')
hold on
ind = cellTypeGroups{2} ;
nCellsCellType{2} = num2str(sum(ind));
histogram(auROCTargetVsFirstBase(ind),auROCBinEdges,...
    'Normalization','probability')
hold on
ind = cellTypeGroups{3};
nCellsCellType{3} = num2str(sum(ind));
histogram(auROCTargetVsFirstBase(ind),auROCBinEdges,...
    'Normalization','probability')
figXAxis([],'auROC',[0 1])
figYAxis([],'probability',pctLim)
vline(0.5,'k--')
legend(cellfun(@(x,y) sprintf('%s, n = %s',x,y),...
    cellType, nCellsCellType,'unif',0),'Location','Northwest')
title({'auROC calculated between target resp (all vis targets) and FIRST base stim resp';...
    'cells are non-overlapping populations'})

subplot 212
ind = (targetCells & ~BLRespCells & ~lateBLRespCells & ~ lateBLSuppCells) & uTestTargetVsFirstBase';
nCellsCellType{1} = num2str(sum(ind));
histogram(auROCTargetVsFirstBase(ind),auROCBinEdges,...
    'Normalization','probability')
hold on
ind = (BLRespCells | lateBLRespCells) & uTestTargetVsFirstBase' & ~targetCells;
nCellsCellType{2} = num2str(sum(ind));
histogram(auROCTargetVsFirstBase(ind),auROCBinEdges,...
    'Normalization','probability')
hold on
ind = lateBLSuppCells & ~targetCells & uTestTargetVsFirstBase';
nCellsCellType{3} = num2str(sum(ind));
histogram(auROCTargetVsFirstBase(ind),auROCBinEdges,...
    'Normalization','probability')
figXAxis([],'auROC',[0 1])
figYAxis([],'probability',pctLim)
vline(0.5,'k--')
title('auROC significantly different than 0.5')
legend(nCellsCellType,'Location','Northwest')

% auROC groups of anticipation & target cells
auROCFirstBaseGroupInd{1} = uTestTargetVsFirstBase' & auROCTargetVsFirstBase' > 0.5 ...
    & (anticipationCells | targetCells);
auROCFirstBaseGroupInd{2} = uTestTargetVsFirstBase' & auROCTargetVsFirstBase' < 0.5 ...
    & (anticipationCells | targetCells);
auROCFirstBaseGroupInd{3} = ~uTestTargetVsFirstBase' ...
    & (anticipationCells | targetCells);

setFigParams4Print('landscape')
figure
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCFirstBaseGroupInd{igroup}),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCFirstBaseGroupInd{igroup}),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',[])
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('auROC groups')

respLim = [-0.06 0.1];
figure
suptitle('cells group by anticipation groups uses in previous figs')
subplot 131
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCFirstBaseGroupInd{igroup} & (BLRespCells | lateBLRespCells)),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCFirstBaseGroupInd{igroup} & (BLRespCells | lateBLRespCells)),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('responsive cells, auROC groups')

subplot 132
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCFirstBaseGroupInd{igroup} & lateBLSuppCells),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCFirstBaseGroupInd{igroup} & lateBLSuppCells),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('suppressed cells, auROC groups')

subplot 133
legData = [];
for igroup = 1:3
    y = mean(visTargetBL(1:frRateHz*2,auROCFirstBaseGroupInd{igroup} & targetCells),2);
    yerr = ste(visTargetBL(1:frRateHz*2,auROCFirstBaseGroupInd{igroup} & targetCells),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
vline(ttRespwinBL,'k:')
legend(legData,auROCGroups,'location','northwest')
title('target cells, auROC groups')
print([fnout '_targetTCAnti&TargetGrouped'],'-dpdf','-fillpage')
% auROC groups, anticipation time-course
visTCnCycles = cell2mat(trialTCEaCycEaCell(nCycles,:,visualTrials));
audTCnCycles = cell2mat(trialTCEaCycEaCell(nCycles,:,auditoryTrials));


tt = ((26:size(visTCnCycles,1))-33).*(1000/frRateHz);
ttAxisTick = -100:500:3000;
ttAxisLabel = ttTargetAxisTick;

%% auroc groups first base:target, anticipation responses

setFigParams4Print('landscape')
suptitle('all anticipation cells')
for igroup = 1:3
    figure(aurocCellGroupsTCFig)
    responseLim = [-0.02 0.1];
    subplot(2,3,igroup)
    ind = (BLRespCells | lateBLRespCells) & auROCFirstBaseGroupInd{igroup};
    y = mean(visTCnCycles(26:end,ind),2);
    yerr = ste(visTCnCycles(26:end,ind),2);
    shadedErrorBar(tt,y,yerr,'k');
    hold on
    y = mean(audTCnCycles(26:end,ind),2);
    yerr = ste(audTCnCycles(26:end,ind),2);
    shadedErrorBar(tt,y,yerr,'c');
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim) 
    title(sprintf('%s , n = %s',auROCGroups{igroup},num2str(sum(ind))))
    
    respLim = [-0.15 0.25];
    figure(aurocCellGroupsScatFig)
    subplot(2,3,igroup)
    [~,p] = ttest(lateWinRespV(ind),lateWinRespA(ind));
    scatter(lateWinRespV(ind),lateWinRespA(ind),'ko');
    hold on
    plot(respLim,respLim,'k--')
    figXAxis([],'visual (dF/F)',respLim)
    figYAxis([],'auditory (dF/F)',respLim)
    figAxForm([])
    title(sprintf('%s (%s), p = %s',auROCGroups{igroup},num2str(sum(ind)),...
        num2str(p)))
    
    figure(aurocCellGroupsTCFig)
    responseLim = [-0.1 0.02];
    subplot(2,3,igroup+3)
    ind = lateBLSuppCells & auROCFirstBaseGroupInd{igroup};
    y = mean(visTCnCycles(26:end,ind),2);
    yerr = ste(visTCnCycles(26:end,ind),2);
    shadedErrorBar(tt,y,yerr,'k');
    hold on
    y = mean(audTCnCycles(26:end,ind),2);
    yerr = ste(audTCnCycles(26:end,ind),2);
    shadedErrorBar(tt,y,yerr,'c');
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim) 
    title(sprintf('%s , n = %s',auROCGroups{igroup},num2str(sum(ind))))
    if igroup == 3
    print([fnout '_antiScataurocGroups'],'-dpdf','-fillpage')
    end   
    
    
    figure(aurocCellGroupsScatFig)
    subplot(2,3,igroup+3)
    [~,p] = ttest(lateWinRespV(ind),lateWinRespA(ind));
    scatter(lateWinRespV(ind),lateWinRespA(ind),'ko');
    hold on
    plot(respLim,respLim,'k--')
    figXAxis([],'visual (dF/F)',respLim)
    figYAxis([],'auditory (dF/F)',respLim)
    figAxForm([])
    title(sprintf('%s (%s), p = %s',auROCGroups{igroup},num2str(sum(ind)),...
        num2str(p)))
    if igroup == 3
    print([fnout '_antiTCaurocGroups'],'-dpdf','-fillpage')
    end   
end

% %% auROC vs orientation
% auROCZeroTargetVsLastBase = auROCTargetVsLastBase-0.5;
% figure;
% subplot 121
% p = polarplot(deg2rad(tuningPeakAllCells((BLRespCells | lateBLRespCells) & isTunedAllCells)),...
%     auROCZeroTargetVsLastBase((BLRespCells | lateBLRespCells) & isTunedAllCells),'ok');
% hold on
% p = polarplot(deg2rad(tuningPeakAllCells((BLRespCells | lateBLRespCells)  & uTestTargetVsLastBase' & isTunedAllCells)),...
%     auROCZeroTargetVsLastBase((BLRespCells | lateBLRespCells) & uTestTargetVsLastBase' & isTunedAllCells),'or');
% title({'auROC by tuning for all tuned cells- selective cells are red';'responsive'})
% 
% subplot 122
% p = polarplot(deg2rad(tuningPeakAllCells(lateBLSuppCells & isTunedAllCells)), ...
%     auROCZeroTargetVsLastBase(lateBLSuppCells & isTunedAllCells),'ok');
% hold on
% p = polarplot(deg2rad(tuningPeakAllCells(lateBLSuppCells  & uTestTargetVsLastBase' & isTunedAllCells)),...
%     auROCZeroTargetVsLastBase(lateBLSuppCells & uTestTargetVsLastBase' & isTunedAllCells),'or');
% title({'auROC by tuning for all tuned cells- selective cells are red';'suppressed'})
% 
% % print([fnout 'selectivityTuningPolarPlot'],'-dpdf','-fillpage')

%% cells selected by auditory trials
audResponsiveCells = logical(cell2mat(responsiveCellsA'));
lateAudRespCells = logical(cell2mat(lateResponsiveCellsA') & ~audResponsiveCells);
lateAudSuppCells = logical(cell2mat(lateSuppressedCellsA') & ~audResponsiveCells);

setFigParams4Print('landscape')
figure
suptitle('auditory responsive cells')
responseLim = [-0.02 0.065];
subplot 311
ind = audResponsiveCells;
y = mean(visTCnCycles(26:end,ind),2);
yerr = ste(visTCnCycles(26:end,ind),2);
shadedErrorBar(tt,y,yerr,'k');
hold on
y = mean(audTCnCycles(26:end,ind),2);
yerr = ste(audTCnCycles(26:end,ind),2);
shadedErrorBar(tt,y,yerr,'c');
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',responseLim) 
title(sprintf('1st stim responsive, n = %s',num2str(sum(ind))))

subplot 312
ind = lateAudRespCells;
y = mean(visTCnCycles(26:end,ind),2);
yerr = ste(visTCnCycles(26:end,ind),2);
shadedErrorBar(tt,y,yerr,'k');
hold on
y = mean(audTCnCycles(26:end,ind),2);
yerr = ste(audTCnCycles(26:end,ind),2);
shadedErrorBar(tt,y,yerr,'c');
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',responseLim) 
title(sprintf('late responsive, n = %s',auROCGroups{igroup},num2str(sum(ind))))

responseLim = [-0.065 0.02];
subplot 313
ind = lateAudSuppCells;
y = mean(visTCnCycles(26:end,ind),2);
yerr = ste(visTCnCycles(26:end,ind),2);
shadedErrorBar(tt,y,yerr,'k');
hold on
y = mean(audTCnCycles(26:end,ind),2);
yerr = ste(audTCnCycles(26:end,ind),2);
shadedErrorBar(tt,y,yerr,'c');
figXAxis([],'time (ms)',[tt(1) tt(end)])
figYAxis([],'dF/F',responseLim) 
title(sprintf('late suppressed, n = %s',auROCGroups{igroup},num2str(sum(ind))))


print([fnout '_antiTCAudResponsiveCellSelection'],'-dpdf','-fillpage')

%% iti F for suppressed cells

itiF = cell2mat(itiFEaExpt');
lateWinF = cell2mat(lateWinFEaExpt');

itiFCellGroup(1) = mean(itiF(BLRespCells | lateBLRespCells));
itiFCellGroup(2) = mean(itiF(lateBLSuppCells));
itiFCellGroup(3) = mean(itiF(anticipationCells));

itiFCellGroup_sem(1) = ste(itiF(BLRespCells | lateBLRespCells),2);
itiFCellGroup_sem(2) = ste(itiF(lateBLSuppCells),2);
itiFCellGroup_sem(3) = ste(itiF(anticipationCells),2);

figure;
h = bar(1:3,itiFCellGroup);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [1 1 1];
hold on
h = errorbar(1:3,itiFCellGroup,itiFCellGroup_sem,'k.');
figXAxis([],'',[], 1:3,{'responsive','suppressed','all'})
figYAxis([],'iti F',[])
figAxForm([]);
title('F during inter-trial interval')
print([fnout '_itiFCellGroups'],'-dpdf','-fillpage')

%%  selectivity of neuron populations by experiment
rCells = cellfun(@(x,y,z) (x & y) | z, responsiveCellsV,...
    responseCutoffMet,lateResponseCellInd,'unif',0);
sCells = lateSuppCellInd;

lateWinSIEaExpt_rCells = cellfun(@(x,y) x(y),lateWinSIEaExpt,rCells,'unif',0);
lateWinSIEaExpt_sCells = cellfun(@(x,y) x(y),lateWinSIEaExpt,sCells,'unif',0);

[~,p] = 

leg = [];
siLim = [-5 5];
nLim = [0 60];
figure
colormap('cool')
subplot 121
h = plot(ones(1,nexp),cellfun(@mean,lateWinSIEaExpt_rCells),'ko');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h = errorbar(1,mean(cellfun(@mean,lateWinSIEaExpt_rCells)),...
    ste(cellfun(@mean,lateWinSIEaExpt_rCells),1),'ko');
h.MarkerFaceColor = 'k';
h = plot(ones(1,nexp)*2,cellfun(@mean,lateWinSIEaExpt_sCells),'ko');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
leg(1) = h;
h = errorbar(2,nanmean(cellfun(@mean,lateWinSIEaExpt_sCells)),...
    ste(cellfun(@mean,lateWinSIEaExpt_sCells),1),'ko');
h.MarkerFaceColor = 'k';
leg(2) = h;
plot(1:2,...
    [cellfun(@mean,lateWinSIEaExpt_rCells) cellfun(@mean,lateWinSIEaExpt_sCells)],...
    'k-')
figXAxis([],'cell group',[0 3],1:2,{'resp','supp'})
figYAxis([],'mean SI',siLim)
hline(0,'k--')
legend(leg,{'ea expt','across expt'})
fig = gca;
fig.Box = 'off';
fig.TickDir = 'out';

subplot 122
h = scatter(ones(1,nexp),cellfun(@mean,lateWinSIEaExpt_rCells),200,...
    cellfun(@sum,rCells),'.');
hold on
h = scatter(ones(1,nexp)*2,cellfun(@mean,lateWinSIEaExpt_sCells),200,...
    cellfun(@sum,sCells),'.');
colorbar
caxis(nLim)
figXAxis([],'cell group',[0 3],1:2,{'resp','supp'})
figYAxis([],'mean SI',siLim)
hline(0,'k--')
hline(0,'k--')
legend(leg,{'ea expt','across expt'})
fig = gca;
fig.Box = 'off';
fig.TickDir = 'out';
title('color is n cells')
print([fnout '_SIeaExptCellGroups'],'-dpdf','-fillpage')