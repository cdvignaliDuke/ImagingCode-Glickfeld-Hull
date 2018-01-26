clear all
close all
ds = 'FSAV_V1_GAD';
cellsOrDendrites = 1;
doRawData = 0;
%%
rc = behavConstsAV;
dataGroup = ds;
eval(dataGroup)
titleStr = ds(6:end);
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];

if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary_cells' ds(5:end) '.mat']));    
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary_dendrites' ds(5:end) '.mat']));    
end        

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
nav = 2;
hitTrials = 1;
missTrials = 2;
nout = 2;
cycLengthFr = mouse(1).expt(1).info.cyc_time;
frRateHz = expt(1).frame_rate;
cellLabelName = expt(1).redChannelLabel;
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

%%
szEaExpt = nan(1,nexp);
exptDepth = nan(1,nexp);
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
        dirTuningMworks = loadMworksFile(expt(ind).SubNum,dt,expt(ind).dirtuning_time);
        szEaExpt(exptN) = dirTuningMworks.gratingDiameterDeg;
        if strcmp(ds(end-4:end),'naive')
            exptIndicator{exptN} = strjoin(expt(ind).indicator);
        end
        exptDepth(exptN) = expt(ind).z;
    end
end
%% responses
basewin = 1:34;
basewin_0 = 32:34;
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
exptSN = nan(1,nexp);
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
nCatchTrials = zeros(nexp,1);
labeledCellsEaExpt = cell(nexp,1);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end
            if mouse(imouse).expt(iexp).info.isCatch
                catchHit = sum(cellfun(@length,mouse(imouse).expt(iexp).align(3)...
                    .av(visualTrials).outcome(3).tcyc));
                catchMiss = sum(cellfun(@length,mouse(imouse).expt(iexp).align(3)...
                    .av(visualTrials).outcome(4).tcyc));
                nCatchTrials(exptN) = catchHit+catchMiss;
            end
        exptName{exptN} = [mouse(imouse).expt(iexp).mouse_name '-' ...
            mouse(imouse).expt(iexp).date];
        exptSN(exptN) = str2num(mouse(imouse).expt(iexp).mouse_name);
        labeledCellsEaExpt{exptN} = mouse(imouse).expt(iexp).cells.isLabeledCell;
        dAV_cycle1 = cat(3,mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{1},...
            mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{1});
        dV_cycle1 = mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{1};
        dA_cycle1 = mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{1};
        responsiveCellsV{exptN} = ttest(squeeze(mean(dAV_cycle1(respwin,:,:),1)),...
            squeeze(mean(dAV_cycle1(basewin,:,:),1)),'dim',2,'tail','right','alpha',cellGroupsAlpha)';
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
% %         pmtResponsive{exptN} = ttest(squeeze(mean(dV_cycle1(respwin-4,:,:),1)),...
% %             squeeze(mean(dV_cycle1(basewin-4,:,:),1)),'dim',2,'tail','left')';
        
        
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
        lateResponsiveCellsA{exptN} = ttest(squeeze(mean(dA_nCycles(respwin(1):end,:,:),1)),...
            squeeze(mean(dA_nCycles(basewin,:,:),1)),'dim',2,'tail','right','alpha',cellGroupsAlpha)';
        lateSuppressedCellsA{exptN} = ttest(squeeze(mean(dA_nCycles(respwin(1):end,:,:),1)),...
            squeeze(mean(dA_nCycles(basewin,:,:),1)),'dim',2,'tail','left','alpha',cellGroupsAlpha)';
        
        
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
                baselineAllTrials = mean(tcThisCycAllOutcome(basewin_0+extraFrames,:,:),1);
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
vline((basewin_0([1 end])-33).*(1000/frRateHz),'k:')
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
vline((basewin_0([1 end])-33).*(1000/frRateHz),'k:')
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

%% expt cell counts table
nRespCellEaExpt = cellfun(@(x,y,z) sum((x & y) | z),responsiveCellsV,...
    responseCutoffMet,lateResponseCellInd);
nSuppCellEaExpt = cellfun(@sum,lateSuppCellInd);
nAllSuppCell = cellfun(@sum,lateSuppressedCellsAV);
nTargetRespCellEaExpt = cellfun(@sum,targetResponsiveCells);
exptCellCountTable = table(nCellEaExpt,nRespCellEaExpt,nSuppCellEaExpt,nAllSuppCell,nTargetRespCellEaExpt,szEaExpt',exptDepth');
exptCellCountTable.Properties.RowNames = exptName';
exptCellCountTable.Properties.VariableNames([6,7]) = {'stimSize','exptDepth'};


%% cell groups
BLRespCells = logical(cell2mat(responsiveCellsV')) & cell2mat(responseCutoffMet');
TarRespCells = logical(cell2mat(targetResponsiveCells'));
lateBLRespCells = cell2mat(lateResponseCellInd');
lateBLSuppCells = cell2mat(lateSuppCellInd');
labeledCells = cell2mat(labeledCellsEaExpt);
%%
tc = squeeze(trialTCEaCycEaCell(nCycles,:,:));
tcAllCells = cell(1,nav);
for iav = 1:nav
    tcAllCells{iav} = cell2mat(tc(:,iav)');
end
%% heatmaps of labeled and unlabeled cells
responseLim = [-0.02 0.1];
ind = labeledCells' & BLRespCells;
nLabeledCells = sum(ind);

visTCAllCells = cell2mat(trialTCEaCycEaCell(nCycles,:,visualTrials));
audTCAllCells = cell2mat(trialTCEaCycEaCell(nCycles,:,auditoryTrials));
AVTCAllCells = cell2mat(longTrialTCAV);

visTCLabeledCells = visTCAllCells(:,ind);
audTCLabeledCells = audTCAllCells(:,ind);
AVTCLabeledCells = AVTCAllCells(:,ind);

[~,sortLabeledCellsInd] = sort(mean(AVTCLabeledCells(end-cycLengthFr*3+1:end,:),1));

visTCHeatmap = flipud(visTCLabeledCells(:,sortLabeledCellsInd)');
audTCHeatmap = flipud(audTCLabeledCells(:,sortLabeledCellsInd)');
AVTCHeatmap = flipud(AVTCLabeledCells(:,sortLabeledCellsInd)');
tt = ((26:size(AVTCHeatmap,2))-33).*(1000/frRateHz);
ttAxisTick = find(ismember(floor(tt),0:500:3000));
ttAxisLabel = floor(tt(ttAxisTick));

respLim = [-0.1 0.1];
figure
colormap(brewermap([],'*RdBu'));

subplot 231
h = imagesc(visTCHeatmap(:,26:end));
figXAxis([],'Time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Visual Trials');

subplot 232
h = imagesc(audTCHeatmap(:,26:end));
figXAxis([],'Time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Auditory Trials');

subplot 233
y = mean(tcAllCells{visualTrials}(:,ind),2);
yerr = ste(tcAllCells{visualTrials}(:,ind),2);
h = shadedErrorBar([],y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcAllCells{auditoryTrials}(:,ind),2);
yerr = ste(tcAllCells{auditoryTrials}(:,ind),2);
h = shadedErrorBar([],y(26:length(y)),yerr(26:length(y)),'c');
figXAxis([],'Time (ms)',[1 length(26:length(y))],ttAxisTick,ttAxisLabel)
figYAxis([],'dF/F',responseLim)
figAxForm([])
colorbar
title('Blk: Visual, Cyan: Auditory');

ind = ~labeledCells' & BLRespCells;
nNonLabeledCells = sum(ind);
suptitle(sprintf('Top: %s+ cells (%s), Bottom: %s- cells (%s)',cellLabelName,...
    num2str(nLabeledCells),cellLabelName,num2str(nNonLabeledCells)));

visTCNonLabeledCells = visTCAllCells(:,ind);
audTCNonLabeledCells = audTCAllCells(:,ind);
AVTCNonLabeledCells = AVTCAllCells(:,ind);

[~,sortNonLabeledCellsInd] = sort(mean(AVTCNonLabeledCells(end-cycLengthFr*3+1:end,:),1));

visTCHeatmap = flipud(visTCNonLabeledCells(:,sortNonLabeledCellsInd)');
audTCHeatmap = flipud(audTCNonLabeledCells(:,sortNonLabeledCellsInd)');
AVTCHeatmap = flipud(AVTCNonLabeledCells(:,sortNonLabeledCellsInd)');

subplot 234
h = imagesc(visTCHeatmap(:,26:end));
figXAxis([],'Time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Visual Trials');

subplot 235
h = imagesc(audTCHeatmap(:,26:end));
figXAxis([],'Time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Auditory Trials');

subplot 236
y = mean(tcAllCells{visualTrials}(:,ind),2);
yerr = ste(tcAllCells{visualTrials}(:,ind),2);
h = shadedErrorBar([],y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcAllCells{auditoryTrials}(:,ind),2);
yerr = ste(tcAllCells{auditoryTrials}(:,ind),2);
h = shadedErrorBar([],y(26:length(y)),yerr(26:length(y)),'c');
figXAxis([],'Time (ms)',[1 length(26:length(y))],ttAxisTick,ttAxisLabel)
figYAxis([],'dF/F',responseLim)
figAxForm([])
colorbar
title('Blk: Visual, Cyan: Auditory');

print([fnout '_HM_LabeledCells_Responsive'],'-dpdf','-fillpage')

ind = labeledCells' & lateBLSuppCells;
nLabeledCells = sum(ind);

visTCAllCells = cell2mat(trialTCEaCycEaCell(nCycles,:,visualTrials));
audTCAllCells = cell2mat(trialTCEaCycEaCell(nCycles,:,auditoryTrials));
AVTCAllCells = cell2mat(longTrialTCAV);

visTCLabeledCells = visTCAllCells(:,ind);
audTCLabeledCells = audTCAllCells(:,ind);
AVTCLabeledCells = AVTCAllCells(:,ind);

[~,sortLabeledCellsInd] = sort(mean(AVTCLabeledCells(end-cycLengthFr*3+1:end,:),1));

visTCHeatmap = flipud(visTCLabeledCells(:,sortLabeledCellsInd)');
audTCHeatmap = flipud(audTCLabeledCells(:,sortLabeledCellsInd)');
AVTCHeatmap = flipud(AVTCLabeledCells(:,sortLabeledCellsInd)');
tt = ((26:size(AVTCHeatmap,2))-33).*(1000/frRateHz);
ttAxisTick = find(ismember(floor(tt),0:500:3000));
ttAxisLabel = floor(tt(ttAxisTick));

respLim = [-0.1 0.1];
figure
colormap(brewermap([],'*RdBu'));

subplot 231
h = imagesc(visTCHeatmap(:,26:end));
figXAxis([],'Time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Visual Trials');

subplot 232
h = imagesc(audTCHeatmap(:,26:end));
figXAxis([],'Time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Auditory Trials');

subplot 233
y = mean(tcAllCells{visualTrials}(:,ind),2);
yerr = ste(tcAllCells{visualTrials}(:,ind),2);
h = shadedErrorBar([],y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcAllCells{auditoryTrials}(:,ind),2);
yerr = ste(tcAllCells{auditoryTrials}(:,ind),2);
h = shadedErrorBar([],y(26:length(y)),yerr(26:length(y)),'c');
figXAxis([],'Time (ms)',[1 length(26:length(y))],ttAxisTick,ttAxisLabel)
figYAxis([],'dF/F',responseLim)
figAxForm([])
colorbar
title('Blk: Visual, Cyan: Auditory');

ind = ~labeledCells' & lateBLSuppCells;
nNonLabeledCells = sum(ind);
suptitle(sprintf('Top: %s+ cells (%s), Bottom: %s- cells (%s)',cellLabelName,...
    num2str(nLabeledCells),cellLabelName,num2str(nNonLabeledCells)));

visTCNonLabeledCells = visTCAllCells(:,ind);
audTCNonLabeledCells = audTCAllCells(:,ind);
AVTCNonLabeledCells = AVTCAllCells(:,ind);

[~,sortNonLabeledCellsInd] = sort(mean(AVTCNonLabeledCells(end-cycLengthFr*3+1:end,:),1));

visTCHeatmap = flipud(visTCNonLabeledCells(:,sortNonLabeledCellsInd)');
audTCHeatmap = flipud(audTCNonLabeledCells(:,sortNonLabeledCellsInd)');
AVTCHeatmap = flipud(AVTCNonLabeledCells(:,sortNonLabeledCellsInd)');

subplot 234
h = imagesc(visTCHeatmap(:,26:end));
figXAxis([],'Time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Visual Trials');

subplot 235
h = imagesc(audTCHeatmap(:,26:end));
figXAxis([],'Time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Auditory Trials');

subplot 236
y = mean(tcAllCells{visualTrials}(:,ind),2);
yerr = ste(tcAllCells{visualTrials}(:,ind),2);
h = shadedErrorBar([],y(26:length(y)),yerr(26:length(y)),'k');
hold on
y = mean(tcAllCells{auditoryTrials}(:,ind),2);
yerr = ste(tcAllCells{auditoryTrials}(:,ind),2);
h = shadedErrorBar([],y(26:length(y)),yerr(26:length(y)),'c');
figXAxis([],'Time (ms)',[1 length(26:length(y))],ttAxisTick,ttAxisLabel)
figYAxis([],'dF/F',responseLim)
figAxForm([])
colorbar
title('Blk: Visual, Cyan: Auditory');

print([fnout '_HM_LabeledCells_Supp'],'-dpdf','-fillpage')
%% cell group TCs
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

% visTCAllCells = cell2mat(trialTCEaCycEaCell(nCycles,:,visualTrials));
% audTCAllCells = cell2mat(trialTCEaCycEaCell(nCycles,:,auditoryTrials));
% AVTCAllCells = cell2mat(longTrialTCAV);

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

%%
figure;
colormap(brewermap([],'*RdBu'));
subplot 221
h = histogram(lateRespSI,binEdgesSI);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [1 1 1];
h = vline(0,'k-');
h.LineWidth = 2;
h = vline(mean(lateRespSI),'r-');
h.LineWidth = 2;
hold on
h = histogram(lateRespSI(logical(lateRespModAV)),binEdgesSI);
h.FaceColor = 'k';
h.EdgeColor = [1 1 1];
h = vline(mean(lateRespSI(logical(lateRespModAV))),'m-');
h.LineWidth = 2;
figXAxis([],'late resp selectivity index (using cycles 5-8)',siLim)
figYAxis([],'n ncells',[])
figAxForm([])
title({sprintf('significantly modulated late (%s/%s)'...
    ,num2str(sum(lateRespModAV)),num2str(length(lateRespSI))),...
    'red: avg SI all cels, pink: avg SI signif cells'})

% visModHeatmap = flipud(AVTCAllCells(:,sortAllCellsInd(lateRespModAV & lateRespSI > 0))');
% audModHeatmap = flipud(AVTCAllCells(:,sortAllCellsInd(lateRespModAV & lateRespSI < 0))');

visModTC = AVTCAllCells(:,lateRespModAV & lateRespSI > 0);
[~,visModSortInd] = sort(mean(visModTC(end-lateWinFrames:end,:),1));
visModHeatmap = flipud(visModTC(:,visModSortInd)');
audModTC = AVTCAllCells(:,lateRespModAV & lateRespSI < 0);
[~,audModSortInd] = sort(mean(audModTC(end-lateWinFrames:end,:),1));
audModHeatmap = flipud(audModTC(:,audModSortInd)');

subplot 222
h = imagesc(visModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Vis Selective Cells (+SI), All Trials');
subplot 224
h = imagesc(audModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('Aud Selective Cells (-SI), All Trials');

print([fnout '_modCellsSelectivity&Heatmaps'],'-dpdf','-fillpage')

vis_visModTC = visTCAllCells(:,lateRespModAV & lateRespSI > 0);
aud_visModTC = audTCAllCells(:,lateRespModAV & lateRespSI > 0);

vis_visModHeatmap = flipud(vis_visModTC(:,visModSortInd)');
aud_visModHeatmap = flipud(aud_visModTC(:,visModSortInd)');

vis_audModTC = visTCAllCells(:,lateRespModAV & lateRespSI < 0);
aud_audModTC = audTCAllCells(:,lateRespModAV & lateRespSI < 0);

vis_audModHeatmap = flipud(vis_audModTC(:,audModSortInd)');
aud_audModHeatmap = flipud(aud_audModTC(:,audModSortInd)');

figure
colormap(brewermap([],'*RdBu'));
subplot 221
h = imagesc(vis_visModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('+SI Cells, Visual Trials');
subplot 222
h = imagesc(aud_visModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('+SI Cells, Auditory Trials');
subplot 223
h = imagesc(vis_audModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('-SI Cells, Visual Trials');
subplot 224
h = imagesc(aud_audModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('-SI Cells, Auditory Trials');
print([fnout '_modCellsSelectivity&Heatmaps_visVSaudTrials'],'-dpdf','-fillpage')

[~,visModSortInd] = sort(mean(AVTCAllCells(end-lateWinFrames:end,lateRespSI > 0),1));
[~,audModSortInd] = sort(mean(AVTCAllCells(end-lateWinFrames:end,lateRespSI < 0),1));

vis_visModTC = visTCAllCells(:,lateRespSI > 0);
aud_visModTC = audTCAllCells(:,lateRespSI > 0);

vis_visModHeatmap = flipud(vis_visModTC(:,visModSortInd)');
aud_visModHeatmap = flipud(aud_visModTC(:,visModSortInd)');

vis_audModTC = visTCAllCells(:,lateRespSI < 0);
aud_audModTC = audTCAllCells(:,lateRespSI < 0);

vis_audModHeatmap = flipud(vis_audModTC(:,audModSortInd)');
aud_audModHeatmap = flipud(aud_audModTC(:,audModSortInd)');

figure
colormap(brewermap([],'*RdBu'));
subplot 221
h = imagesc(vis_visModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('+SI Cells, Visual Trials');
subplot 222
h = imagesc(aud_visModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('+SI Cells, Auditory Trials');
subplot 224
h = imagesc(vis_audModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('-SI Cells, Visual Trials');
subplot 223
h = imagesc(aud_audModHeatmap(:,26:end));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colorbar
caxis(respLim)
title('-SI Cells, Auditory Trials');
print([fnout '_allCellsSelectivity&Heatmaps_visVSaudTrials'],'-dpdf','-fillpage')

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
fitOriAllCells = cell2mat(oriVonMisesFitAllCells);

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
%%
targetCells = cell2mat(targetResponsiveCells');
%% target hit vs miss
targetTCsHits = cell2mat(visTargetTCsHM(:,hitTrials)');
targetTCsHitsBL = targetTCsHits - mean(targetTCsHits(basewinTarget,:));
targetTCsMiss = cell2mat(visTargetTCsHM(:,missTrials)');
targetTCsMissBL = targetTCsMiss - mean(targetTCsMiss(basewinTarget,:));

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
vline(ttRespwinBL,'k:')
    yerr = ste(visTargetBL(1:frRateHz*2,auROCGroupInd{igroup} & cellTypeGroups{1}),2);
    h = shadedErrorBar(ttTarget,y,yerr,auROCColors{igroup});
    legData(igroup) = h.mainLine;
    hold on
end
figXAxis([],'time(ms)',[-1000 1000],ttTargetAxisTick,ttTargetAxisLabel)
figYAxis([],'dF/F',respLim)
figAxForm([])
vline(ttRespwinTarget,'k:')
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

setFigParams4Print('portrait')
binEdgesSI = -10:0.5:10;
siLim = [-10 10];
figure
subplot 211
ind = audResponsiveCells | lateAudRespCells;
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
ind = lateAudSuppCells;
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
% print([fnout '_SICellGroupsLateWin'],'-dpdf','-fillpage')
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

%% first response for aud and vis selected cells
responseCutoffA = cellfun(@(x) mean(x(respwin,:),1) > 0.002,...
    squeeze(tcEaCellEaCycAllOutcome(1,:,auditoryTrials)),'unif',0);
respCellA = cell2mat(cellfun(@(x,y) x & y, responsiveCellsA',responseCutoffA,'unif',0));
respCellV = cell2mat(cellfun(@(x,y) x & y, responsiveCellsV,responseCutoffMet,'unif',0)');

respLim = [-0.005 0.04];
figure;
suptitle(sprintf('cells selected from (%s cells overlap) ...',num2str(sum(respCellA & respCellV))))
for iplot = 1:4
    visTC = cell2mat(squeeze(tcEaCellEaCycAllOutcome(iplot,:,visualTrials)));
    audTC = cell2mat(squeeze(tcEaCellEaCycAllOutcome(iplot,:,auditoryTrials)));
    tt = ((26:size(visTC,1))-33).*(1000/frRateHz);
    subplot(2,4,iplot)
    y = mean(visTC(26:end,respCellV),2);
    yerr = ste(visTC(26:end,respCellV),2);
    shadedErrorBar(tt,y,yerr,'k');
    hold on
    y = mean(audTC(26:end,respCellV),2);
    yerr = ste(audTC(26:end,respCellV),2);
    shadedErrorBar(tt,y,yerr,'c');
    if iplot == 1
        title(sprintf('visual trials, n = %s, stim %s',num2str(sum(respCellV)),...
            num2str(iplot)))
    else
        title(sprintf('stim %s',num2str(iplot)))
    end
    figXAxis([],'time(ms)',[],[-200:200:200],[-200:200:200])
    figYAxis([],'dF/F',respLim)
    figAxForm([],0)
    subplot(2,4,iplot+4)
    y = mean(visTC(26:end,respCellA),2);
    yerr = ste(visTC(26:end,respCellA),2);
    shadedErrorBar(tt,y,yerr,'k');
    hold on
    y = mean(audTC(26:end,respCellA),2);
    yerr = ste(audTC(26:end,respCellA),2);
    shadedErrorBar(tt,y,yerr,'c');
    if iplot == 1
        title(sprintf('auditory trials, n = %s, stim %s',num2str(sum(respCellA)),...
            num2str(iplot)))
    else
        title(sprintf('stim %s',num2str(iplot)))
    end
    figXAxis([],'time(ms)',[],[-200:200:200],[-200:200:200])
    figYAxis([],'dF/F',respLim)
end

print([fnout '_tcFirstStimAud&VisSelectCells'],'-dpdf','-fillpage')