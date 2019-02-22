clear all
close all
ds = 'FSAV_attentionV1';
cellsOrDendrites = 1;
%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        [titleStr '_anticipation_']); 
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        [titleStr '_den_anticipation_']); 
end

if strcmp(ds,'FSAV_attentionV1')
    load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        'V1_100ms_naive_anticipation_imgStats'))
    imgStats_naive = imgStats;
    clear imgStats
end

%%

nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 8;
% lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;
minTargetRT = (nVisDelayFr_target+respwin_target(1)-nBaselineFr)./frameRateHz*1000;

oriBinSize = 45;
orientations = 0:oriBinSize:(180-oriBinSize);
oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
nOri = length(orientations);

nMovWin = 15;
movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelFr_target = 30:(30+nMovWin-1);
% movWinLabelMs = 

trOutType = {'h';'m';'fa';'cr'};
trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
%% pool experiment data

antiDataExpt = struct;
oriTuningExpt = struct;
targetDataExpt = struct;
decodeDataExpt = struct;
decodeDataExpt.av(visualTrials).name = 'Visual';
decodeDataExpt.av(auditoryTrials).name = 'Auditory';
nTargets = 2; %sum(unique(targetInd) > 1);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if imouse == 1 && iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        
        d = mouse(imouse).expt(iexp);
        
        cycLengthFr = d.info.cycTimeFrames;
        nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);
        
        maxCycles = max(cat(2,d.av(visualTrials).align(alignStart).nCycles,...
            d.av(auditoryTrials).align(alignStart).nCycles));
        
        tc_AV = [];
%         cycTC = cell(2,nCycles);
        cycTC = cell(2,maxCycles);
        longTC = cell(2,1);
        for iav = 1:2
            dd = d.av(iav);
            trOut = [];
            trStim = [];
            trResp = [];
            trTC = [];
            trResp_movWin = [];
            for ialign = 2:4
                ddd = dd.align(ialign);
%                 if iav 
%                     fprintf('%s/%s auditory trials < 200 ms\n', ...
%                         num2str(sum(ddd.reactTime < 200)),num2str(length(ddd.reactTime)))
%                 end
                tc = ddd.respTC;
                if ialign == 4
                    tc_temp = circshift(tc,-1,1);
                    tc_temp(end,:,:) = nan;
                    trTC = cat(3,trTC,tc);
                    trResp = cat(2,trResp,squeeze(mean(tc(respwin_target,:,:),1)...
                        - mean(tc(basewin_0_target,:,:),1)));
                else
                    trResp = cat(2,trResp,squeeze(mean(tc(respwin,:,:),1)...
                        - mean(tc(basewin_0,:,:),1)));
                    trTC = cat(3,trTC,tc);
                end
                nfr = length(respwin)-1;
                nMovWin = 15;
                trResp_movWin_temp = nan(size(tc,2),size(tc,3),nMovWin);
                if ialign == 4
                    tempwinstart = 28;
                else
                    tempwinstart = 29;
                end
                for iwin = 1:nMovWin
                    tempwin = tempwinstart+(iwin-1):(tempwinstart+(iwin-1)+nfr);
                    if ialign == 4
                        trResp_movWin_temp(:,:,iwin) = squeeze(mean(tc(tempwin,:,:),1)) - ...
                            squeeze(mean(tc(basewin_0_target,:,:),1));
                    else
                        trResp_movWin_temp(:,:,iwin) = squeeze(mean(tc(tempwin,:,:),1)) - ...
                            squeeze(mean(tc(basewin_0,:,:),1));
                    end
                end
                trResp_movWin = cat(2,trResp_movWin,trResp_movWin_temp);
                if ialign == alignCR
                    trOut = cat(2,trOut,repmat({'cr'},[1,length(ddd.outcome)]));
                else
                    trOut = cat(2,trOut,ddd.outcome);
                end
                if isempty(ddd.ori) & isempty(ddd.amp)
                    trStim = cat(2,trStim,zeros(1,length(ddd.outcome)));
                elseif iav == 1
                    trStim = cat(2,trStim,ddd.ori);
                elseif iav == 2
                    trStim = cat(2,trStim,ddd.amp);
                end
            end
            trOut(strcmp(trOut,'success')) = {'h'};
            trOut(strcmp(trOut,'ignore')) = {'m'};
            trOut(strcmp(trOut,'failure')) = {'fa'};
            
            decodeDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            decodeDataExpt(exptN).av(iav).outcome = trOut;
            decodeDataExpt(exptN).av(iav).stim = trStim;
            decodeDataExpt(exptN).av(iav).resp = trResp;
            decodeDataExpt(exptN).av(iav).movWinResp = trResp_movWin;
            
            if iav == 1
                binEdges = oriBins;
            elseif iav == 2
                binEdges = ampBins;
            end
            trStimID = discretize(trStim,binEdges);
            rng(0)
            matchTrInd = getMatchedOutcomeTrialIndex(trOut,trStimID,minTrN);
            trOutTC = cell(1,length(trOutTypeName));
            for iout = 1:length(trOutType)
                if iout == 1 || iout == 2
                    for istim = 2:3
                        ind = trStimID == istim;
                        if istim == 2
                            ind2 = ((iout-1)*3)+2;
                        elseif istim == 3
                            ind2 = ((iout-1)*3)+3;
                        end
                        ind3 = ind & strcmp(trOut,trOutType(iout));
                        if sum(ind3) == 0
                            trOutTC{ind2} = nan(size(trTC,1),size(trTC,2),1);
                        else
                            trOutTC{ind2} = trTC(:,:,ind3);
                        end
%                         fprintf(num2str(ind2))
                    end
                    ind2 = ((iout-1)*3)+1;
                    trOutTC{ind2} = trTC(:,:,strcmp(trOut,trOutType(iout)) & matchTrInd);
%                     fprintf(num2str(ind2))
                else
%                     fprintf(num2str(iout+4))
                    trOutTC{iout+4} = trTC(:,:,strcmp(trOut,trOutType(iout)));
                end
            end
            targetDataExpt(exptN).av(iav).trOutTC = trOutTC;
            
            
            de = d.av(iav).align(alignStart);
            if strcmp(ds,'FSAV_attentionV1')
                hits = strcmp(de.outcome,'success');
            elseif strcmp(ds,'FSAV_V1_100ms_naive')
                hits = true(1,length(de.outcome));
            end
            misses = strcmp(de.outcome,'ignore');
            tc_AV = cat(3,tc_AV,de.respTC(:,:,hits)); % make hits or misses
%             for icyc = 1:nCycles
            for icyc = 1:maxCycles
                tc = de.respTC(:,:,de.nCycles >= icyc & (hits | misses));
                cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                cycTC{iav,icyc} = tc(...
                    (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
            end
            longTC{iav} = de.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                de.nCycles >= nCycLongTC & (hits | misses));
            
            de = d.av(iav).align(alignTarget);
            if iav == 1
                binEdges = oriBins;
                targets = de.ori;
            elseif iav == 2
                binEdges = ampBins;
                targets = de.amp;
            end
            tc = de.respTC(1:(nBaselineFr*2),:,:);
            targetInd = discretize(targets,binEdges);
            targetDataExpt(exptN).av(iav).tc = cell(1,2);
            for itar = 1:nTargets
                ind = targetInd == itar+1;
                targetDataExpt(exptN).av(iav).tc{itar} = tc(:,:,ind);
                targetDataExpt(exptN).av(iav).targets{itar} = targets(ind);
            end
        end
        antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
        antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
        antiDataExpt(exptN).longTC = longTC;
        antiDataExpt(exptN).cycTC = cycTC;
        targetDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
        do = d.oriTuning;
        if size(do.oriResp,2) == 8
            oriTuningExpt(exptN).oriResp = do.oriResp(:,1:2:8);
            oriTuningExpt(exptN).oriRespErr = do.oriRespSem(:,1:2:8);
        elseif size(do.oriResp,2) ~= 4
            error('error in orientations used for passivie tuning')
        else
            oriTuningExpt(exptN).oriResp = do.oriResp;
            oriTuningExpt(exptN).oriRespErr = do.oriRespSem;
        end
        oriTuningExpt(exptN).fit = do.oriFit;
        oriTuningExpt(exptN).isTuned = do.oriFitReliability < tuningReliabilityThresh;
        [~,oriTuningExpt(exptN).fitPeak] = max(do.oriFit,[],1);
        [~,oriPref] = histc(oriTuningExpt(exptN).fitPeak,oriBinEdges);
        oriPref(oriPref == length(orientations)+1 | oriPref == length(oriBinEdges) == 1) = 1;
        oriTuningExpt(exptN).oriPref = oriPref;
        oriTuningExpt(exptN).tuningReliability = do.oriFitReliability;
    end
end

shortCycExptInd = cell2mat({antiDataExpt.exptCycLengthFr}) == 11;
isShortCycExpt = [];
respCellsExpt = struct;
for iexp = 1:nexp
    firstTCAV = cat(3,antiDataExpt(iexp).cycTC{visualTrials,1},...
        antiDataExpt(iexp).cycTC{auditoryTrials,1});
    
    lateTC_vis = [];
    lateTC_aud = [];
    lateCycles = 5:size(antiDataExpt(iexp).cycTC,2);
    for icyc = 1:length(lateCycles)
        lateTC_vis = cat(3,lateTC_vis,...
            antiDataExpt(iexp).cycTC{visualTrials,lateCycles(icyc)});
        lateTC_aud = cat(3,lateTC_aud,...
            antiDataExpt(iexp).cycTC{auditoryTrials,lateCycles(icyc)});
    end
    lateTCAV = cat(3,lateTC_vis,lateTC_aud);
    
    longTCAV = cat(3,antiDataExpt(iexp).longTC{visualTrials,1},...
        antiDataExpt(iexp).longTC{auditoryTrials,1});
    
    firstRespCells = ttest(...
        squeeze(mean(firstTCAV(respwin,:,:),1)),...
        squeeze(mean(firstTCAV(basewin_0,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha);
    lateRespCells = ttest(...
        squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
        squeeze(mean(longTCAV(basewin,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha);
    lateSuppCells = ttest(...
        squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
        squeeze(mean(longTCAV(basewin,:,:),1)),...
        'dim',2,'tail','left','alpha',cellGroupsAlpha);
    lateCycRespCells = ttest(...
        squeeze(mean(lateTCAV(respwin,:,:),1)),...
        squeeze(mean(lateTCAV(basewin_0,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha);
    
    lateCycRespCutoffPass = mean(mean(lateTCAV(respwin,:,:),3),1) > minRespThreshold_decode;
    
    eaTarRespCells = sum(cell2mat(cellfun(@(x) ...
        ttest(squeeze(mean(x(respwin_target,:,:),1)),...
        squeeze(mean(x(basewin_0_target,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha),...
        targetDataExpt(iexp).av(visualTrials).tc,'unif',0)),2) > 0;
    
    
    allTargetTC = [];
    for itar = 1:length(targetDataExpt(iexp).av(visualTrials).tc)
        allTargetTC = cat(3,allTargetTC,...
            targetDataExpt(iexp).av(visualTrials).tc{itar});
    end
    allTarRespCells = ttest(squeeze(mean(allTargetTC(respwin_target,:,:),1)),...
        squeeze(mean(allTargetTC(basewin_0_target,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha);
    
    eaTarRespCutoffPass = sum(cell2mat(cellfun(@(x) mean(mean(x(respwin_target,:,:),3),1) - mean(mean(x(basewin_0_target,:,:),3),1),...
        targetDataExpt(iexp).av(visualTrials).tc,'unif',0)') > minRespThreshold_decode,1) > 0 ...
        | (mean(mean(allTargetTC(respwin_target,:,:),3),1) - mean(mean(allTargetTC(basewin_0_target,:,:),3),1)) > minRespThreshold_decode; 
        
    respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
    respCellsExpt(iexp).firstRespCells = firstRespCells;
    respCellsExpt(iexp).lateRespCells = lateRespCells;
    respCellsExpt(iexp).lateSuppCells = lateSuppCells;
    respCellsExpt(iexp).lateCycRespCells = lateCycRespCells;
    respCellsExpt(iexp).targetRespCells = eaTarRespCells | allTarRespCells;
%     respCellsExpt(iexp).decodeAnalysisCells = ...
%         (firstRespCells & mean(mean(firstTCAV(respwin,:,:),3),1)' > minRespThreshold)...
%         | ((eaTarRespCells | allTarRespCells) & eaTarRespCutoffPass');
    respCellsExpt(iexp).decodeAnalysisCells = ...
        (lateCycRespCells & lateCycRespCutoffPass')...
        | ((eaTarRespCells | allTarRespCells) & eaTarRespCutoffPass');
    if shortCycExptInd(iexp)
        isShortCycExpt = cat(1,isShortCycExpt,true(length(firstRespCells),1));
    else
        isShortCycExpt = cat(1,isShortCycExpt,false(length(firstRespCells),1));
    end
        
end

decodeAnalysis = struct;
decodeAnalysis.av(visualTrials).name = 'Visual';
decodeAnalysis.av(auditoryTrials).name = 'Auditory';
for iexp = 1:nexp
    rng(0) % cells randomized and trials randomized
    cellInd = respCellsExpt(iexp).decodeAnalysisCells & ...
        oriTuningExpt(iexp).tuningReliability' <= tuningReliabilityThresh_decode;
    decodeAnalysis(iexp).nCells = sum(cellInd);
    if sum(cellInd) > maxCellN
        ind = find(cellInd);
        cellSampleID = randsample(ind,maxCellN);
        cellInd = false(length(cellInd),1);
        cellInd(cellSampleID) = true;
    end
    decodeAnalysis(iexp).nCellsSelected = sum(cellInd);
    decodeAnalysis(iexp).cellInd = cellInd;
    respOther = cell(1,2);
    trOutOther = cell(1,2);
    detectGLMOther = cell(1,2);
    targetGLMOther = cell(1,2);
    for iav = 1:2
        trOut = decodeDataExpt(iexp).av(iav).outcome;
%         [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut);        
        respAllCells = zscore(decodeDataExpt(iexp).av(iav).resp');
        
        if iav == 1
            trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
        elseif iav == 2
            trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
        end
        nStimPerBin = histcounts(trStimID);
        minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
        respStimSort = cell(1,nStimBins);
        trOutStimSort = cell(1,nStimBins);
        for istim = 1:nStimBins
            ind = find(trStimID == istim);
            if length(ind) >= minTrN_mdl
                if istim == 1
                    matchTrialsInd = [];
                    if sum(nStimPerBin >= minTrN_mdl) == 2
                        n = minBinN;
                    elseif minBinN == nStimPerBin(istim)
                        error('not enough FA/CR trials')
                    else
                        n = (nStimBins-1).*minBinN;
                        if n > length(ind)
                            error('not enough FA/CR trials')
                        end
                    end
                    indSample = randsample(ind,n);
                    matchTrialsInd = cat(2,matchTrialsInd,indSample);
                else
                    indSample = randsample(ind,minBinN);
                    matchTrialsInd = cat(2,matchTrialsInd,indSample);
                end
                respStimSort{istim} = respAllCells(indSample,cellInd);
                trOutStimSort{istim} = trOut(indSample);
            end
        end
        nMatchedTrials = cumsum(cellfun(@length,trOutStimSort));
        for istim = 1:nStimBins
            if istim == 1
                stimSortInd = cell(1,nStimBins);
                stimSortInd{istim} = 1:nMatchedTrials;
            else
                stimSortInd{istim} = ...
                    (nMatchedTrials(istim-1)+1):nMatchedTrials(istim);
            end
        end
        
%         matchTrialsInd = 1:size(respAllCells,2);
%         resp = zscore(respAllCells(cellInd,matchTrialsInd)');
        resp = respAllCells(matchTrialsInd,cellInd);
%         [c,score,latent] = pca(decodeDataExpt(iexp).av(iav).resp(:,matchTrialsInd)');
%         resp_pcs = zscore(score(:,1:10));
%         [detectTrIndAll, targetTrIndAll] = getStimAndBehaviorYs(trOut);
        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));
        
        detectCorr = corr(detectTrInd,resp);
        targetCorr = corr(targetTrInd,resp);
        
        C = eye(size(resp,2));
        p=1;
        [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
        [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
        
        fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
        detectWeight = detectGLM.beta(2:end);
        targetWeight = targetGLM.beta(2:end);
        
        dv_detect = mean(detectTrInd);
        dv_target = mean(targetTrInd);
        
        pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp,detectTrInd,dv_detect);
        pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);
        
        pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
        pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
        
        pctCorrDetect_xStim_train = nan(1,nStimBins);
        pctCorrDetect_xStim_ho = nan(1,nStimBins);
        pctCorrTarget_xStim_train = nan(1,nStimBins);
        pctCorrTarget_xStim_ho = nan(1,nStimBins);
        for istim = 1:nStimBins
            if nStimPerBin(istim) >= minTrN_mdl
                [detectStimInd, targetStimInd] = getStimAndBehaviorYs(...
                    trOutStimSort{istim});
                pctCorrDetect_xStim_train(istim) = getPctCorr_trainData(...
                    detectGLM,respStimSort{istim},detectStimInd,dv_detect);
                pctCorrDetect_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                    resp,detectTrInd,stimSortInd{istim},dv_detect);
                pctCorrTarget_xStim_train(istim) = getPctCorr_trainData(...
                    targetGLM,respStimSort{istim},targetStimInd,dv_target);
                pctCorrTarget_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                    resp,targetTrInd,stimSortInd{istim},dv_target);
            end
        end
       
        fprintf('Expt %s, starting resp win analysis\n',num2str(iexp))
        nwins = size(decodeDataExpt(iexp).av(iav).movWinResp,3);
        pctCorrectDetect_train_respwin = nan(1,nwins);
        pctCorrectDetect_ho_respwin = nan(1,nwins);
        pctCorrectTarget_train_respwin = nan(1,nwins);
        pctCorrectTarget_ho_respwin = nan(1,nwins);
        for iwin = 1:nwins
            rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
            r = rAllCells(matchTrialsInd,cellInd);
            [~,~,detectGLM_temp] = glmfit(r*C,detectTrInd,'binomial');
            [~,~,targetGLM_temp] = glmfit(r*C,targetTrInd,'binomial');
            
            pctCorrectDetect_train_respwin(iwin) = getPctCorr_trainData(...
                detectGLM_temp,r,detectTrInd,dv_detect);
            pctCorrectDetect_ho_respwin(iwin) = getPctCorr_hoData(r,detectTrInd,dv_detect);
            pctCorrectTarget_train_respwin(iwin) = getPctCorr_trainData(...
                targetGLM_temp,r,targetTrInd,dv_target);
            pctCorrectTarget_ho_respwin(iwin) = getPctCorr_hoData(r,targetTrInd,dv_target);
        end
        
        respOther{iav} = resp;
        trOutOther{iav} = trOut(matchTrialsInd);
        detectGLMOther{iav} = detectGLM;
        targetGLMOther{iav} = targetGLM;        
        
        decodeAnalysis(iexp).av(iav).dvDetect = dv_detect;
        decodeAnalysis(iexp).av(iav).dvTarget = dv_target;
        decodeAnalysis(iexp).av(iav).correlationDetect = detectCorr;
        decodeAnalysis(iexp).av(iav).correlationTarget = targetCorr;
        decodeAnalysis(iexp).av(iav).weightDetect = detectWeight;
        decodeAnalysis(iexp).av(iav).weightTarget = targetWeight;
        decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_train = pctCorrectDetect_train;
        decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_holdout = pctCorrectDetect_ho;
        decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_train = pctCorrectTarget_train;
        decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout = pctCorrectTarget_ho;
        decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_train = pctCorrDetect_xStim_train;
        decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_holdout = pctCorrDetect_xStim_ho;
        decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_train = pctCorrTarget_xStim_train;
        decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_holdout = pctCorrTarget_xStim_ho;
        decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_train = pctCorrectDetect_train_respwin;
        decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_holdout = pctCorrectDetect_ho_respwin;
        decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_train = pctCorrectTarget_train_respwin;
        decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_holdout = pctCorrectTarget_ho_respwin;
    end
    fprintf('Expt %s: Testing opposite model...\n',num2str(iexp))
    for iav = 1:2
        if iav == 1
            otherAV = 2;
        elseif iav == 2
            otherAV = 1;
        end
        resp = respOther{otherAV};
        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});
        
        dv_detect = mean(detectTrInd);
        dv_target = mean(targetTrInd);
        
        pctCorrectDetect = getPctCorr_trainData(detectGLMOther{iav},resp,detectTrInd,dv_detect);
        pctCorrectTarget = getPctCorr_trainData(targetGLMOther{iav},resp,targetTrInd,dv_target);
        
        decodeAnalysis(iexp).av(iav).pctCorrectDetect_otherAV = pctCorrectDetect;
        decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = pctCorrectTarget;
    end
end


nBoot = 1000;
dcBootstrap = struct;
for iexp = 1:nexp
    fprintf('Expt %s\n',num2str(iexp))
    rng(0)
    for iav = 1:2
        trOut = decodeDataExpt(iexp).av(iav).outcome;
        respAllCells = zscore(decodeDataExpt(iexp).av(iav).resp');

        if iav == 1
            trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
        elseif iav == 2
            trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
        end
        nStimPerBin = histcounts(trStimID);
        minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
        for istim = 1:nStimBins
            ind = find(trStimID == istim);
            if length(ind) >= minTrN_mdl
                if istim == 1
                    matchTrialsInd = [];
                    if sum(nStimPerBin >= minTrN_mdl) == 2
                        n = minBinN;
                    elseif minBinN == nStimPerBin(istim)
                        error('not enough FA/CR trials')
                    else
                        n = (nStimBins-1).*minBinN;
                        if n > length(ind)
                            error('not enough FA/CR trials')
                        end
                    end
                    indSample = randsample(ind,n);
                    matchTrialsInd = cat(2,matchTrialsInd,indSample);
                else
                    indSample = randsample(ind,minBinN);
                    matchTrialsInd = cat(2,matchTrialsInd,indSample);
                end
            end
        end
        for iboot = 1:nBoot
            fprintf('.')

            nTotalCells = length(respCellsExpt(iexp).decodeAnalysisCells);
            cellInd = randsample(nTotalCells,maxCellN);

            resp = respAllCells(matchTrialsInd,cellInd);
            [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

            detectCorr = corr(detectTrInd,resp);
            targetCorr = corr(targetTrInd,resp);

            C = eye(size(resp,2));
            p=1;
            [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
            [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');

            detectWeight = detectGLM.beta(2:end);
            targetWeight = targetGLM.beta(2:end);
            
            dcBootstrap(iexp).av(iav).cellInd{iboot} = cellInd;
            dcBootstrap(iexp).av(iav).targetCorrelation{iboot} = targetCorr;
            dcBootstrap(iexp).av(iav).detectCorrelation{iboot} = detectCorr;
            dcBootstrap(iexp).av(iav).targetWeight{iboot} = targetWeight;
            dcBootstrap(iexp).av(iav).detectWeight{iboot} = detectWeight;
        end
    end
    dcBootstrap(iexp).nTotalCells = nTotalCells;
end

targetWeightBoot = cell(1,2);
detectWeightBoot = cell(1,2);
for iexp = 1:nexp
    n = dcBootstrap(iexp).nTotalCells;
    for iav = 1:2
        d = dcBootstrap(iexp).av(iav);
        bootMatrix_w_target = nan(nBoot,n);
        bootMatrix_w_detect = nan(nBoot,n);
        for iboot = 1:nBoot
            ind = d.cellInd{iboot};
            bootMatrix_w_target(iboot,ind) = d.targetWeight{iboot};
            bootMatrix_w_detect(iboot,ind) = d.detectWeight{iboot};
        end
        targetWeightBoot{iav} = cat(2,targetWeightBoot{iav},nanmean(bootMatrix_w_target));
        detectWeightBoot{iav} = cat(2,detectWeightBoot{iav},nanmean(bootMatrix_w_detect));
    end
end


antiAnalysis = struct;
antiAnalysis.longTC = cell(1,3);
antiAnalysis.longTCErr = cell(1,3);
antiAnalysis.lateCycTC = cell(1,3);
antiAnalysis.firstCycTC = cell(1,3);
antiAnalysis.lateCycSI = [];
antiAnalysis.lateCycFracSI = [];
targetAnalysis = struct;
targetAnalysis.tc = cell(3,2);
targetAnalysis.targets = cell(1,2);
targetAnalysis.targetsAuROC = [];
targetAnalysis.targetsAuROCTest = [];
targetAnalysis.trOutTC = cell(length(trOutTypeName),2);
firstStimAuROC = [];
lateStimAuROC = [];
firstStimAuROCTest = [];
lateStimAuROCTest = [];
taskTuningTest = [];
taskTuningPref = [];
% taskTuningTest_temp = [];
% taskTuningPref_temp = [];
for iexp = 1:nexp
    longTC_vis = antiDataExpt(iexp).longTC{visualTrials};
    longTC_aud = antiDataExpt(iexp).longTC{auditoryTrials};
    
    nCells = size(longTC_vis,2);
    
    firstCycResp_vis = antiDataExpt(iexp).cycTC{visualTrials,1};
    firstCycResp_aud = antiDataExpt(iexp).cycTC{auditoryTrials,1};
    
    cycTC_vis = antiDataExpt(iexp).cycTC(visualTrials,:);
    cycTC_aud = antiDataExpt(iexp).cycTC(auditoryTrials,:);
    
    lateCycTC_vis = [];
    lateCycTC_aud = [];
    lateCycles = 5:size(antiDataExpt(iexp).cycTC,2);
    for icyc = 1:length(lateCycles)
        lateCycTC_vis = cat(3,lateCycTC_vis,...
            cycTC_vis{lateCycles(icyc)} - mean(cycTC_vis{lateCycles(icyc)}(basewin_0,:,:),1));
        lateCycTC_aud = cat(3,lateCycTC_aud,...
            cycTC_aud{lateCycles(icyc)} - mean(cycTC_aud{lateCycles(icyc)}(basewin_0,:,:),1));
    end
    
    lateCycSI = getSelectivityIndex(squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
        squeeze(mean(lateCycTC_aud(respwin,:,:),1))');
    v=mean(mean(lateCycTC_vis(respwin,:,:),3),1);
    v(v<0) = 0;
    a=mean(mean(lateCycTC_aud(respwin,:,:),3),1);
    a(a<0) = 0;
    lateCycFracSI = (v-a)./v;
    lateCycFracSI_aud = (a-v)./a;
    lateCycFracSI(a > v) = -lateCycFracSI_aud(a > v);
    
    firstCycTC_vis = antiDataExpt(iexp).cycTC{visualTrials,1};
    firstCycTC_aud = antiDataExpt(iexp).cycTC{auditoryTrials,1};
    
    firstVisResp = squeeze(mean(firstCycTC_vis(respwin,:,:),1));
    lateVisResp = squeeze(mean(lateCycTC_vis(respwin,:,:),1));
    
    tarVisResp = cellfun(@(x) ...
        squeeze(mean(x(respwin_target,:,:),1) - mean(x(basewin_0_target,:,:),1)),...
        cat(2,targetDataExpt(iexp).av(visualTrials).tc,...
        {cat(3,targetDataExpt(iexp).av(visualTrials).tc{1},...
        targetDataExpt(iexp).av(visualTrials).tc{2})}),'unif',0);
    auroc_first = nan(nCells,3);
    auroc_first_test = nan(nCells,3);
    auroc_late = nan(nCells,3);
    auroc_late_test = nan(nCells,3);
    auroc_target = nan(nCells,1);
    auroc_target_test = nan(nCells,1);
    for icell = 1:nCells
        auroc_target(icell) = roc_gh(tarVisResp{1}(icell,:),tarVisResp{2}(icell,:));
        auroc_target_test(icell) = ranksum(...
            tarVisResp{1}(icell,:),tarVisResp{2}(icell,:)) < 0.05;
        
        auroc_first(icell,:) = cellfun(@(x) ...
            roc_gh(firstVisResp(icell,:),x(icell,:)),tarVisResp);
        auroc_first_test(icell,:) = cellfun(@(x) ...
            ranksum(firstVisResp(icell,:),x(icell,:)),tarVisResp);
        auroc_late(icell,:) = cellfun(@(x) ...
            roc_gh(lateVisResp(icell,:),x(icell,:)),tarVisResp);
        auroc_late_test(icell,:) = cellfun(@(x) ...
            ranksum(lateVisResp(icell,:),x(icell,:)),tarVisResp);
    end
    
    stimResp4anova = cat(2,firstVisResp,tarVisResp{1},tarVisResp{2});
    stimID4anova = cat(2,ones(1,size(firstVisResp,2)),...
        ones(1,size(tarVisResp{1},2))*2,ones(1,size(tarVisResp{2},2))*3);
    taskStimAnova = nan(1,nCells);
    for icell = 1:nCells
        taskStimAnova(icell) = anova1(stimResp4anova(icell,:),stimID4anova,'off');
    end
    [~,taskTuningID] = max(cat(2,mean(firstVisResp,2),...
        mean(tarVisResp{1},2),mean(tarVisResp{2},2)),[],2);
    
    
    antiAnalysis.longTC{visualTrials} = cat(2,antiAnalysis.longTC{visualTrials},...
        mean(longTC_vis,3));
    antiAnalysis.longTC{auditoryTrials} = cat(2,antiAnalysis.longTC{auditoryTrials},...
        mean(longTC_aud,3));
    antiAnalysis.longTC{allTrialsInd} = cat(2,antiAnalysis.longTC{allTrialsInd},...
        mean(cat(3,longTC_vis,longTC_aud),3));
    antiAnalysis.longTCErr{allTrialsInd} = cat(2,antiAnalysis.longTCErr{allTrialsInd},...
        ste(cat(3,longTC_vis,longTC_aud),3));
    antiAnalysis.lateCycTC{visualTrials} = cat(2,antiAnalysis.lateCycTC{visualTrials},...
        mean(lateCycTC_vis,3));
    antiAnalysis.lateCycTC{auditoryTrials} = cat(2,antiAnalysis.lateCycTC{auditoryTrials},...
        mean(lateCycTC_aud,3));
    antiAnalysis.lateCycTC{allTrialsInd} = cat(2,antiAnalysis.lateCycTC{allTrialsInd},...
        mean(cat(3,lateCycTC_vis,lateCycTC_aud),3));
    antiAnalysis.lateCycSI = cat(2,antiAnalysis.lateCycSI,lateCycSI);
    antiAnalysis.lateCycFracSI = cat(2,antiAnalysis.lateCycFracSI,lateCycFracSI);
    antiAnalysis.firstCycTC{visualTrials} = cat(2,antiAnalysis.firstCycTC{visualTrials},...
        mean(firstCycTC_vis,3));
    antiAnalysis.firstCycTC{auditoryTrials} = cat(2,antiAnalysis.firstCycTC{auditoryTrials},...
        mean(firstCycTC_aud,3));
    antiAnalysis.firstCycTC{allTrialsInd} = cat(2,antiAnalysis.firstCycTC{allTrialsInd},...
        mean(cat(3,firstCycTC_vis,firstCycTC_aud),3));
    
    targetAnalysis.tc{1,visualTrials} = cat(2,targetAnalysis.tc{1,visualTrials},...
        mean(targetDataExpt(iexp).av(visualTrials).tc{1},3));
    targetAnalysis.tc{2,visualTrials} = cat(2,targetAnalysis.tc{2,visualTrials},...
        mean(targetDataExpt(iexp).av(visualTrials).tc{2},3));
    targetAnalysis.tc{3,visualTrials} = cat(2,targetAnalysis.tc{3,visualTrials},...
        mean(cat(3,targetDataExpt(iexp).av(visualTrials).tc{1},...
        targetDataExpt(iexp).av(visualTrials).tc{2}),3));
    targetAnalysis.tc{1,auditoryTrials} = cat(2,targetAnalysis.tc{1,auditoryTrials},...
        mean(targetDataExpt(iexp).av(auditoryTrials).tc{1},3));
    targetAnalysis.tc{2,auditoryTrials} = cat(2,targetAnalysis.tc{2,auditoryTrials},...
        mean(targetDataExpt(iexp).av(auditoryTrials).tc{2},3));
    targetAnalysis.tc{3,auditoryTrials} = cat(2,targetAnalysis.tc{3,auditoryTrials},...
        mean(cat(3,targetDataExpt(iexp).av(auditoryTrials).tc{1},...
        targetDataExpt(iexp).av(auditoryTrials).tc{2}),3));
    targetAnalysis.trOutTC(:,visualTrials) = cellfun(@(x,y) cat(2,x,mean(y(1:(nBaselineFr*2),:,:),3)),...
        targetAnalysis.trOutTC(:,visualTrials),...
        targetDataExpt(iexp).av(visualTrials).trOutTC','unif',0);
    targetAnalysis.trOutTC(:,auditoryTrials) = cellfun(@(x,y) cat(2,x,mean(y(1:(nBaselineFr*2),:,:),3)),...
        targetAnalysis.trOutTC(:,auditoryTrials),...
        targetDataExpt(iexp).av(auditoryTrials).trOutTC','unif',0);
    
    targetAnalysis.targets{visualTrials} = cat(2,targetAnalysis.targets{visualTrials},...
        cat(1,ones(1,nCells).*mean(targetDataExpt(iexp).av(visualTrials).targets{1}),...
        ones(1,nCells).*mean(targetDataExpt(iexp).av(visualTrials).targets{2})));
    targetAnalysis.targetsAuROC = cat(1,targetAnalysis.targetsAuROC,auroc_target);
    targetAnalysis.targetsAuROCTest = cat(1,targetAnalysis.targetsAuROCTest,auroc_target_test);
    taskTuningTest = cat(1,taskTuningTest,taskStimAnova' < 0.05);
    taskTuningPref = cat(1,taskTuningPref,taskTuningID);
    
    firstStimAuROC = cat(1,firstStimAuROC,auroc_first);
    lateStimAuROC = cat(1,lateStimAuROC,auroc_late);
    firstStimAuROCTest = cat(1,firstStimAuROCTest,auroc_first_test);
    lateStimAuROCTest = cat(1,lateStimAuROCTest,auroc_late_test);
    
%     ind = taskStimAnova' > 0.05 & ttest(firstVisResp,tarVisResp{3},'dim',2);
%     taskTuningID_temp = taskTuningID;
%     taskTuningID_temp(ind) = 4;
%     taskTuningPref_temp =cat(1,taskTuningPref_temp,taskTuningID_temp);
end
% antiAnalysis.adapt = cellfun(@(x,y) ...
%     mean(x(respwin,:),1)./mean(y(respwin,:),1),...
%     antiAnalysis.lateCycTC,antiAnalysis.firstCycTC,'unif',0);

cellInfo = struct;
cellInfo.firstRespCells = logical(cell2mat({respCellsExpt.firstRespCells}'));
cellInfo.lateRespCells = logical(cell2mat({respCellsExpt.lateRespCells}'));
cellInfo.lateSuppCells = logical(cell2mat({respCellsExpt.lateSuppCells}'));
cellInfo.lateCycRespCells = logical(cell2mat({respCellsExpt.lateCycRespCells}'));
cellInfo.minRespCells = (mean(antiAnalysis.firstCycTC{visualTrials}(respwin,:),1) > ...
    minRespThreshold & mean(antiAnalysis.firstCycTC{auditoryTrials}(respwin,:),1) > ...
    minRespThreshold)';
cellInfo.targetRespCells = logical(cell2mat({respCellsExpt.targetRespCells}'));
cellInfo.isShortCycExpt = isShortCycExpt;
cellInfo.isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
cellInfo.oriResp = cell2mat({oriTuningExpt.oriResp}');
cellInfo.oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
cellInfo.oriFit = cell2mat({oriTuningExpt.fit})';
cellInfo.oriPref = cell2mat({oriTuningExpt.oriPref})';
cellInfo.taskTuningPref = taskTuningPref;
% cellInfo.taskTuningPref = taskTuningPref_temp;
cellInfo.taskTuningTest = taskTuningTest;

auroc_first = nan(length(cellInfo.firstRespCells),1);
auroc_late = nan(length(cellInfo.firstRespCells),1);
auroc_first_test = nan(length(cellInfo.firstRespCells),1);
auroc_late_test = nan(length(cellInfo.firstRespCells),1);
for icell = 1:length(cellInfo.firstRespCells)
    if ~cellInfo.targetRespCells(icell)
        auroc_first(icell) = firstStimAuROC(icell,3);
        auroc_late(icell) = lateStimAuROC(icell,3);
        auroc_first_test(icell) = firstStimAuROCTest(icell,3) < 0.05;
        auroc_late_test(icell) = lateStimAuROCTest(icell,3) < 0.05;
    elseif ~targetAnalysis.targetsAuROCTest(icell)
        auroc_first(icell) = firstStimAuROC(icell,3);
        auroc_late(icell) = lateStimAuROC(icell,3);
        auroc_first_test(icell) = firstStimAuROCTest(icell,3) < 0.05;
        auroc_late_test(icell) = lateStimAuROCTest(icell,3) < 0.05;
    elseif targetAnalysis.targetsAuROCTest(icell) && targetAnalysis.targetsAuROC(icell) < 0.5
        auroc_first(icell) = firstStimAuROC(icell,1);
        auroc_late(icell) = lateStimAuROC(icell,1);
        auroc_first_test(icell) = firstStimAuROCTest(icell,1) < 0.05;
        auroc_late_test(icell) = lateStimAuROCTest(icell,1) < 0.05;
    elseif targetAnalysis.targetsAuROCTest(icell) && targetAnalysis.targetsAuROC(icell) >= 0.5
        auroc_first(icell) = firstStimAuROC(icell,2);
        auroc_late(icell) = lateStimAuROC(icell,2);
        auroc_first_test(icell) = firstStimAuROCTest(icell,2) < 0.05;
        auroc_late_test(icell) = lateStimAuROCTest(icell,2) < 0.05;
    end
end
cellInfo.firstStimAuROC = auroc_first;
cellInfo.lateStimAuROC = auroc_late;
cellInfo.firstStimAuROCTest = auroc_first_test;
cellInfo.lateStimAuROCTest = auroc_late_test;

cellInfo.dcModelCells = logical(cell2mat({decodeAnalysis.cellInd}'));
cellInfo.targetWeight = targetWeightBoot;
cellInfo.detectWeight = detectWeightBoot;

save([fnout 'imgAnalysisData'],...
    'decodeAnalysis','antiAnalysis','targetAnalysis','cellInfo')
%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.005 0.015];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];
siOriLim = [-3 3];
oriBarLim_win = [0 0.08];
oriBarLim_resp = [0 0.04];
oriLim_taskResp = [-0.005 0.035];
oriNLim = [0 120];
oriTCLim = [-0.005 0.08];
targetTCLim = [-0.015 0.08];
outTCLim = [-0.005 0.04];
firstTCLim = [-0.005 0.04];
% adaptLim = [0 1];

tcStartFrame = 26;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = size(antiAnalysis.longTC{1,1},1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = size(antiAnalysis.lateCycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT_target = (...
    [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
    .*(1000/frameRateHz);
baseWinTT = (...
    [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);

weightLim = [-3 3];
binnedWeightLim = [-0.4 0.4];
weightLimSum = [-0.6 0.6];
siLimSum = [-1 1];
%% plot anticipation analysis (Figure 2)

setFigParams4Print('landscape')

% heatmaps
lateWinRespAll = mean(antiAnalysis.longTC{allTrialsInd}(lateWinFr,:),1);
[~,lateWinSortInd] = sort(lateWinRespAll);
hm = flipud(antiAnalysis.longTC{allTrialsInd}(:,lateWinSortInd)');

figure
colormap(brewermap([],'*RdBu'));
% subplot 121
imagesc(hm(:,tcStartFrame:end))
hold on
if strcmp(ds,'FSAV_attentionV1')
    exCellInd = [exampleCell_1,exampleCell_2];
    exCellMat = zeros(1,length(cellInfo.firstRespCells));
    exCellMat(exCellInd) = 1;
    exCellSortInd = find(flip(exCellMat(lateWinSortInd)));
    hline(exCellSortInd,'k-')
end
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabelFr_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title('All Trials, All Cells')
print([fnout 'heatmapAllTrialsAllCells'],'-dpdf','-fillpage')

figure
colormap gray
% subplot 122
ind = cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells;
indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
imagesc(flipud(indBW(lateWinSortInd)))
figYAxis([],'Cell #',[])
figAxForm
colorbar
title('Anti. Resp. Cells')
print([fnout 'cellIDforHeatmap'],'-dpdf','-fillpage')

% time-courses and quantification scatters
setFigParams4Print('portrait')
figure
subplot 311
ind = cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells;
for iav = 1:2
    y = mean(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & ...
        cellInfo.isShortCycExpt),2);
    yerr = ste(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & cellInfo.isShortCycExpt),2);
    hold on
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',respTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
figAxForm([],0)
title(sprintf('First Stim Responsive Cells (%s/%s)',...
    num2str(sum((...
    cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)...
    & cellInfo.isShortCycExpt)),...
    num2str(sum(cellInfo.isShortCycExpt))))

subplot 323
x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(scatLim_win,scatLim_win,'k--')
plot(scatLim_win,[0 0],'k:')
plot([0 0],scatLim_win,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim_win)
figYAxis([],'Auditory (dF/F)',scatLim_win)
figAxForm
title(sprintf('Late Window, All Resp. Cells (%s/%s), p = %s',...
    num2str(sum(cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))


subplot 325
for iav = 1:2
    y = mean(antiAnalysis.lateCycTC{iav}...
        ((tcStartFrame:end),(cellInfo.lateCycRespCells)),2);
    yerr = ste(antiAnalysis.lateCycTC{iav}...
        ((tcStartFrame:end),(cellInfo.lateCycRespCells)),2);
    hold on
    shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',cycTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title(sprintf('Late Cycle Resp. Cells (%s/%s)',...
    num2str(sum(cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells))))

subplot 326
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(scatLim_cyc,scatLim_cyc,'k--')
plot(scatLim_cyc,[0 0],'k:')
plot([0 0],scatLim_cyc,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim_cyc)
figYAxis([],'Auditory (dF/F)',scatLim_cyc)
figAxForm
title(sprintf('Late Cycles, Late Cyc. Resp. Cells (%s/%s), p = %s',...
    num2str(sum(cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))

print([fnout 'tcLongTrialsAndLateCycWithQuant_LateCycRespCells'],'-dpdf','-fillpage')

% example cells time-courses
if strcmp(ds,'FSAV_attentionV1')
    setFigParams4Print('landscape')

    figure
    subplot 221
    y = antiAnalysis.longTC{allTrialsInd}(tcStartFrame:end,exampleCell_1);
    yerr = antiAnalysis.longTCErr{allTrialsInd}(tcStartFrame:end,exampleCell_1);
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    hold on
    figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
    figYAxis([],'dF/F',exCellTCLim)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Example Cell #%s',num2str(exampleCell_1)))

    subplot 222
    y = cellInfo.oriResp(exampleCell_1,:);
    yerr = cellInfo.oriRespErr(exampleCell_1,:);
    errorbar(orientations,y,yerr,'.')
    hold on
    x = 0:180;
    y = cellInfo.oriFit(exampleCell_1,:);
    plot(x,y,'-')
    figXAxis([],'Orienation (deg)',[-10 190])
    figYAxis([],'dF/F',oriRespLim)
    figAxForm
    title(sprintf('Passive Ori. Tuning, Cell #%s',num2str(exampleCell_1)));

    subplot 223
    y = antiAnalysis.longTC{allTrialsInd}(tcStartFrame:end,exampleCell_2);
    yerr = antiAnalysis.longTCErr{allTrialsInd}(tcStartFrame:end,exampleCell_2);
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    hold on
    figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
    figYAxis([],'dF/F',exCellTCLim)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Example Cell #%s',num2str(exampleCell_2)))

    subplot 224
    y = cellInfo.oriResp(exampleCell_2,:);
    yerr = cellInfo.oriRespErr(exampleCell_2,:);
    errorbar(orientations,y,yerr,'.')
    hold on
    x = 0:180;
    y = cellInfo.oriFit(exampleCell_2,:);
    plot(x,y,'-')
    figXAxis([],'Orienation (deg)',[-10 190])
    figYAxis([],'dF/F',oriRespLim)
    figAxForm
    title(sprintf('Passive Ori. Tuning, Cell #%s',num2str(exampleCell_2)));

    print([fnout 'exampleCellsTCWithTuning'],'-dpdf','-fillpage')
end

% selectivity
figure
y = antiAnalysis.lateCycSI(cellInfo.lateCycRespCells);
h = cdfplot(y);
hold on;
vline(mean(y),'k-')
figXAxis([],'Selectivity Index',siLim)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
title(sprintf('Late Cyc. Resp. Cells, mean = %s, ste = %s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant'))))

print([fnout 'selectivityCDF'],'-dpdf')

%% create a structure with stats from analyses
imgStats = struct;
imgStats.lateWinTimeMs = round(lateWinTT,4,'significant');

imgStats.nCells.totalCells = length(cellInfo.firstRespCells);
imgStats.nCells.firstResp = sum(cellInfo.firstRespCells);
imgStats.nCells.lateCycResp = sum(cellInfo.lateCycRespCells);
imgStats.nCells.lateCycOnlyResp = sum(...
    ~cellInfo.firstRespCells & cellInfo.lateCycRespCells);
imgStats.nCells.lateTCResp = sum(cellInfo.lateRespCells);
imgStats.nCells.lateSupp = sum(cellInfo.lateSuppCells);
imgStats.nCells

x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
imgStats.av(visualTrials).lateWinResp = mean(x);
imgStats.av(visualTrials).lateWinRespErr = ste(x,2);
imgStats.av(auditoryTrials).lateWinResp = mean(y);
imgStats.av(auditoryTrials).lateWinRespErr = ste(y,2);
[~,imgStats.lateWinTest] = ttest(x,y);

fprintf('Visual - Mean/Err Late Win: %s/%s\n',...
    num2str(round(imgStats.av(visualTrials).lateWinResp,2,'significant')),...
    num2str(round(imgStats.av(visualTrials).lateWinRespErr,2,'significant')))
fprintf('Auditory - Mean/Err Late Win: %s/%s\n',...
    num2str(round(imgStats.av(auditoryTrials).lateWinResp,2,'significant')),...
    num2str(round(imgStats.av(auditoryTrials).lateWinRespErr,2,'significant')))

x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
imgStats.av(visualTrials).lateCycResp = mean(x);
imgStats.av(visualTrials).lateCycRespErr = ste(x,2);
imgStats.av(auditoryTrials).lateCycResp = mean(y);
imgStats.av(auditoryTrials).lateCycRespErr = ste(y,2);
[~,imgStats.lateCycTest] = ttest(x,y);

fprintf('Visual - Mean/Err Late Cyc: %s/%s\n',...
    num2str(round(imgStats.av(visualTrials).lateCycResp,2,'significant')),...
    num2str(round(imgStats.av(visualTrials).lateCycRespErr,2,'significant')))
fprintf('Auditory - Mean/Err Late Cyc: %s/%s\n',...
    num2str(round(imgStats.av(auditoryTrials).lateCycResp,2,'significant')),...
    num2str(round(imgStats.av(auditoryTrials).lateCycRespErr,2,'significant')))

imgStats.allRespCellsSI = antiAnalysis.lateCycSI(...
    cellInfo.firstRespCells | cellInfo.lateRespCells);

fprintf('Mean/Err SI: %s/%s\n',...
    num2str(round(mean(imgStats.allRespCellsSI,2),2,'significant')),...
    num2str(round(ste(imgStats.allRespCellsSI,2),2,'significant')))

if strcmp(ds,'FSAV_attentionV1')
    [~,imgStats.respCellsSITestVsNaive] = kstest2(imgStats.allRespCellsSI, ...
        imgStats_naive.allRespCellsSI,'tail','smaller');
    fprintf('Kolmogorov-Smirnov Test SI, behavior vs. naive, p = %s\n',...
        num2str(round(imgStats.respCellsSITestVsNaive,2,'significant')))
end

save([fnout 'imgStats'],'imgStats')

%% plot orientation analysis (Figure 3)
% adaptAnalysis = struct;
% ind = cellInfo.firstRespCells & cellInfo.minRespCells;
% adaptAnalysis.firstRespCells = cellfun(@(x) mean(x(ind)),antiAnalysis.adapt);
% adaptAnalysis.firstRespCellsErr = cellfun(@(x) ste(x(ind),2),antiAnalysis.adapt);
% [~,adaptAnalysis.testAV] = ttest(antiAnalysis.adapt{visualTrials}(ind),...
%     antiAnalysis.adapt{auditoryTrials}(ind),'tail','right');
% ind2 = ind & cellInfo.isTuned;
% adaptAnalysis.oriGroupsTest = anova1(antiAnalysis.adapt{3}(ind2),...
%     cellInfo.oriPref(ind2),'off');

% get first, late win, and late cycle response for passive ori tuning groups
oriGroups = struct;
oriGroups.n = nan(1,nOri);
oriGroups.nShortCyc = nan(1,nOri);
oriGroups.nTarOrDist = nan(1,nOri);
oriGroups.nTar = nan(1,nOri);
oriGroups.lateCycSI = nan(1,nOri);
oriGroups.lateCycSIErr = nan(1,nOri);
oriGroups.firstResp = nan(2,nOri);
oriGroups.firstRespErr = nan(2,nOri);
oriGroups.lateWin = nan(2,nOri);
oriGroups.lateWinErr = nan(2,nOri);
oriGroups.lateCycResp = nan(2,nOri);
oriGroups.lateCycRespErr = nan(2,nOri);
oriGroups.firstTC = cell(2,nOri);
oriGroups.firstTCErr = cell(2,nOri);
oriGroups.longTC = cell(2,nOri);
oriGroups.longTCErr = cell(2,nOri);
oriGroups.cycTC = cell(2,nOri);
oriGroups.cycTCErr = cell(2,nOri);
oriGroups.lateCycTest = nan(1,nOri);
oriGroups.lateWinTest = nan(1,nOri);
oriGroups.targetTC = cell(1,nOri);
oriGroups.targetTCErr = cell(1,nOri);
oriGroups.targetTuningResp = cell(1,nOri);
oriGroups.targetTuningRespErr = cell(1,nOri);
oriGroups.targetTuningStim = cell(1,nOri);
oriGroups.targetTuningStimErr = cell(1,nOri);
oriGroups.firstStimRespForTargetAnalysis = nan(1,nOri);
oriGroups.firstStimRespErrForTargetAnalysis = nan(1,nOri);
% adaptAnalysis.oriGroupsN = nan(1,nOri);
% adaptAnalysis.oriGroups = nan(3,nOri);
% adaptAnalysis.oriGroupsErr = nan(3,nOri);
for iori = 1:nOri
    ind = cellInfo.isTuned & cellInfo.oriPref == iori &...
        (cellInfo.lateCycRespCells);
%     ind2 = cellInfo.isTuned & cellInfo.oriPref == iori & ...
%         cellInfo.firstRespCells & cellInfo.minRespCells;
    ind_tarAndDist = cellInfo.isTuned & cellInfo.oriPref == iori &...
        (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    ind_tar = cellInfo.isTuned & cellInfo.oriPref == iori & cellInfo.targetRespCells;
    oriGroups.n(iori) = sum(ind);
    oriGroups.nShortCyc(iori) = sum(ind & cellInfo.isShortCycExpt);
    oriGroups.nTarOrDist(iori) = sum(ind_tarAndDist);
    oriGroups.nTar(iori) = sum(ind_tar);
%     adaptAnalysis.oriGroupsN(iori) = sum(ind2);
    oriGroups.lateCycSI(iori) = mean(antiAnalysis.lateCycSI(ind));
    oriGroups.lateCycSIErr(iori) = ste(antiAnalysis.lateCycSI(ind),2);
    [~,oriGroups.lateWinTest(iori)] = ttest(mean(antiAnalysis.longTC{visualTrials}...
        (lateWinFr,ind & cellInfo.isShortCycExpt),1),...
        mean(antiAnalysis.longTC{auditoryTrials}...
        (lateWinFr,ind & cellInfo.isShortCycExpt),1));
    [~,oriGroups.lateCycTest(iori)] = ttest(mean(antiAnalysis.lateCycTC{visualTrials}...
        (respwin,ind),1),mean(antiAnalysis.lateCycTC{auditoryTrials}...
        (respwin,ind),1));
    
    oriGroups.targetTC{iori} = cell2mat(cellfun(@(x) mean(x(:,ind_tar),2),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0)');
    oriGroups.targetTCErr{iori} = cell2mat(cellfun(@(x) ste(x(:,ind_tar),2),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0)');
    
    oriGroups.targetTuningResp{iori} = cat(1,mean(mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind_tarAndDist),1)),cellfun(@(x) ...
        mean(mean(x(respwin_target,ind_tarAndDist),1) - mean(x(basewin_0,ind_tarAndDist),1)),...
        targetAnalysis.tc(1:2,visualTrials)));
    oriGroups.targetTuningRespErr{iori} = cat(1,ste(mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind_tarAndDist),1),2),cellfun(@(x) ...
        ste(mean(x(respwin_target,ind_tarAndDist),1) - mean(x(basewin_0,ind_tarAndDist),1),2),...
        targetAnalysis.tc(1:2,visualTrials)));
    oriGroups.firstStimRespForTargetAnalysis(iori) = mean(mean(...
            antiAnalysis.firstCycTC{visualTrials}(respwin,ind_tarAndDist),1));
    oriGroups.firstStimRespErrForTargetAnalysis(iori) = ste(mean(...
            antiAnalysis.firstCycTC{visualTrials}(respwin,ind_tarAndDist),1),2);
    oriGroups.targetTuningStim{iori} = cat(1,0,...
        mean(targetAnalysis.targets{visualTrials}(:,ind_tarAndDist),2));
    oriGroups.targetTuningStimErr{iori} = cat(1,0,...
        ste(targetAnalysis.targets{visualTrials}(:,ind_tarAndDist),2));
        
    for iav = 1:2
        oriGroups.firstResp(iav,iori) = mean(mean(...
            antiAnalysis.firstCycTC{iav}(respwin,ind),1));
        oriGroups.firstRespErr(iav,iori) = ste(mean(...
            antiAnalysis.firstCycTC{iav}(respwin,ind),1),2);
        
        oriGroups.firstTC(iav,iori) = {mean(...
            antiAnalysis.firstCycTC{iav}(:,ind),2)};
        oriGroups.firstTCErr(iav,iori) = {ste(...
            antiAnalysis.firstCycTC{iav}(:,ind),2)};
        
        oriGroups.lateWin(iav,iori) = mean(mean(...
            antiAnalysis.longTC{iav}(lateWinFr,ind),1));
        oriGroups.lateWinErr(iav,iori) = ste(mean(...
            antiAnalysis.longTC{iav}(lateWinFr,ind),1),2);
        
        oriGroups.lateCycResp(iav,iori) = mean(mean(...
            antiAnalysis.lateCycTC{iav}(respwin,ind),1));
        oriGroups.lateCycRespErr(iav,iori) = ste(mean(...
            antiAnalysis.lateCycTC{iav}(respwin,ind),1),2);
        
        oriGroups.longTC{iav,iori} = mean(...
            antiAnalysis.longTC{iav}(:,ind & cellInfo.isShortCycExpt),2);
        oriGroups.longTCErr{iav,iori} = ste(...
            antiAnalysis.longTC{iav}(:,ind & cellInfo.isShortCycExpt),2);
        
        oriGroups.cycTC{iav,iori} = mean(...
            antiAnalysis.lateCycTC{iav}(:,ind),2);
        oriGroups.cycTCErr{iav,iori} = ste(...
            antiAnalysis.lateCycTC{iav}(:,ind),2);
%         adaptAnalysis.oriGroups(iav,iori) = mean(...
%             antiAnalysis.adapt{iav}(ind2),2);
%         adaptAnalysis.oriGroupsErr(iav,iori) = ste(...
%             antiAnalysis.adapt{iav}(ind2),2);
%         if iav == 2
%             adaptAnalysis.oriGroups(3,iori) = mean(...
%                 antiAnalysis.adapt{3}(ind2),2);
%             adaptAnalysis.oriGroupsErr(3,iori) = ste(...
%                 antiAnalysis.adapt{3}(ind2),2);
%         end
        
        
    end
end
ind = cellInfo.isTuned &...
    (cellInfo.lateCycRespCells);
lateCycRespDiffData = ...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1) - ...
    mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,ind),1);
oriGroups.lateCycRespTest = anova1(lateCycRespDiffData,cellInfo.oriPref(ind),'off');
oriGroups.lateCycSITest = anova1(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind),'off');

figure
suptitle('Tuned, First or Late Distractor Responsive Neurons')
for iori = 1:nOri
    subplot(2,2,iori)
    for iav = 1:2
        y = oriGroups.cycTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.cycTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',cycTCLim)  
    vline(respWinTT,'k--')
    hline(0,'k:')
    figAxForm
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
    
end

print([fnout 'tuningAnalysis_oriGroupsLateCycTC'],'-dpdf','-fillpage')    


figure
suptitle('First Cyc; Tuned, First or Late Distractor Responsive Neurons')
for iori = 1:nOri
    subplot(2,2,iori)
    for iav = 1:2
        y = oriGroups.firstTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.firstTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',firstTCLim)  
    vline(respWinTT,'k--')
    hline(0,'k:')
    figAxForm
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
end
print([fnout 'tuningAnalysis_oriGroupsFirstTC'],'-dpdf','-fillpage')  



setFigParams4Print('portrait')
figure;
suptitle('Tuned, First or Late Distractor Responsive Neurons, Short Cyc Expts')
for iori = 1:nOri
    subplot(4,2,iori)
    for iav = 1:2
        y = oriGroups.longTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.longTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
    figYAxis([],'dF/F',oriTCLim)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nShortCyc(iori)),num2str(round(oriGroups.lateWinTest(iori),2,'significant'))))
    
%     subplot(4,2,iori)
%     for iav = 1:2
%         y = oriGroups.cycTC{iav,iori}(tcStartFrame:end);
%         yerr = oriGroups.cycTCErr{iav,iori}(tcStartFrame:end);
%         hold on
%         shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
%     end
%     figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
%     figYAxis([],'dF/F',oriTCLim)  
%     vline(respWinTT,'k--')
%     hline(0,'k:')
%     figAxForm
%     title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
%         num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
    subplot(4,2,iori+4)
    for itar = 1:2
        y = oriGroups.targetTC{iori}(tcStartFrame:end,itar);
        yerr = oriGroups.targetTCErr{iori}(tcStartFrame:end,itar);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,hiLoColor{itar});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    title(sprintf('Pref %s, n=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nTar(iori))))
    
end

print([fnout 'tuningAnalysis_oriGroupsTC'],'-dpdf','-fillpage')    

setFigParams4Print('portrait')
figure;
suptitle(sprintf('Tuned, Late Distractor-Responsive Neurons, %s/%s',...
    num2str(sum(cellInfo.isTuned & (cellInfo.lateCycRespCells))),...
    num2str(length(cellInfo.isTuned))))
xsub = [-0.25 +0.25];
subplot 421
h = bar(oriGroups.firstResp','group','EdgeColor','none','BarWidth',1);
legend({'Visual','Auditory'})
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title('First Stim')
subplot 422
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.firstResp(iav,:),...
        oriGroups.firstRespErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title('First Stim')

subplot 423
h = bar(oriGroups.lateWin','group','EdgeColor','none','BarWidth',1);
legend({'Visual','Auditory'})
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_win)
figAxForm
title('Late Window')
subplot 424
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.lateWin(iav,:),...
        oriGroups.lateWinErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_win)
figAxForm
title('Late Window')

subplot 425
h = bar(oriGroups.lateCycResp','group','EdgeColor','none','BarWidth',1);
legend({'Visual','Auditory'})
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title('Late Stim Resp')
subplot 426
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.lateCycResp(iav,:),...
        oriGroups.lateCycRespErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycRespTest,2,'significant'))))

subplot 427
h = bar(oriGroups.n,'EdgeColor','none','BarWidth',0.5);
hold on
for iori = 1:nOri
    text(iori,oriGroups.n(iori)+1,num2str(oriGroups.n(iori)))
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriNLim)
figAxForm

subplot 427
h = bar(oriGroups.n,'EdgeColor','none','BarWidth',0.5);
hold on
for iori = 1:nOri
    text(iori,oriGroups.n(iori)+1,num2str(oriGroups.n(iori)))
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'N Cells',oriNLim)
figAxForm

subplot 428
h = bar(oriGroups.lateCycSI,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:nOri,oriGroups.lateCycSI,oriGroups.lateCycSIErr,'.')
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'Selectivity',siOriLim)
figAxForm
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycSITest,2,'significant'))))

print([fnout 'tuningAnalysis_lateDistRespCells'],'-dpdf','-fillpage')

% task tuning
setFigParams4Print('landscape')    
figure
suptitle(sprintf('Tuned and Distractor or Target Resp. Neurons %s/%s',...
    num2str(sum(cellInfo.isTuned & ...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells))),...
    num2str(length(cellInfo.isTuned))))
for iori = 1:nOri
    subplot(2,2,iori)
    x = oriGroups.targetTuningStim{iori};
    xerr = oriGroups.targetTuningStimErr{iori};
    y = oriGroups.targetTuningResp{iori};
    yerr = oriGroups.targetTuningRespErr{iori};
    errorbar(0,oriGroups.firstStimRespForTargetAnalysis(iori),...
        oriGroups.firstStimRespErrForTargetAnalysis(iori),'.')
    hold on
    errorbar(x,y,yerr,yerr,xerr,xerr,'.')
    figXAxis([],'Task Ori. (deg)',[-10 100],x,oriBins([1,3,4]))
    figYAxis([],'dF/F',oriLim_taskResp)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Pref %s, n=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nTarOrDist(iori))))

end

print([fnout 'tuningAnalysis_oriGroupsTaskTuning'],'-dpdf','-fillpage')   
%%
aurocGroups = struct;
aurocGroups.name = {'Dist.';'Tar';'n.d.'};
aurocGroups.cmp(1).name = 'first:target';
aurocGroups.cmp(2).name = 'late:target';
for icmp = 1:2
    aurocGroups.cmp(icmp).n = nan(1,3);
    aurocGroups.cmp(icmp).n4SI = nan(1,3);
    aurocGroups.cmp(icmp).targetTC = cell(1,3);
    aurocGroups.cmp(icmp).targetTCErr = cell(1,3); 
    aurocGroups.cmp(icmp).targetTC_distOnly = cell(1,3);
    aurocGroups.cmp(icmp).targetTCErr_distOnly = cell(1,3);    
    aurocGroups.cmp(icmp).firstTC = cell(1,3);
    aurocGroups.cmp(icmp).firstTCErr = cell(1,3);
    aurocGroups.cmp(icmp).lateCycSI = nan(1,3);
    aurocGroups.cmp(icmp).lateCycSIErr = nan(1,3);
    aurocGroups.cmp(icmp).lateCycSI_temp = nan(1,2);
    aurocGroups.cmp(icmp).lateCycSIErr_temp = nan(1,2);
end
for igrp = 1:3
    if igrp == 3
        ind_first = ~cellInfo.firstStimAuROCTest & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
        ind_late = ~cellInfo.lateStimAuROCTest & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    elseif igrp == 1
        ind_first = cellInfo.firstStimAuROCTest & cellInfo.firstStimAuROC < 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
        ind_late = cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC < 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    elseif igrp == 2
        ind_first = cellInfo.firstStimAuROCTest & cellInfo.firstStimAuROC > 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
        ind_late = cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC > 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    end
    aurocGroups.cmp(1).n(igrp) = sum(ind_first);
    aurocGroups.cmp(2).n(igrp) = sum(ind_late);
    aurocGroups.cmp(1).n4SI(igrp) = sum(ind_first &...
        (cellInfo.lateCycRespCells));
    aurocGroups.cmp(2).n4SI(igrp) = sum(ind_late &...
        (cellInfo.lateCycRespCells));
    aurocGroups.cmp(1).targetTC{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first),2);
    aurocGroups.cmp(1).targetTCErr{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first),2);
    aurocGroups.cmp(2).targetTC{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late),2);
    aurocGroups.cmp(2).targetTCErr{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late),2);
    aurocGroups.cmp(1).targetTC_distOnly{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first & (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).targetTCErr_distOnly{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first & (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).targetTC_distOnly{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late & (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).targetTCErr_distOnly{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late & (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).firstTC{igrp} = mean(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_first &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).firstTCErr{igrp} = ste(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_first &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).firstTC{igrp} = mean(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).firstTCErr{igrp} = ste(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).lateCycSI(igrp) = mean(antiAnalysis.lateCycSI(ind_first &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).lateCycSIErr(igrp) = ste(antiAnalysis.lateCycSI(ind_first &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).lateCycSI(igrp) = mean(antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).lateCycSIErr(igrp) = ste(antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.lateCycRespCells)),2);
end
ind = (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
[~,sortInd] = sort(cellInfo.firstStimAuROC(ind));
tc = targetAnalysis.tc{3,visualTrials}(:,ind)';
aurocGroups.cmp(1).hm = flipud(tc(sortInd,:));
[~,sortInd] = sort(cellInfo.lateStimAuROC(ind));
tc = targetAnalysis.tc{3,visualTrials}(:,ind)';
aurocGroups.cmp(2).hm = flipud(tc(sortInd,:));

setFigParams4Print('portrait')
figure
suptitle('Target Aligned Resp., All Task Resp. Cells')
colormap(brewermap([],'*RdBu'))
for icmp = 1:2
    subplot(3,2,icmp)
    imagesc(aurocGroups.cmp(1).hm)
    caxis(hmLim)
    colorbar
    figXAxis([],'Time (ms)',[],ttLabelFr_target,ttLabel_target)
    figYAxis([],'Cell # (auROC Sorted)',[])
    figAxForm
    title(aurocGroups.cmp(icmp).name)
    subplot(3,2,icmp+2)
    L = [];
    for igrp = 1:3
        y = aurocGroups.cmp(icmp).targetTC{igrp};
        yerr = aurocGroups.cmp(icmp).targetTCErr{igrp};
        hold on
        h = shadedErrorBar_chooseColor(tt_targetTC,y,yerr,aurocColor{igrp});
        L(igrp) = h.mainLine;
    end
    figXAxis([],'Time (ms)',[tt_targetTC(1) tt_targetTC(end)],...
        ttLabel_target,ttLabel_target)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    vline(preTargetStimLabel,'k:')
    figAxForm
    title(sprintf('%s, n=%s',aurocGroups.cmp(icmp).name,...
        num2str(sum(aurocGroups.cmp(icmp).n))))
    legend(L,strcat(aurocGroups.name,repmat(' n=',[3,1]),...
        cellfun(@num2str,num2cell(aurocGroups.cmp(icmp).n)','unif',0)))

    subplot(3,2,icmp+4)
    for igrp = 1:3
        y = aurocGroups.cmp(icmp).firstTC{igrp}(tcStartFrame:end);
        yerr = aurocGroups.cmp(icmp).firstTCErr{igrp}(tcStartFrame:end);
        hold on
        h = shadedErrorBar_chooseColor(tt_cycTC,y,yerr,aurocColor{igrp});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    title(sprintf('%s, First Stim Resp.',aurocGroups.cmp(icmp).name))
end
print([fnout 'tuningAnalysis_aurocGroupsTC'],'-dpdf','-fillpage')  

setFigParams4Print('landscape')
figure
suptitle('Distractor Resp. Cells')
for icmp = 1:2
    subplot(2,2,icmp)
    L = [];
    for igrp = 1:3
        y = aurocGroups.cmp(icmp).targetTC_distOnly{igrp};
        yerr = aurocGroups.cmp(icmp).targetTCErr_distOnly{igrp};
        hold on
        h = shadedErrorBar_chooseColor(tt_targetTC,y,yerr,aurocColor{igrp});
        L(igrp) = h.mainLine;
    end
    figXAxis([],'Time (ms)',[tt_targetTC(1) tt_targetTC(end)],...
        ttLabel_target,ttLabel_target)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    vline(preTargetStimLabel,'k:')
    figAxForm
    title(sprintf('%s, n=%s',aurocGroups.cmp(icmp).name,...
        num2str(sum(aurocGroups.cmp(icmp).n4SI))))
    legend(L,strcat(aurocGroups.name,repmat(' n=',[3,1]),...
        cellfun(@num2str,num2cell(aurocGroups.cmp(icmp).n4SI)','unif',0)))
    
    subplot(2,2,icmp+2)
    bar(aurocGroups.cmp(icmp).lateCycSI,'EdgeColor','none','BarWidth',0.5);
    hold on
    errorbar(1:3,aurocGroups.cmp(icmp).lateCycSI,aurocGroups.cmp(icmp).lateCycSIErr,'.')
    figXAxis([],'Pref. Task Stim. (by auROC)',[0 4],1:3,aurocGroups.name)
    figYAxis([],'Selectivity',siOriLim)
    figAxForm
    for igrp = 1:3
        text(igrp,aurocGroups.cmp(icmp).lateCycSI(igrp)+.2,num2str(aurocGroups.cmp(icmp).n4SI(igrp)));
    end
    title(aurocGroups.cmp(icmp).name)
end
print([fnout 'tuningAnalysis_aurocGroupsSI'],'-dpdf','-fillpage')
%%
taskGroups = struct;
taskGroups.type = {'Dist.';'Hard Tar';'Easy Tar';'n.p.'};
taskGroups.name = 'One-Way ANOVA Across Task Stim';
taskGroups.n = nan(1,4);
taskGroups.n4SI = nan(1,4);
taskGroups.targetTC = cell(1,4);
taskGroups.targetTCErr = cell(1,4);    
taskGroups.firstTC = cell(1,4);
taskGroups.firstTCErr = cell(1,4);
taskGroups.lateCycSI = nan(1,4);
taskGroups.lateCycSIErr = nan(1,4);
for igrp = 1:4
    if igrp == 4
        ind = ~cellInfo.taskTuningTest & ...
            (cellInfo.lateCycRespCells ...
            | cellInfo.targetRespCells);
    else
        ind = cellInfo.taskTuningTest & cellInfo.taskTuningPref == igrp & ...
            (cellInfo.lateCycRespCells ...
            | cellInfo.targetRespCells);
    end
    taskGroups.n(igrp) = sum(ind);
    taskGroups.n4SI(igrp) = sum(ind &...
        (cellInfo.lateCycRespCells));
    taskGroups.targetTC{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind),2);
    taskGroups.targetTCErr{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind),2);
    taskGroups.firstTC{igrp} = mean(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind),2);
    taskGroups.firstTCErr{igrp} = ste(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind),2);
    taskGroups.lateCycSI(igrp) = mean(antiAnalysis.lateCycSI(ind &...
        (cellInfo.lateCycRespCells)),2);
    taskGroups.lateCycSIErr(igrp) = ste(antiAnalysis.lateCycSI(ind &...
        (cellInfo.lateCycRespCells)),2);
end

setFigParams4Print('portrait')
figure
suptitle('Target Aligned Resp., All Task Resp. Cells')
colormap(brewermap([],'*RdBu'))

subplot(2,2,1)
L = [];
for igrp = 1:4
    y = taskGroups.targetTC{igrp};
    yerr = taskGroups.targetTCErr{igrp};
    hold on
    h = shadedErrorBar_chooseColor(tt_targetTC,y,yerr,taskTuneColor{igrp});
    L(igrp) = h.mainLine;
end
figXAxis([],'Time (ms)',[tt_targetTC(1) tt_targetTC(end)],...
    ttLabel_target,ttLabel_target)
figYAxis([],'dF/F',targetTCLim)  
hline(0,'k:')
vline(respWinTT,'k--')
vline(preTargetStimLabel,'k:')
figAxForm
title(sprintf('Task Stim Pref. Group, n=%s',num2str(sum(taskGroups.n))))
legend(L,strcat(taskGroups.type,repmat(' n=',[4,1]),...
    cellfun(@num2str,num2cell(taskGroups.n)','unif',0)))

subplot(2,2,2)
for igrp = 1:4
    y = taskGroups.firstTC{igrp}(tcStartFrame:end);
    yerr = taskGroups.firstTCErr{igrp}(tcStartFrame:end);
    hold on
    h = shadedErrorBar_chooseColor(tt_cycTC,y,yerr,taskTuneColor{igrp});
end
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',targetTCLim)  
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title('Task Stim Pref. Groups, First Stim Resp.')
subplot(2,2,3)
bar(taskGroups.lateCycSI,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:4,taskGroups.lateCycSI,taskGroups.lateCycSIErr,'.')
figXAxis([],'Pref. Task Stim. (by ANOVA)',[0 5],1:4,taskGroups.type)
figYAxis([],'Selectivity',siOriLim)
figAxForm
for igrp = 1:4
    text(igrp,taskGroups.lateCycSI(igrp)+.2,num2str(taskGroups.n4SI(igrp)));
end
title('Dist. Resp. Cells')

print([fnout 'tuningAnalysis_taskTuningGroups'],'-dpdf')
%%

respBinEdges = [-1, -0.002:0.002:0.006, 1];
nTotal = nan(1,4);
figure
subplot 131
% ind = (cellInfo.lateCycRespCells) & cellInfo.targetRespCells;
ind = (cellInfo.lateCycRespCells);
ind2 = ind & cellInfo.firstStimAuROC < 0.5 & cellInfo.firstStimAuROCTest;
nTotal(1) = sum(ind2);
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1);
y = antiAnalysis.lateCycSI(ind2);
plot(x,y,'o')
c1 = corrcoef(x,y);
subplot 132
x = binnedMean(mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges,minTrN);
[y, yerr] = binnedMean(antiAnalysis.lateCycSI(ind2),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges,minTrN);
errorbar(x,y,yerr,'.-')
subplot 133
n = histcounts(...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges);
plot(x,n,'.-')
subplot 131
hold on
ind2 = ind & cellInfo.firstStimAuROC > 0.5 & cellInfo.firstStimAuROCTest;
nTotal(2) = sum(ind2);
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1);
y = antiAnalysis.lateCycSI(ind2);
plot(x,y,'o')
c2 = corrcoef(x,y);
subplot 132
hold on
x = binnedMean(mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges,minTrN);
[y, yerr] = binnedMean(antiAnalysis.lateCycSI(ind2),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges,minTrN);
errorbar(x,y,yerr,'.-')
subplot 133
hold on
n = histcounts(...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges);
plot(x,n,'.-')
subplot 131
hold on
ind2 = ind & ~cellInfo.firstStimAuROCTest;
nTotal(3) = sum(ind2);
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1);
y = antiAnalysis.lateCycSI(ind2);
plot(x,y,'o')
c3 = corrcoef(x,y);
subplot 132
hold on
x = binnedMean(mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges,minTrN);
[y, yerr] = binnedMean(antiAnalysis.lateCycSI(ind2),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges,minTrN);
errorbar(x,y,yerr,'.-')
subplot 133
hold on
n = histcounts(...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges);
plot(x,n,'.-')
subplot 132
x = binnedMean(mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),respBinEdges,minTrN);
[y, yerr] = binnedMean(antiAnalysis.lateCycSI(ind),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),respBinEdges,minTrN);
errorbar(x,y,yerr,'k.-')
subplot 133
hold on
n = histcounts(...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),respBinEdges);
plot(x,n,'k.-')
nTotal(4) = sum(ind);
subplot 131
figXAxis([],'Late Stim Resp by Task Pref.',[])
figYAxis([],'Selectivity',[])
figAxForm
title('Dist. Resp. Neurons')
legend(strcat({'Dist; corr=';'Tar; corr=';'NP; corr='},...
    {num2str(c1(1,2));num2str(c2(1,2));num2str(c3(1,2))}),'location','southoutside')
subplot 132
figXAxis([],'Late Stim Resp by Task Pref.',[])
figYAxis([],'Selectivity',[])
figAxForm
legend({'Dist';'Tar';'NP;';'All'},'location','southoutside')
subplot 133
figXAxis([],'Late Stim Resp by Task Pref.',[])
figYAxis([],'N Cells',[])
figAxForm
legend(strcat({'Dist; n=';'Tar; n=';'NP; n=';'All; n='},...
    cellfun(@num2str,num2cell(nTotal),'unif',0)'),'location','southoutside')

print([fnout 'tuningAnalysis_taskTuningGroups_SIbyResp'],'-dpdf','-fillpage')

respBinEdges = [-1, -0.002:0.002:0.006, 1];
nTotal = nan(1,5);
cOri = nan(1,4);
ind = (cellInfo.lateCycRespCells) & cellInfo.isTuned;
figure
for iori = 1:nOri
    subplot 131
    hold on
    ind2 = ind & cellInfo.oriPref == iori;
    nTotal(1) = sum(ind2);
    x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1);
    y = antiAnalysis.lateCycSI(ind2);
    plot(x,y,'o')
    c = corrcoef(x,y);
    subplot 132
    hold on
    x = binnedMean(mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),...
        mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges,minTrN);
    [y, yerr] = binnedMean(antiAnalysis.lateCycSI(ind2),...
        mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges,minTrN);
    errorbar(x,y,yerr,'.-')
    subplot 133
    hold on
%     n = histcounts(discretize(...
%         mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges));

    n = histcounts(...
        mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind2),1),respBinEdges);
    plot(x,n,'.-')
    nTotal(iori) = sum(ind2);
    cOri(iori) = c(1,2);
end
subplot 132
hold on
x = binnedMean(mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),respBinEdges,minTrN);
[y, yerr] = binnedMean(antiAnalysis.lateCycSI(ind),...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),respBinEdges,minTrN);
errorbar(x,y,yerr,'k.-')
subplot 133
hold on
n = histcounts(...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),respBinEdges);
plot(x,n,'k.-')

nTotal(5) = sum(ind);
subplot 131
figXAxis([],'Late Stim Resp by Task Pref.',[])
figYAxis([],'Selectivity',[])
figAxForm
title('Dist. Resp. Neurons')
legend(strcat(cellfun(@num2str,num2cell(orientations),'unif',0)',repmat('; corr=',[4,1]),...
    cellfun(@num2str,num2cell(cOri),'unif',0)'),'location','southoutside')
subplot 132
figXAxis([],'Late Stim Resp by Task Pref.',[])
figYAxis([],'Selectivity',[])
figAxForm
legend(cellfun(@num2str,num2cell(orientations),'unif',0)','location','southoutside')
subplot 133
figXAxis([],'Late Stim Resp by Task Pref.',[])
figYAxis([],'N Cells',[])
figAxForm
legend(strcat([cellfun(@num2str,num2cell(orientations),'unif',0)';'All'],...
    repmat({'; n='},[5,1]),...
    cellfun(@num2str,num2cell(nTotal),'unif',0)'),'location','southoutside')

print([fnout 'tuningAnalysis_oriGroups_SIbyResp'],'-dpdf','-fillpage')


%%
imgStats.cellGroups(1).name = 'ori';
imgStats.cellGroups(2).name = 'task';


ind = cellInfo.lateCycRespCells;
auROCpref = nan(1,length(ind));
auROCpref(cellInfo.firstStimAuROC < 0.5 & cellInfo.firstStimAuROCTest) = 1;
auROCpref(cellInfo.firstStimAuROC > 0.5 & cellInfo.firstStimAuROCTest) = 2;
auROCpref(~cellInfo.firstStimAuROCTest) = 3;
imgStats.cellGroups(2).test = anova1(antiAnalysis.lateCycSI(ind),...
    auROCpref(ind));

%% decoding analysis (Figure 4)

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));

targetRespAll = mean(targetAnalysis.tc{allTrialsInd,visualTrials}...
    (respwin_target,:),1);
[~,targetSortInd] = sort(targetRespAll);
subplot 131
% lateCycRespAll = mean(antiAnalysis.lateCycTC{allTrialsInd}(respwin,:),1);
% [~,lateCycSortInd] = sort(lateCycRespAll);
hm = flipud(antiAnalysis.lateCycTC{visualTrials}...
    (tcStartFrame:cycTCEndFr,targetSortInd)');
imagesc(hm)
caxis(hmLim)
colorbar
figXAxis([],'Time (ms)',[],ttLabelFr_cyc,ttLabel_cyc)
figYAxis([],'Cell # (Target Resp Sorted)',[])
figAxForm
title('Late Dist. Resp, All Cells')

subplot 132
% targetRespAll = mean(targetAnalysis.tc{allTrialsInd,visualTrials}...
%     (respwin_target,:),1);
hm = flipud(targetAnalysis.tc{1,visualTrials}...
    ((tcStartFrame-1):(cycTCEndFr-1),targetSortInd)');
imagesc(hm)
caxis(hmLim)
colorbar
figXAxis([],'Time (ms)',[],ttLabelFr_cyc,ttLabel_cyc)
figYAxis([],'Cell # (Target Resp Sorted)',[])
figAxForm
title('Hard Target Resp, All Cells')

subplot 133
hm = flipud(targetAnalysis.tc{2,visualTrials}...
    ((tcStartFrame-1):(cycTCEndFr-1),targetSortInd)');
imagesc(hm)
caxis(hmLim)
colorbar
figXAxis([],'Time (ms)',[],ttLabelFr_cyc-1,ttLabel_cyc)
figYAxis([],'Cell # (Target Resp Sorted)',[])
figAxForm
title('Easy Target Resp, All Cells')

print([fnout 'heatmapVisTrialsAllCells_dist&target'],'-dpdf','-fillpage')

figure
colormap gray
suptitle('Responsive Cells IDs, sorted by target response')

subplot 131
ind = cellInfo.lateCycRespCells;
indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
imagesc(flipud(indBW(targetSortInd)))
title(sprintf('Late Cyc. Resp. Cells n=%s',num2str(sum(ind))))
figYAxis([],'Cell # (Target Resp Sorted)',[])
figAxForm

subplot 132
ind = cellInfo.targetRespCells;
indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
imagesc(flipud(indBW(targetSortInd)))
title(sprintf('Target Resp. Cells n=%s',num2str(sum(ind))))
figYAxis([],'Cell # (Target Resp Sorted)',[])
figAxForm

print([fnout 'cellIDforHeatmap_dist&target'],'-dpdf','-fillpage')

figure
subplot 121
ind = strcmp(trOutTypeName,'H-All');
x = tt_cycTC;
y = nanmean(targetAnalysis.trOutTC{ind,visualTrials}(tcStartFrame:end,cellInfo.lateCycRespCells),2);
yerr = ste(targetAnalysis.trOutTC{ind,visualTrials}(tcStartFrame:end,cellInfo.lateCycRespCells),2);
shadedErrorBar_chooseColor(x,y,yerr,trOutColor{strcmp(trOutType,'h')});
hold on
ind = strcmp(trOutTypeName,'M-All');
y = nanmean(targetAnalysis.trOutTC{ind,visualTrials}(tcStartFrame:end,cellInfo.lateCycRespCells),2);
yerr = ste(targetAnalysis.trOutTC{ind,visualTrials}(tcStartFrame:end,cellInfo.lateCycRespCells),2);
shadedErrorBar_chooseColor(x,y,yerr,trOutColor{strcmp(trOutType,'m')});
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',outTCLim)  
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
subplot 122
ind = strcmp(trOutTypeName,'FA');
x = tt_cycTC;
y = nanmean(targetAnalysis.trOutTC{ind,visualTrials}(tcStartFrame:end,cellInfo.lateCycRespCells),2);
yerr = ste(targetAnalysis.trOutTC{ind,visualTrials}(tcStartFrame:end,cellInfo.lateCycRespCells),2);
shadedErrorBar_chooseColor(x,y,yerr,trOutColor{strcmp(trOutType,'fa')});
hold on
ind = strcmp(trOutTypeName,'CR');
y = nanmean(targetAnalysis.trOutTC{ind,visualTrials}(tcStartFrame:end,cellInfo.lateCycRespCells),2);
yerr = ste(targetAnalysis.trOutTC{ind,visualTrials}(tcStartFrame:end,cellInfo.lateCycRespCells),2);
shadedErrorBar_chooseColor(x,y,yerr,trOutColor{strcmp(trOutType,'cr')});
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',outTCLim)  
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm

dcModel = struct;
dcModel(1).cellTypeName = {'Dist. Only';'Dist. & Tar.'};
dcModel(1).name = 'Target';
dcModel(2).name = 'Detect';
cellTypeIDAllCells = sum(cat(2,(cellInfo.lateCycRespCells & ~cellInfo.targetRespCells).*1,...
    (cellInfo.lateCycRespCells & cellInfo.targetRespCells).*2,...
    (~cellInfo.lateCycRespCells & cellInfo.targetRespCells).*3),2);
dcModel(1).cellTypeID = cellTypeIDAllCells(cellInfo.dcModelCells);
for imod = 1:2
    dcModel(imod).av(1).name = 'Visual';
    dcModel(imod).av(2).name = 'Auditory';
    for iav = 1:2
        dcModel(imod).av(iav).correlation = [];
        dcModel(imod).av(iav).weight = [];
        dcModel(imod).av(iav).pctCorrect = nan(4,nexp);
        dcModel(imod).av(iav).testPerformance = nan(4,nexp);
        dcModel(imod).av(iav).pctCorrect_movRespWin = nan(nMovWin,nexp);
        dcModel(imod).av(iav).pctCorrect_otherAV = nan(1,nexp);
    end
end
for iav = 1:2
    for iexp = 1:nexp
        d = decodeAnalysis(iexp).av(iav);
        
        dcModel(1).av(iav).correlation = cat(2, ...
            dcModel(1).av(iav).correlation,d.correlationTarget);
        dcModel(1).av(iav).weight = cat(1, ...
            dcModel(1).av(iav).weight,d.weightTarget);
        dcModel(1).av(iav).pctCorrect(:,iexp) = cat(2,...
            d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
        dcModel(1).av(iav).testPerformance(:,iexp) = cat(2,...
            d.pctCorrectXStimTarget_train,d.pctCorrectAllTarget_train) - ...
            cat(2,...
            d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
        dcModel(1).av(iav).pctCorrect_movRespWin(:,iexp) = d.pctCorrectTargetMovRespWin_holdout;
        dcModel(1).av(iav).pctCorrect_otherAV(iexp) = d.pctCorrectTarget_otherAV;
        
        dcModel(2).av(iav).correlation = cat(2, ...
            dcModel(2).av(iav).correlation,d.correlationDetect);
        dcModel(2).av(iav).weight = cat(1, ...
            dcModel(2).av(iav).weight,d.weightDetect);
        dcModel(2).av(iav).pctCorrect(:,iexp) = cat(2,...
            d.pctCorrectXStimDetect_holdout,d.pctCorrectAllDetect_holdout);
        dcModel(2).av(iav).testPerformance(:,iexp) = cat(2,...
            d.pctCorrectXStimDetect_train,d.pctCorrectAllDetect_train) - ...
            cat(2,...
            d.pctCorrectXStimDetect_holdout,d.pctCorrectAllDetect_holdout);
        dcModel(2).av(iav).pctCorrect_movRespWin(:,iexp) = d.pctCorrectDetectMovRespWin_holdout;
        dcModel(2).av(iav).pctCorrect_otherAV(iexp) = d.pctCorrectDetect_otherAV;
    end
end

stimLabel = {'0';'HT';'ET';'All'};
mdlXStimAlpha = 0.05./3;

figure;
for imod = 1:2
    if imod == 1
        iplot = 0;
    else
        iplot = 2;
    end
    for iav = 1:2
        y = dcModel(imod).av(iav).pctCorrect;
        yerr = ste(dcModel(imod).av(iav).pctCorrect,2);
        stimInd = sum(~isnan(y),2) > minTrN;
        y = y(stimInd,:);
        yerr = yerr(stimInd,:);
        x = 1:sum(stimInd);
        [mdlTest,mdlP] = ttest(y,0.5,'dim',2,'alpha',mdlXStimAlpha);
        subplot(2,2,iplot+iav)
        for iexp = 1:nexp
%             ind = ~isnan(y(:,iexp));
            hold on
            plot(1:length(x),y(:,iexp),'k.-');
        end
        for itest = 1:length(x)
            if mdlTest(itest) == 1 
                text(x(itest),1,sprintf('*%s',...
                    num2str(round(mdlP(itest),2,'significant'))))
            else
                text(x(itest),1,sprintf('%s',...
                    num2str(round(mdlP(itest),2,'significant'))))
            end
        end
        errorbar(1:length(x),nanmean(y,2),yerr,'.-','MarkerSize',15)
        figXAxis([],'Stim',[0 length(x)+1],x,stimLabel(stimInd))
        hline(0.5,'k:')
        figYAxis([],'Fraction Correct',[0 1])
        figAxForm
        title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
            dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
    end
end

print([fnout 'decodeModelPctCorrectXStim'],'-dpdf','-fillpage')

figure;
suptitle('Performace of trained model on hold-out trials')
for imod = 1:2
    if imod == 1
        iplot = 0;
    else
        iplot = 2;
    end
    for iav = 1:2
        y = dcModel(imod).av(iav).testPerformance;
        yerr = ste(dcModel(imod).av(iav).testPerformance,2);
        stimInd = sum(~isnan(y),2) > minTrN;
        y = y(stimInd,:);
        yerr = yerr(stimInd,:);
        x = 1:sum(stimInd);
        [mdlTest,mdlP] = ttest(y,0,'dim',2,'alpha',mdlXStimAlpha);
        subplot(2,2,iplot+iav)
        for iexp = 1:nexp
            ind = ~isnan(y(:,iexp));
            hold on
            plot(x(ind),y(ind,iexp),'k.-');
        end
        for itest = 1:length(x)
            if mdlTest(itest) == 1
                text(x(itest),0.5,sprintf('*%s',...
                    num2str(round(mdlP(itest),2,'significant'))))
            else
                text(x(itest),0.5,sprintf('%s',...
                    num2str(round(mdlP(itest),2,'significant'))))
            end
        end
        errorbar(x,nanmean(y,2),yerr,'.-','MarkerSize',15)
        figXAxis([],'Stim',[0 length(x)+1],1:length(x),stimLabel(stimInd))
        hline(0,'k:')
        figYAxis([],'Train - HO Frac. Correct',[-0.5 0.5])
        figAxForm
        title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
            dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
    end
end

print([fnout 'decodeModelPerformanceXStim'],'-dpdf','-fillpage')

figure
suptitle('Percent Correct of responses analyzed at different time windows around stimulus; red is used analysis window; blue is min RT')

winInd = [1:3,5:nwins];
x = movWinLabelMs(winInd);
for imod = 1:2
    if imod == 1
        iplot = 0;
    else
        iplot = 2;
    end
    for iav = 1:2
        y = dcModel(imod).av(iav).pctCorrect_movRespWin(winInd,:);
        yerr = ste(dcModel(imod).av(iav).pctCorrect_movRespWin(winInd,:),2);
%         [mdlTest,mdlP] = ttest(y,0,'dim',2,'alpha',mdlXStimAlpha);
        subplot(2,2,iplot+iav)
        for iexp = 1:nexp
            ind = ~isnan(y(:,iexp));
            hold on
            plot(x(ind),y(ind,iexp),'k.-');
        end
        errorbar(x,nanmean(y,2),yerr,'.-','MarkerSize',15)
        figXAxis([],'Resp. Win. Center (Relative to Stim)',[x(1) x(end)],x(1:2:length(x)),round(x(1:2:length(x))))
        hline(0.5,'k:')
        figYAxis([],'Frac. Correct',[0 1])
        figAxForm
        title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
            dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
        vline(respWinTT,'r--')
        vline(baseWinTT,'k--')
        vline(0,'k:')
        if iav == 1
            vline(visRTwindow(1),'b--')
        else
            vline(audRTwindow(1),'b--')
        end
    end
end

print([fnout 'decodeModelPerformanceXRespWindow'],'-dpdf','-fillpage')

for iav = 1:2
    figure
    suptitle(sprintf('%s,Task correlation and weights of cells used in model',...
        dcModel(1).av(iav).name))
    
    weightInd = dcModel(imod).av(iav).weight > 0;
    binnedTargetWeight = nan(1,2);
    binnedTargetWeightErr = nan(1,2);
    binnedDetectWeight = nan(1,2);
    binnedDetectWeightErr = nan(1,2);
    binnedDetectWeightTest = nan(1,2);
    binnedTargetWeight(1) = mean(dcModel(1).av(iav).weight(~weightInd));
    binnedTargetWeight(2) = mean(dcModel(1).av(iav).weight(weightInd));
    binnedTargetWeightErr(1) = ste(dcModel(1).av(iav).weight(~weightInd),1);
    binnedTargetWeightErr(2) = ste(dcModel(1).av(iav).weight(weightInd),1);
    binnedDetectWeight(1) = mean(dcModel(2).av(iav).weight(~weightInd));
    binnedDetectWeight(2) = mean(dcModel(2).av(iav).weight(weightInd));
    binnedDetectWeightErr(1) = ste(dcModel(2).av(iav).weight(~weightInd),1);
    binnedDetectWeightErr(2) = ste(dcModel(2).av(iav).weight(weightInd),1);
    [~,binnedDetectWeightTest(1)] = ttest(dcModel(2).av(iav).weight(~weightInd));
    [~,binnedDetectWeightTest(2)] = ttest(dcModel(2).av(iav).weight(weightInd));
    
    for imod = 1:2
        x = dcModel(imod).av(iav).correlation;
        y = dcModel(imod).av(iav).weight;

        subplot(2,2,imod)
        plot(x,y,'.')
        figXAxis([],'Correlation',[-1 1])
        figYAxis([],'Weight',weightLim)
        figAxForm
        vline(0,'k:')
        hline(0,'k:')
        title(sprintf('%s Model',dcModel(imod).name))
    end
    
    subplot(2,2,3)
    plot(dcModel(1).av(iav).weight,dcModel(2).av(iav).weight,'.')
    figXAxis([],'Target Weight',weightLim)
    figYAxis([],'Detect Weight',weightLim)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    
    subplot(2,2,4)
    hold on
    for i = 1:2
        if i == 1
            ind = ~weightInd;
        elseif i == 2
            ind = weightInd;
        end
        y = dcModel(2).av(iav).weight(ind);
        x = binnedTargetWeight(i).*(ones(sum(ind),1));
        plot(x,y,'.')
        text(binnedTargetWeight(i),weightLim(end),num2str(round(...
            binnedDetectWeightTest(i),2,'significant')))
    end
    errorbar(binnedTargetWeight,binnedDetectWeight,binnedDetectWeightErr,'.-')
    figXAxis([],'Target Weight',binnedWeightLim)
    figYAxis([],'Detect Weight',weightLim)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    
    dcModel(1).av(iav).binnedWeight = binnedTargetWeight;
    dcModel(1).av(iav).binnedWeightErr = binnedTargetWeightErr;
    dcModel(2).av(iav).binnedWeight = binnedDetectWeight;
    dcModel(2).av(iav).binnedWeightErr = binnedDetectWeightErr;
    dcModel(2).av(iav).binnedWeightTest = binnedDetectWeightTest;
    
    print([fnout 'decodeModelWeights_' dcModel(1).av(iav).name],'-dpdf','-fillpage')
end

figure
for imod = 1:2
    for iav = 1:2
        if imod == 1
            subplot(2,2,iav)
        elseif imod == 2
            subplot(2,2,iav+2)
        end
        if iav == visualTrials
            otherAV = 'Auditory';
        else
            otherAV = 'Visual';
        end
        x = 1:2;
        y = [dcModel(imod).av(iav).pctCorrect(end,:)', dcModel(imod).av(iav).pctCorrect_otherAV'];
        hold on
        plot(x,y,'k')
        errorbar(x,mean(y,1),ste(y,1),'.')
        [~,p] = ttest(y,0.5);
        [~,p2] = ttest(y(:,1),y(:,2));
        
        figXAxis([],'Test Trials',[0 3],x,{dcModel(imod).av(iav).name,otherAV})
        figYAxis([],'Fraction Correct',[0 1])
        figAxForm
        hline(0.5,'k:')
        title(sprintf('%s Model, Train %s',dcModel(imod).name,dcModel(imod).av(iav).name))
        for i = 1:2
            if p(i) < 0.05/7
                text(i,1,sprintf('*%s',num2str(round(p(i),2,'significant'))))
            else
                text(i,1,num2str(round(p(i),2,'significant')))
            end
        end
        if p2 < 0.05/4
            text(1.5,0,sprintf('*%s',num2str(round(p2,2,'significant'))))
        else
            text(1.5,0,sprintf('%s',num2str(round(p2,2,'significant'))))
        end
    end
end

print([fnout 'decodeModelPctCorrectOtherAVModel'],'-dpdf','-fillpage')

figure
for imod = 1:2
    subplot(1,2,imod)
    x = dcModel(imod).av(visualTrials).weight;
    y = dcModel(imod).av(auditoryTrials).weight;
    hold on
    plot(x,y,'.')
    mdl = fitlm(x,y);
    rsq = mdl.Rsquared.Ordinary;
    yfit = predict(mdl,x);
    plot(x,yfit,'-')
    figXAxis([],'Visual Weight',weightLim)
    figYAxis([],'Auditory Weight',weightLim)
    figAxForm
    title(sprintf('%s Model, Rsq = %s',dcModel(imod).name,num2str(round(rsq,2,'significant'))))    
end

print([fnout 'decodeModelCompareAVWeights'],'-dpdf','-fillpage')

%% Figure 5, Effect of attention on weights

figure;
suptitle('Target Model Attention - Late Resp. Cells')
for iav = 1:2
    subplot(2,2,iav)
    ind = cellInfo.lateCycRespCells;
    ind1 = cellInfo.lateCycRespCells & ~cellInfo.targetRespCells;
    ind2 = cellInfo.lateCycRespCells & cellInfo.targetRespCells;
    
    x = antiAnalysis.lateCycSI(ind1);
    y = cellInfo.targetWeight{iav}(ind1);    
    plot(x,y,'.')
    hold on
    x = antiAnalysis.lateCycSI(ind2);
    y = cellInfo.targetWeight{iav}(ind2);    
    plot(x,y,'.')
    
    figXAxis([],'Selectivity',siLim)
    figYAxis([],'Weight',weightLim)
    figAxForm
    vline(0,'k:')
    hline(0,'k:')
    title(sprintf('Train %s',dcModel(1).av(iav).name))
    
    subplot(2,2,iav+2)
    
    x = mean(antiAnalysis.lateCycSI(ind1));
    xerr = ste(antiAnalysis.lateCycSI(ind1),2);
    y = nanmean(cellInfo.targetWeight{iav}(ind1));
    yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
    errorbar(x,y,yerr,yerr,xerr,xerr,'.')
    hold on
    x = mean(antiAnalysis.lateCycSI(ind2));
    xerr = ste(antiAnalysis.lateCycSI(ind2),2);
    y = nanmean(cellInfo.targetWeight{iav}(ind2));
    yerr = ste(cellInfo.targetWeight{iav}(ind2),2);    
    errorbar(x,y,yerr,yerr,xerr,xerr,'.')
    
    figXAxis([],'Selectivity',siLimSum)
    figYAxis([],'Weight',weightLimSum)
    figAxForm
    vline(0,'k:')
    hline(0,'k:')
    title(sprintf('Train %s',dcModel(1).av(iav).name))
    legend(dcModel(1).cellTypeName,'location','southwest')
end

print([fnout 'targetModelWeightXAttn'],'-dpdf','-fillpage')

figure;
suptitle('Detect Model Attention - Late Resp. Cells')
for iav = 1:2
    subplot(2,2,iav)
    ind = cellInfo.lateCycRespCells;
    ind1 = cellInfo.lateCycRespCells & ~cellInfo.targetRespCells;
    ind2 = cellInfo.lateCycRespCells & cellInfo.targetRespCells;
    
    x = antiAnalysis.lateCycSI(ind1);
    y = cellInfo.targetWeight{iav}(ind1);    
    plot(x,y,'.')
    hold on
    x = antiAnalysis.lateCycSI(ind2);
    y = cellInfo.detectWeight{iav}(ind2);    
    plot(x,y,'.')
    
    figXAxis([],'Selectivity',siLim)
    figYAxis([],'Weight',weightLim)
    figAxForm
    vline(0,'k:')
    hline(0,'k:')
    title(sprintf('Train %s',dcModel(1).av(iav).name))
    
    subplot(2,2,iav+2)
    
    x = mean(antiAnalysis.lateCycSI(ind1));
    xerr = ste(antiAnalysis.lateCycSI(ind1),2);
    y = nanmean(cellInfo.detectWeight{iav}(ind1));
    yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
    errorbar(x,y,yerr,yerr,xerr,xerr,'.')
    hold on
    x = mean(antiAnalysis.lateCycSI(ind2));
    xerr = ste(antiAnalysis.lateCycSI(ind2),2);
    y = nanmean(cellInfo.detectWeight{iav}(ind2));
    yerr = ste(cellInfo.detectWeight{iav}(ind2),2);    
    errorbar(x,y,yerr,yerr,xerr,xerr,'.')
    
    figXAxis([],'Selectivity',siLimSum)
    figYAxis([],'Weight',weightLimSum)
    figAxForm
    vline(0,'k:')
    hline(0,'k:')
    title(sprintf('Train %s',dcModel(1).av(iav).name))
    legend(dcModel(1).cellTypeName,'location','southwest')
end

print([fnout 'detectModelWeightXAttn'],'-dpdf','-fillpage')
