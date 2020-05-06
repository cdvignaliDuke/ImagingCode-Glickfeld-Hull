maxPCs = 20;
iexp=2;
rng(0)
cellInd = respCellsExpt(iexp).lateRespCells|respCellsExpt(iexp).lateCycRespCells|respCellsExpt(iexp).targetRespCells;
if sum(cellInd) > maxPCs
    nPCs = maxPCs;
else
    nPCs = sum(cellInd);
end
if isfield(decodeDataExpt(iexp).av(visualTrials),'catchOutcome')
    respAllCells_av_withInv = zscore(cat(2,...
        decodeDataExpt(iexp).av(visualTrials).resp,...
        decodeDataExpt(iexp).av(auditoryTrials).resp,...
        decodeDataExpt(iexp).av(visualTrials).catchResp)');
    nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
    naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
    ninv = size(decodeDataExpt(iexp).av(visualTrials).catchResp,2);
    respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
    respAllCells_aud = respAllCells_av_withInv((nvis+1):(nvis+naud),:);
    catchResp = respAllCells_av_withInv((end-ninv+1):end,:);
    
    [coeffAllCells_vis,scoresAllCells_vis,latentAllCells_vis] = pca(respAllCells_vis(:,cellInd));
    [coeffAllCells_aud,scoresAllCells_aud,latentAllCells_aud] = pca(respAllCells_aud(:,cellInd));
    [coeffAllCells_catch,scoresAllCells_catch,latentAllCells_catch] = pca(catchResp(:,cellInd));
    
else
    respAllCells_av_withInv = zscore(cat(2,...
        decodeDataExpt(iexp).av(visualTrials).resp,...
        decodeDataExpt(iexp).av(auditoryTrials).resp)');
    nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
    naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
    respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
    respAllCells_aud = respAllCells_av_withInv((nvis+1):end,:);
    
    [coeffAllCells_vis,scoresAllCells_vis,latentAllCells_vis] = pca(respAllCells_vis(:,cellInd));
    [coeffAllCells_aud,scoresAllCells_aud,latentAllCells_aud] = pca(respAllCells_aud(:,cellInd));
    
end
for iav = 1:2
    trOut = decodeDataExpt(iexp).av(iav).outcome;
    if iav == 1
        respAllCells = scoresAllCells_vis;
        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
    elseif iav == 2
        respAllCells = scoresAllCells_aud;
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
            respStimSort{istim} = respAllCells(indSample,:);
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
    trSampleIndOther{iav} = matchTrialsInd;
%                 dcInfo(iexp).av(iav).resp = respAllCells(matchTrialsInd,:);
%                 dcInfo(iexp).av(iav).trOut = trOut(matchTrialsInd);
%         matchTrialsInd = 1:size(respAllCells,2);
%         resp = zscore(respAllCells(cellInd,matchTrialsInd)');
    decodeAnalysis(iexp).av(iav).respAllCells = respAllCells(matchTrialsInd,:);
    decodeAnalysis(iexp).av(iav).trOut = trOut(matchTrialsInd);
    resp = respAllCells(matchTrialsInd,1:nPCs);
    
     [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

    detectCorr = corr(detectTrInd,resp);
    targetCorr = corr(targetTrInd,resp);

    C = eye(size(resp,2));
    p=1;
    [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
    [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');

    fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
    detectWeight_scores = detectGLM.beta(2:end);
    targetWeight_scores = targetGLM.beta(2:end);
        
    dv_detect = mean(detectTrInd);
    dv_target = mean(targetTrInd);
    
    pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp,detectTrInd,dv_detect);
    pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);
    
    pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
    pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
 
    detectWeight = coeffAllCells_vis(:,1:nPCs)*detectWeight_scores;
    targetWeight = coeffAllCells_vis(:,1:nPCs)*targetWeight_scores;
end