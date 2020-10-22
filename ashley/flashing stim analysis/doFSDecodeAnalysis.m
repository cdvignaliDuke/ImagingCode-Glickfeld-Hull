function decodeAnalysis = doFSDecodeAnalysis(datasetStr,cellsOrDendrites,doAuditory,doBehavior,tuningThresh,nCellBins,oriBinSize,train_stim,doMatch,doHalfHoldout);
%datasetStr: name of dataset
%cellsOrDendrites: 1 == cells; 2 == dendrites
%doAuditory: 1 = auditory trials
%doBehavior: 1 = behaving mouse
%tuningThresh: size of 90% confidence interval for ori est
%nCellBins: number of groups with 100 max and 8 min cells
%oriBinSize: 0 = use all stimuli; >0 sets bin size
%train_stim: list of stimuli used to train decoder (eg: [0 90])
%doMatch: 1 = match for trial number
ds = datasetStr;

rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
tuningReliabilityThresh_decode = tuningThresh;
eval(ds)

mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

if expt(1).folder(1:5) == 'ashle'
    titleStr = ds(6:end);
    if cellsOrDendrites == 1
        load(fullfile(rc.caOutputDir,ds,...
            [mouse_str '_decodeStruct_cells' ds(5:end) '.mat']));
        fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeAnalysis_cells' ds(5:end) '.mat']); 
    elseif cellsOrDendrites == 2
        load(fullfile(rc.caOutputDir,ds,...
            [mouse_str '_decodeStruct_dend' ds(5:end) '.mat']));
        fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeAnalysis_dend' ds(5:end) '.mat']);
    end
elseif expt(1).folder(1:5) == 'linds'
    titleStr = ds;
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_decodeStruct_cells_' ds '.mat']));
    fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeAnalysis_cells_' ds '_halfHoldout_noThresh.mat']); 
end

if doBehavior
    trOutType = {'h';'m';'fa';'cr'};
    trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
end

nexp = size(expt,2);
decodeAnalysis = struct;
decodeAnalysis.av(visualTrials).name = 'Visual';
if doAuditory
    nav = 2;
    decodeAnalysis.av(auditoryTrials).name = 'Auditory';
else
    nav = 1;
end

if oriBinSize>0
    orientations = 0:oriBinSize:(180-oriBinSize);
    oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
    nOri = length(orientations);
end

if nCellBins > 1;
    maxCellNum = 48;
    minCellNum = 1;
    SPO = (nCellBins-1)./log2(maxCellNum./minCellNum);
    cellBins = zeros(1,nCellBins);
    for i = 1:nCellBins
        cellBins(:,i) = round(maxCellNum/(2.^((i-1)/SPO)));
    end
elseif nCellBins == 1
    cellBins = nCellBins;
    doHalfHoldout = 1;
elseif nCellBins == 0
    cellBins = 15;
end


for iexp = 1:nexp
    rng(0) % cells randomized and trials randomized
    cellIndAll = decodeDataExpt(iexp).tuningReliability' <= tuningReliabilityThresh_decode;
    decodeAnalysis(iexp).nCells = sum(cellIndAll);
    decodeAnalysis(iexp).cellBins = cellBins;
    for ig = 1:nCellBins
        fprintf('Expt %s, hold-out analysis: %s cells\n',num2str(iexp),num2str(cellBins(ig)))
        maxCellN = cellBins(ig);
        if maxCellN > 8
            doPCA = 1;
            nPCs = 8;
        else
            doPCA = 0;
        end
        if sum(cellIndAll)<maxCellN
            decodeAnalysis(iexp).cellBin(ig).nCellsSelected = NaN;
            continue
        end
        decodeAnalysis(iexp).cellBin(ig).nCellsSelected = maxCellN;
        if maxCellN == 1
            nBoot = sum(cellIndAll);
        else
            nBoot = round(100./(maxCellN./decodeAnalysis(iexp).nCells));
        end
        %nBoot = 1;
        decodeAnalysis(iexp).cellBin(ig).nBoot = nBoot;
        for iboot = 1:nBoot
            if rem(iboot,20) == 1
                fprintf('%s/%s bootstraps complete\n',num2str(iboot-1),num2str(nBoot))
            end
            ind = find(cellIndAll);
            if cellBins == 1
                cellSampleID = ind(iboot);
            else
                cellSampleID = randsample(ind,maxCellN);
            end
            cellInd = false(length(cellIndAll),1);
            cellInd(cellSampleID) = true;
            decodeAnalysis(iexp).cellBin(ig).cellInd(:,iboot) = cellInd;
            if nav == 2
                if isfield(decodeDataExpt(iexp).av(visualTrials),'catchOutcome')
                    if ~doPCA
                        respAllCells_av_withInv = zscore(cat(2,...
                        decodeDataExpt(iexp).av(visualTrials).resp,...
                        decodeDataExpt(iexp).av(auditoryTrials).resp,...
                        decodeDataExpt(iexp).av(visualTrials).catchResp)');
                    else
                        respAllCells_av_withInv = cat(2,...
                        decodeDataExpt(iexp).av(visualTrials).resp,...
                        decodeDataExpt(iexp).av(auditoryTrials).resp,...
                        decodeDataExpt(iexp).av(visualTrials).catchResp)';
                    end
                    nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                    naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                    ninv = size(decodeDataExpt(iexp).av(visualTrials).catchResp,2);
                    respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                    respAllCells_aud = respAllCells_av_withInv((nvis+1):(nvis+naud),:);
                    catchResp = respAllCells_av_withInv((end-ninv+1):end,:);
                else
                    if ~doPCA
                        respAllCells_av_withInv = zscore(cat(2,...
                            decodeDataExpt(iexp).av(visualTrials).resp,...
                            decodeDataExpt(iexp).av(auditoryTrials).resp)');
                    else
                        respAllCells_av_withInv = cat(2,...
                            decodeDataExpt(iexp).av(visualTrials).resp,...
                            decodeDataExpt(iexp).av(auditoryTrials).resp)';
                    end
                        nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                        naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                        respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                        respAllCells_aud = respAllCells_av_withInv((nvis+1):end,:);
                end
            else
                if ~doPCA
                    respAllCells_av_withInv = zscore(decodeDataExpt(iexp).av(visualTrials).resp');
                else
                    respAllCells_av_withInv = decodeDataExpt(iexp).av(visualTrials).resp';
                end
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
            end
            for iav = 1:nav
                trOut = decodeDataExpt(iexp).av(iav).outcome;
                stims = unique(decodeDataExpt(iexp).av(iav).stim);
                decodeAnalysis(iexp).stims = stims;
                
                if iav == 1 & oriBinSize>0
                    respAllCells = respAllCells_vis;
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                    respAllCells_train = respAllCells;
                    trStimID_train = trStimID;
                    trOut_train = trOut;
                    ind_train = 1:length(stims);
                    trID_train = 1:length(trOut);
                elseif iav == 1
                    respAllCells = respAllCells_vis;
                    [oris x trStimID] = unique(decodeDataExpt(iexp).av(iav).stim);
                    trStimID = trStimID';
                    if ~isempty(train_stim)
                        ind_train = find(ismember(stims, train_stim));
                    else
                        ind_train = 1:length(stims);
                    end
                    ind = [];
                    for is = 1:length(ind_train)
                    	ind = [ind find(trStimID == ind_train(is))];
                    end
                    trID_train = ind;
                    trStimID_train = trStimID(ind);
                    respAllCells_train = respAllCells(ind,:);
                    trOut_train = trOut(ind);
                elseif iav == 2
                    respAllCells = respAllCells_aud;
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                    respAllCells_train = respAllCells;
                    trStimID_train = trStimID;
                    trOut_train = trOut;
                    trID_train = 1:length(trOut);
                    ind_train = 1:length(stims);
                end
                
                nStimBins_all = length(unique(trStimID));
                nStimBins = length(unique(trStimID_train));
                nStimPerBin_all = histcounts(trStimID);
                nStimPerBin = histcounts(trStimID_train);
                minBinN = min(nStimPerBin_all(nStimPerBin_all >= minTrN_mdl));
                respStimSort_train = cell(1,nStimBins);
                trOutStimSort_train = cell(1,nStimBins);
                respStimSort = cell(1,nStimBins);
                trOutStimSort = cell(1,nStimBins);
                trStimSort = cell(1,nStimBins_all);
                trStimSort_test = cell(1,nStimBins_all);
                trStimSort_train = cell(1,nStimBins_all);
                trInd_train_all = [];
                trInd_train = [];
                trInd_test_all = [];
                for istim = 1:length(stims)
                    ind1 = find(trStimID == istim);
                    trStimSort{istim} = ind1;
                    if find(ind_train==istim)
                        if ~doHalfHoldout
                            trInd_train_all = [trInd_train_all ind1];
                        else
                            ind_temp = ind1(randsample(1:length(ind1),ceil(length(ind1)/2)));
                            trInd_train_all = [trInd_train_all ind_temp];
                            trInd_test_all = [trInd_test_all setdiff(ind1,ind_temp)];
                            trStimSort_train{istim} = ind_temp;
                            trStimSort_test{istim} = setdiff(ind1,ind_temp);
                        end
                        ind2 = trID_train(find(trStimID_train == istim));
                        n = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                        if length(ind2) >= minTrN_mdl
                            if istim == 1
                                matchTrialsInd = [];
                                indSample = randsample(ind2,n);
                                matchTrialsInd = cat(2,matchTrialsInd,indSample);
                            else
                                indSample = randsample(ind2,minBinN);
                                matchTrialsInd = cat(2,matchTrialsInd,indSample);
                            end
                            trInd_train = [trInd_train ind2];
                            respStimSort_train{istim} = respAllCells(indSample,cellInd);
                            trOutStimSort_train{istim} = trOut(indSample);
                        end
                    end
                    respStimSort{istim} = respAllCells(ind1,cellInd);
                    trOutStimSort{istim} = trOut(ind1);
                end
                nMatchedTrials = cumsum(cellfun(@length,trOutStimSort_train(ind_train)));
                for istim = 1:nStimBins
                    if istim == 1
                        stimSortInd = cell(1,nStimBins);
                        stimSortInd{istim} = 1:nMatchedTrials;
                    else
                        stimSortInd{istim} = ...
                            (nMatchedTrials(istim-1)+1):nMatchedTrials(istim);
                    end
                end
                if doMatch
                    trInd_train_all = matchTrialsInd;
                end
                decodeAnalysis(iexp).cellBin(ig).av(iav).respAllCells(:,:,iboot) = respAllCells(trInd_train_all,cellInd);
                decodeAnalysis(iexp).cellBin(ig).av(iav).trOut(:,iboot) = trOut(trInd_train_all)';
                decodeAnalysis(iexp).cellBin(ig).av(iav).nTrPerStim = nStimPerBin_all;
                resp = respAllCells(:,cellInd);
                if doPCA
                    [resp_pca, pca_trial, latent] = pca(resp);
                    decodeAnalysis(iexp).cellBin(ig).av(iav).latent(:,iboot) = latent;
                    decodeAnalysis(iexp).cellBin(ig).av(iav).resp(:,:,iboot) = pca_trial;
                    if nPCs>0
                        resp = pca_trial(:,1:nPCs);
                    else
                        resp = pca_trial;
                    end
                end
                
                
                [detectTrInd_train, targetTrInd_train] = getStimAndBehaviorYs(trOut(trInd_train_all));
                [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut);

                C = eye(size(resp,2));
                
                targetCorr = corr(targetTrInd,resp);

                [~,~,targetGLM] = glmfit(resp(trInd_train_all,:)*C,targetTrInd(trInd_train_all,:),'binomial');
                targetWeight = targetGLM.beta;
                dv_target = mean(targetTrInd(trInd_train_all));
                if ~doHalfHoldout
                    pctCorrectTarget_train = getPctCorr_trainData_LG(targetWeight,resp(trInd_train_all,:),targetTrInd(trInd_train_all,:),dv_target);
                    [pctCorrectTarget_ho] = getPctCorr_hoData_train(resp,targetTrInd,trInd_train_all,trInd_train_all,dv_target);
                else
                    [pctCorrectTarget_train MITarget_train] = getPctCorr_trainData_LG(targetWeight,resp(trInd_train_all,:),targetTrInd(trInd_train_all,:),dv_target);
                    [pctCorrectTarget_ho MITarget_ho] = getPctCorr_trainData_LG(targetWeight,resp(trInd_test_all,:),targetTrInd(trInd_test_all,:),dv_target);
                end                    
                pctCorrTarget_xStim_train = nan(1,nStimBins_all);
                pctCorrTarget_xStim_ho = nan(1,nStimBins_all);
                MITarget_xStim_train = nan(1,nStimBins_all);
                MITarget_xStim_ho = nan(1,nStimBins_all);
                for istim = 1:nStimBins_all
                    if nStimPerBin_all(istim) >= minTrN_mdl
                        if ~doHalfHoldout
                            pctCorrTarget_xStim_ho(istim) = getPctCorr_hoData_train(...
                                resp,targetTrInd,trStimSort{istim},trInd_train_all,dv_target);
                            if find(train_stim == istim)
                                pctCorrTarget_xStim_train(istim) = getPctCorr_trainData_LG(...
                                    targetWeight,resp(trStimSort{istim},:),targetTrInd(trStimSort{istim},:),dv_target);
                            end
                        else
                            [pctCorrTarget_xStim_train(istim) MITarget_xStim_train(istim)] = getPctCorr_trainData_LG(...
                                    targetWeight,resp(trStimSort_train{istim},:),targetTrInd(trStimSort_train{istim},:),dv_target);
                            [pctCorrTarget_xStim_ho(istim) MITarget_xStim_ho(istim)] = getPctCorr_trainData_LG(...
                                    targetWeight,resp(trStimSort_test{istim},:),targetTrInd(trStimSort_test{istim},:),dv_target);
                        end
                    end
                end
                pctCorrTarget_xStim_ho(1) = 1-pctCorrTarget_xStim_ho(1);
                pctCorrTarget_xStim_train(1) = 1-pctCorrTarget_xStim_train(1);
                
                decodeAnalysis(iexp).cellBin(ig).av(iav).dvTarget(:,iboot) = dv_target;
                decodeAnalysis(iexp).cellBin(ig).av(iav).correlationTarget(:,iboot) = targetCorr;
                decodeAnalysis(iexp).cellBin(ig).av(iav).weightTarget(:,iboot) = targetWeight;
                decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectAllTarget_train(:,iboot) = pctCorrectTarget_train;
                decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectAllTarget_holdout(:,iboot) = pctCorrectTarget_ho;
                decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectXStimTarget_train(:,iboot) = pctCorrTarget_xStim_train;
                decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectXStimTarget_holdout(:,iboot) = pctCorrTarget_xStim_ho;
                decodeAnalysis(iexp).cellBin(ig).av(iav).MIAllTarget_train(:,iboot) = MITarget_train;
                decodeAnalysis(iexp).cellBin(ig).av(iav).MIAllTarget_holdout(:,iboot) = MITarget_ho;
                decodeAnalysis(iexp).cellBin(ig).av(iav).MIXStimTarget_train(:,iboot) = MITarget_xStim_train;
                decodeAnalysis(iexp).cellBin(ig).av(iav).MIXStimTarget_holdout(:,iboot) = MITarget_xStim_ho;
                
                if doBehavior
                    detectCorr = corr(detectTrInd,resp);
                    [~,~,detectGLM] = glmfit(resp(trInd_train_all,:)*C,detectTrInd(trInd_train_all,:),'binomial');
                    detectWeight = detectGLM.beta(2:end);
                    dv_detect = mean(detectTrInd(trInd_train_all,:));
                    pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp(trInd_train_all,:),detectTrInd(trInd_train_all,:),dv_detect);
                    pctCorrectDetect_ho = getPctCorr_hoData_train(resp(trInd_train_all,:),detectTrInd(trInd_train_all,:),trInd_train_all,trInd_train_all,dv_detect);
                    pctCorrDetect_xStim_train = nan(1,nStimBins_all);
                    pctCorrDetect_xStim_ho = nan(1,nStimBins_all);
                    for istim = 1:nStimBins_all
                        if nStimPerBin_all(istim) >= minTrN_mdl
                                pctCorrDetect_xStim_ho(istim) = getPctCorr_hoData_train(...
                                    resp,detectTrInd,trStimSort{istim},trInd_train_all,dv_detect);
                                if find(train_stim == istim)
                                    pctCorrDetect_xStim_train(istim) = getPctCorr_trainData(...
                                        detectGLM,resp(trStimSort{istim},:),detectTrInd(trStimSort{istim},:),dv_detect);
                                end
                        end
                    end
                    pctCorrDetect_xStim_ho(1) = 1-pctCorrDetect_xStim_ho(1);
                    pctCorrDetect_xStim_train(1) = 1-pctCorrDetect_xStim_train(1);

                    decodeAnalysis(iexp).cellBin(ig).av(iav).dvDetect(:,iboot) = dv_detect;
                    decodeAnalysis(iexp).cellBin(ig).av(iav).correlationDetect(:,iboot) = detectCorr;
                    decodeAnalysis(iexp).cellBin(ig).av(iav).weightDetect(:,iboot) = detectWeight;
                    decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectAllDetect_train(:,iboot) = pctCorrectDetect_train;
                    decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectAllDetect_holdout(:,iboot) = pctCorrectDetect_ho;
                    decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectXStimDetect_train(:,iboot) = pctCorrDetect_xStim_train;
                    decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectXStimDetect_holdout(:,iboot) = pctCorrDetect_xStim_ho;
                end
            end
        end
    end
end
save(fnout, 'decodeAnalysis');
        
        