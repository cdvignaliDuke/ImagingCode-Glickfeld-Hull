function [stimChoiceCorr,trialsUsed,modelWeights,pctCorr_fixedVar,...
    pctCorr_fixedVar_pcsOnly,pctCorr_fixedVarOnly,pctCorr_pcs,modelWeights_pcs] = ...
    getSubSampledModel_boot(fixedPredictorVar,currentModelVar,...
    resp_withFixedVar,nBoot,nSamp)

    % subsample trials and find stimulus:choice correlation for each bootstrap
    stimChoiceCorr = nan(1,nBoot);
    trialsUsed = cell(1,nBoot);
    for i = 1:nBoot
        trialsUsed{i} = randsample(length(currentModelVar),nSamp);
        stimChoiceCorr(i) = corr(currentModelVar(trialsUsed{i}), fixedPredictorVar(trialsUsed{i}));
    end
    
    % run model for each bootstrap
    modelWeights = cell(1,nBoot);
    modelWeights_pcs = cell(1,nBoot);
    pctCorr_pcs = nan(1,nBoot);
    pctCorr_fixedVar = nan(1,nBoot);
    pctCorr_fixedVar_pcsOnly = nan(1,nBoot);
    pctCorr_fixedVarOnly = nan(1,nBoot);
    for i = 1:nBoot
        r = resp_withFixedVar(trialsUsed{i},1:(end-1));
        C = eye(size(r,2));
        [~,~,dc_boot] = glmfit(r*C,currentModelVar(trialsUsed{i}),'binomial');
        modelWeights_pcs{i} = dc_boot.beta;
        pctCorr_pcs(i) = getPctCorr_hoData(r,currentModelVar(trialsUsed{i}),...
            mean(currentModelVar(trialsUsed{i})));
        
        r = resp_withFixedVar(trialsUsed{i},:);
        C = eye(size(r,2));
        [~,~,dc_boot] = glmfit(r*C,currentModelVar(trialsUsed{i}),'binomial');
        modelWeights{i} = dc_boot.beta;
        pctCorr_fixedVar(i) = getPctCorr_hoData(r,currentModelVar(trialsUsed{i}),...
            mean(currentModelVar(trialsUsed{i})));
        pctCorr_fixedVar_pcsOnly(i) = getPctCorr_hoData_choosePCs(...
            r,currentModelVar(trialsUsed{i}),mean(currentModelVar(trialsUsed{i})),1:(size(r,2)-1));
        pctCorr_fixedVarOnly(i) = getPctCorr_hoData(...
            r(:,end),currentModelVar(trialsUsed{i}),mean(currentModelVar(trialsUsed{i})));
    end

end