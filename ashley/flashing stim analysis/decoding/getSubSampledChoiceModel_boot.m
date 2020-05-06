function [stimChoiceCorr,trialsUsed,choiceModelWeights,pctCorr_choice_fixedStim,...
    pctCorr_choice_fixedStim_pcsOnly,pctCorr_choice_fixedStimOnly] = ...
    getSubSampledChoiceModel_boot(stimVar,choiceVar,resp_withFixedStim,nBoot,nSamp)

    % subsample trials and find stimulus:choice correlation for each bootstrap
    stimChoiceCorr = nan(1,nBoot);
    trialsUsed = cell(1,nBoot);
    for i = 1:nBoot
        trialsUsed{i} = randsample(length(choiceVar),nSamp);
        stimChoiceCorr(i) = corr(choiceVar(trialsUsed{i}), stimVar(trialsUsed{i}));
    end
    
    % run model for each bootstrap
    choiceModelWeights = cell(1,nBoot);
    pctCorr_choice_fixedStim = nan(1,nBoot);
    pctCorr_choice_fixedStim_pcsOnly = nan(1,nBoot);
    pctCorr_choice_fixedStimOnly = nan(1,nBoot);
    for i = 1:nBoot
        r = resp_withFixedStim(trialsUsed{i},:);
        C = eye(size(r,2));
        [~,~,dc_boot] = glmfit(r*C,choiceVar(trialsUsed{i}),'binomial');
        choiceModelWeights{i} = dc_boot.beta;
        pctCorr_choice_fixedStim(i) = getPctCorr_hoData(r,choiceVar(trialsUsed{i}),...
            mean(choiceVar(trialsUsed{i})));
        pctCorr_choice_fixedStim_pcsOnly(i) = getPctCorr_hoData_choosePCs(...
            r,choiceVar(trialsUsed{i}),mean(choiceVar(trialsUsed{i})),1:(size(r,2)-1));
        pctCorr_choice_fixedStimOnly(i) = getPctCorr_hoData(...
            r(:,end),choiceVar(trialsUsed{i}),mean(choiceVar(trialsUsed{i})));
    end

end