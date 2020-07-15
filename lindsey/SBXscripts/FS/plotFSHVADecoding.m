ds = 'FS_HVA';
eval(ds)
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];
rc = behavConstsAV;
load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_decodeStruct_cells_' ds '.mat']));

g = 7;


area_list = {'LM','AL','PM'};
narea = length(area_list);
eval(ds)
fnout = fullfile(rc.caOutputDir,ds);
nexp = size(decodeAnalysis,2);
iav = 1;
psychfit = struct([]);
psychfit_thresh = cell(1,narea);
psychfit_slope = cell(1,narea);
psychfit_fit = cell(1,narea);
options = struct;
options.sigmoidName = 'norm';
options.expType = 'YesNo';
priorWidth = @(x) (x>=20).*(x<=80);
options.priors{2} = priorWidth;

pctCorrect = cell(1,narea);
for iexp = 1:nexp
    iarea = find(strcmp(area_list,decodeDataExpt(iexp).exptArea));
    [oris x trStimID] = unique(decodeDataExpt(iexp).av(iav).stim);
    nTrials = zeros(1,size(oris,2));
    for i = 1:length(oris)
    	nTrials(1,i) = length(find(trStimID==i));
    end
    psychfit_thresh_temp = nan(g,1);
    pctCorrect_temp = nan(1,g);
    
    for ig = 1:g
        if ~isnan(decodeAnalysis(iexp).cellBin(ig).nCellsSelected)
            holdout = decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectXStimTarget_holdout;
            holdout(1,:) = 1-holdout(1,:);
            nCorrect = round(nTrials' .* mean(holdout,2));
            pctCorrect_temp(1,ig) = sum(nCorrect(:,:),1)./sum(nTrials(:,:),2);
        end
    end
    pctCorrect{iarea} = [pctCorrect{iarea}; pctCorrect_temp];
end
    
figure;
for iarea = 1:narea
    subplot(2,1,1)
    c = defaultPlotColors(iarea);
    plot(repmat(decodeAnalysis(1).cellBins',[1 size(pctCorrect{iarea},1)]), pctCorrect{iarea}', 'Color', c);
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('Fraction Correct')
    ylim([0.5 1])
    subplot(2,1,2)
    c = defaultPlotColors(iarea);
    errorbar(decodeAnalysis(1).cellBins', nanmean(pctCorrect{iarea}',2), nanstd(pctCorrect{iarea}',[],2)./sqrt(sum(~isnan(pctCorrect{iarea}'),2)), 'Color', c);
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('Fraction Correct')
    ylim([0.5 1])
end
print(fullfile(fnout,['allAreas_fractCorrect_noThresh.pdf']), '-dpdf','-fillpage')


for ig = 1:g
    start = [0 0 0];
    figure;
    for iexp = 1:nexp
        iarea = find(strcmp(area_list,decodeDataExpt(iexp).exptArea));
        [oris x trStimID] = unique(decodeDataExpt(iexp).av(iav).stim);
        nTrials = zeros(1,size(oris,2));
        for i = 1:length(oris)
            nTrials(1,i) = length(find(trStimID==i));
        end
        data = zeros(size(decodeAnalysis(iexp).stims,2),3);
        data(:,1) = decodeAnalysis(iexp).stims';
        options.stimulusRange = [min(data(:,1),[],1) max(data(:,1),[],1)];
        if ~isnan(decodeAnalysis(iexp).cellBin(ig).nCellsSelected)
                data(:,2) = round(nTrials' .* mean(decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectXStimTarget_holdout,2));
                data(:,3) = nTrials';
                [max_val max_ind] = max(data(:,2)./data(:,3),[],1);
                if max_ind==4
                    LR = 1-(mean(data(4:5,2)./data(4:5,3),1));
                else
                    LR = 1-(data(end,2)./data(end,3));
                end
                options.fixedPars = [NaN;NaN;LR;(data(1,2)./data(1,3));NaN];
                fitout = psignifit(data,options);
                subplot(3,3,iarea+(narea.*start(iarea)))
                plotPsych(fitout);
                if start(iarea) == 0
                title(decodeDataExpt(iexp).exptArea)
                end
                xlabel('Orientation (deg)')
                ylabel('Fraction Yes')
                start(iarea) = 1+start(iarea);
        end
    end
    suptitle([num2str(decodeAnalysis(iexp).cellBins(ig)) ' cells'])
    print(fullfile(fnout,['allAreas_' num2str(decodeAnalysis(iexp).cellBins(ig)) 'cells_psychfits_noThresh.pdf']), '-dpdf','-fillpage')
end



for iexp = 1:nexp
    iarea = find(strcmp(area_list,decodeDataExpt(iexp).exptArea));
    [oris x trStimID] = unique(decodeDataExpt(iexp).av(iav).stim);
    nTrials = zeros(1,size(oris,2));
    for i = 1:length(oris)
    	nTrials(1,i) = length(find(trStimID==i));
    end
    data = zeros(size(decodeAnalysis(iexp).stims,2),3);
    data(:,1) = decodeAnalysis(iexp).stims';
    options.stimulusRange = [min(data(:,1),[],1) max(data(:,1),[],1)];
    psychfit_thresh_temp = nan(g,1);
    psychfit_slope_temp = nan(g,1);
    psychfit_fit_temp = nan(g,5);
    figure;
    start = 1;
    for ig = 1:g
        if ~isnan(decodeAnalysis(iexp).cellBin(ig).nCellsSelected)
            data(:,2) = round(nTrials' .* mean(decodeAnalysis(iexp).cellBin(ig).av(iav).pctCorrectXStimTarget_holdout,2));
            data(:,3) = nTrials';
            [max_val max_ind] = max(data(:,2)./data(:,3),[],1);
            if max_ind==4
                LR = 1-(mean(data(4:5,2)./data(4:5,3),1));
            else
                LR = 1-(data(end,2)./data(end,3));
            end
            options.fixedPars = [NaN;NaN;LR;(data(1,2)./data(1,3));NaN];
            fitout = psignifit(data,options);
            psychfit_fit_temp(ig,:) = fitout.Fit';
            psychfit_thresh_temp(ig,:) = getThreshold(fitout, 0.5,1);
            psychfit_slope_temp(ig,:) = getSlopePC(fitout, 0.5,1);
            if psychfit_fit_temp(ig,3)>=psychfit_fit_temp(ig,4)
                psychfit_thresh_temp(ig,:) = 90;
                psychfit_slope_temp(ig,:) = NaN;
            end
            subplot(3,3,start)
            plotPsych(fitout);
            title([num2str(decodeAnalysis(iexp).cellBins(ig)) ' cells'])
            xlabel('Orientation (deg)')
            ylabel('Fraction Yes')
            start = 1+start;
        end
    end
    psychfit_thresh{iarea} = [psychfit_thresh{iarea} psychfit_thresh_temp];
    psychfit_slope{iarea} = [psychfit_slope{iarea} psychfit_slope_temp];
    psychfit_fit{iarea} = cat(3,psychfit_fit{iarea}, psychfit_fit_temp);
    suptitle([decodeDataExpt(iexp).exptName '_' decodeDataExpt(iexp).exptArea])
    print(fullfile(fnout,[decodeDataExpt(iexp).exptName '_' decodeDataExpt(iexp).exptArea '_psychfitsnoThresh.pdf']), '-dpdf','-fillpage')
end
    
figure;
for iarea = 1:narea
    c = defaultPlotColors(iarea);
    subplot(2,2,1)
    plot(repmat(decodeAnalysis(1).cellBins',[1 size(psychfit_thresh{iarea},2)]), psychfit_thresh{iarea}, 'Color', c);
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('Threshold (deg)')
    ylim([0 100])
    subplot(2,2,2)
    plot(repmat(decodeAnalysis(1).cellBins',[1 size(psychfit_slope{iarea},2)]), psychfit_slope{iarea}, 'Color', c);
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('Slope (HR/deg)')
    ylim([0 .05])
    subplot(2,2,3)
    LR =  squeeze(psychfit_fit{iarea}(:,3,:));
    plot(repmat(decodeAnalysis(1).cellBins',[1 size(LR,2)]), LR, 'Color', c);
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('Lapse rate')
    ylim([0 .5])
    subplot(2,2,4)
    FR =  squeeze(psychfit_fit{iarea}(:,4,:));
    plot(repmat(decodeAnalysis(1).cellBins',[1 size(FR,2)]), FR, 'Color', c);
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('FA rate')
    ylim([0 .5])
end
suptitle('LM- blue; AL- red; PM- yellow')
print(fullfile(fnout,'allExpt_psychfitSummary_noThresh.pdf'), '-dpdf','-fillpage')

    
figure;
for iarea = 1:narea
    c = defaultPlotColors(iarea);
    subplot(2,2,1)
    errorbar(decodeAnalysis(1).cellBins', nanmean(psychfit_thresh{iarea},2), std(psychfit_thresh{iarea},[],2)./sqrt(size(~isnan(sum(psychfit_thresh{iarea},2)),2)),'Color', c)
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('Threshold (deg)')
    ylim([0 100])
    subplot(2,2,2)
    errorbar(decodeAnalysis(1).cellBins', nanmean(psychfit_slope{iarea},2), std(psychfit_slope{iarea},[],2)./sqrt(size(~isnan(sum(psychfit_slope{iarea},2)),2)),'Color', c)
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('Slope (HR/deg)')
    ylim([0 .05])
    subplot(2,2,3)
    LR =  squeeze(psychfit_fit{iarea}(:,3,:));
    errorbar(decodeAnalysis(1).cellBins', nanmean(LR,2), std(LR,[],2)./sqrt(size(~isnan(sum(LR(:,2),2)),2)),'Color', c)
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('Lapse rate')
    ylim([0 .5])
    subplot(2,2,4)
    FR =  squeeze(psychfit_fit{iarea}(:,4,:));
    errorbar(decodeAnalysis(1).cellBins', nanmean(FR,2), nanstd(FR,[],2)./sqrt(size(~isnan(sum(FR(:,2),2)),2)),'Color', c)
    hold on
    xlabel('# Cells')
    xlim([0 100])
    ylabel('FA rate')
    ylim([0 .5])
end
print(fullfile(fnout,'avg_psychfitSummary_noThresh.pdf'), '-dpdf','-fillpage')
