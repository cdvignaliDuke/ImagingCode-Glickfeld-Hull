function learningSummaryWF(mouse)

rc = behavConstsWF(mouse);
eval(['i' mouse '_paths'])
load(fullfile(rc.structOutput, [mouse '_dSStruct.mat']))
load(fullfile(rc.structOutput, [mouse '_' roi_date '_roi_masks.mat']))
roi = struct;
for id = 1:length(dS)
    ntrials = length(dS(id).trialOutcomeCell);
    successIx = find(strcmp(dS(id).trialOutcomeCell,'success'));
    failureIx = find(strcmp(dS(id).trialOutcomeCell,'failure'));
    missedIx = find(strcmp(dS(id).trialOutcomeCell,'ignore'));
    dS(id).F = mean(mean(dS(id).pressAlign(10:40,:,successIx),1),3);
    dS(id).pressAlignDF = bsxfun(@minus, dS(id).pressAlign, dS(id).F);
    dS(id).pressAlignDFoverF = bsxfun(@rdivide, dS(id).pressAlignDF, dS(id).F);
    dS(id).targetAlignDF = bsxfun(@minus, dS(id).targetAlign, dS(id).F);
    dS(id).targetAlignDFoverF = bsxfun(@rdivide, dS(id).targetAlignDF, dS(id).F);
    dS(id).releaseAlignDF = bsxfun(@minus, dS(id).releaseAlign, dS(id).F);
    dS(id).releaseAlignDFoverF = bsxfun(@rdivide, dS(id).releaseAlignDF, dS(id).F);
    
    nROI = size(dS(id).F,2);
    figure
    for iR = 1:nROI
        subplot(3,3,iR)
        plot(mean(dS(id).pressAlignDFoverF(:,iR,successIx),3), 'k')
        hold on
        plot(mean(dS(id).pressAlignDFoverF(:,iR,failureIx),3), 'r')
        hold on
        plot(mean(dS(id).pressAlignDFoverF(:,iR,missedIx),3), 'b')
        title(area_list(iR,:))
        vline(dS(id).eventBufferFrames)
        hline(0)
        xlim([0 dS(id).eventBufferFrames*2])
        roi(iR).pressAlignSIx(:,id) = squeeze(mean(dS(id).pressAlignDFoverF(:,iR,successIx),3));
        roi(iR).pressAlignFIx(:,id) = squeeze(mean(dS(id).pressAlignDFoverF(:,iR,failureIx),3));
        roi(iR).pressAlignMIx(:,id) = squeeze(mean(dS(id).pressAlignDFoverF(:,iR,missedIx),3));
    end
    alignYaxes(gcf)
    suptitle(['Align to press: ' mouse ' ' dS(id).date])
    if ~exist(fullfile(rc.structOutput, dS(id).date))
        mkdir(fullfile(rc.structOutput, dS(id).date))
    end
    print(fullfile(rc.structOutput, dS(id).date, [mouse '_' dS(id).date '_pressAlignROIs.pdf']), '-dpdf')
    
    figure
    for iR = 1:nROI
        subplot(3,3,iR)
        plot(mean(dS(id).releaseAlignDFoverF(:,iR,successIx),3), 'k')
        hold on
        plot(mean(dS(id).releaseAlignDFoverF(:,iR,failureIx),3), 'r')
        title(area_list(iR,:))
        vline(dS(id).eventBufferFrames)
        hline(0)
        xlim([0 dS(id).eventBufferFrames*2])
        roi(iR).releaseAlignSIx(:,id) = squeeze(mean(dS(id).releaseAlignDFoverF(:,iR,successIx),3));
        roi(iR).releaseAlignFIx(:,id) = squeeze(mean(dS(id).releaseAlignDFoverF(:,iR,failureIx),3));
    end
    alignYaxes(gcf)
    suptitle(['Align to release: ' mouse ' ' dS(id).date])
    print(fullfile(rc.structOutput, dS(id).date, [mouse '_' dS(id).date '_releaseAlignROIs.pdf']), '-dpdf')
    
    figure
    for iR = 1:nROI
        subplot(3,3,iR)
        plot(mean(dS(id).targetAlignDFoverF(:,iR,successIx),3), 'k')
        hold on
        plot(mean(dS(id).targetAlignDFoverF(:,iR,missedIx),3), 'b')
        title(area_list(iR,:))
        vline(dS(id).eventBufferFrames)
        hline(0)
        xlim([0 dS(id).eventBufferFrames*2])
        roi(iR).targetAlignSIx(:,id) = squeeze(mean(dS(id).targetAlignDFoverF(:,iR,successIx),3));
        roi(iR).targetAlignMIx(:,id) = squeeze(mean(dS(id).targetAlignDFoverF(:,iR,missedIx),3));
    end
    alignYaxes(gcf)
    suptitle(['Align to target: ' mouse ' ' dS(id).date])
    print(fullfile(rc.structOutput, dS(id).date, [mouse '_' dS(id).date '_targetAlignROIs.pdf']), '-dpdf')
end

