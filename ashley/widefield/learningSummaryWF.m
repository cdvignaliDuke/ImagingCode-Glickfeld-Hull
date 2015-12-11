function learningSummaryWF(mouse,doPlot)
    rc = behavConstsWF(mouse);
    eval(['i' mouse '_paths'])
    load(fullfile(rc.structOutput, [mouse '_dSStruct.mat']))
    load(fullfile(rc.structOutput, [mouse '_' roi_date '_roi_masks.mat']))
    roi = struct;
    for id = 1:length(dS)
        ntrials = length(dS(id).trialOutcomeCell);
        dS(id).successIx = find(strcmp(dS(id).trialOutcomeCell,'success'));
        dS(id).earlySIx = intersect(dS(id).successIx, find(dS(id).reactTimesMs<250));
        dS(id).goodSIx = intersect(dS(id).successIx, find(dS(id).reactTimesMs>250 & dS(id).reactTimesMs<1500));
        dS(id).lateSIx = intersect(dS(id).successIx, find(dS(id).reactTimesMs>800));
        dS(id).failureIx = find(strcmp(dS(id).trialOutcomeCell,'failure'));
        dS(id).missedIx = find(strcmp(dS(id).trialOutcomeCell,'ignore'));
        dS(id).F = mean(dS(id).pressAlign(10:40,:,:,:),1);
        dS(id).pressAlignDF = bsxfun(@minus, dS(id).pressAlign, dS(id).F);
        dS(id).pressAlignDFoverF = bsxfun(@rdivide, dS(id).pressAlignDF, dS(id).F);
        dS(id).pressAlignDFoverF_npsub = squeeze(dS(id).pressAlignDFoverF(:,:,1,:)-dS(id).pressAlignDFoverF(:,:,2,:));
        dS(id).targetAlignDF = bsxfun(@minus, dS(id).targetAlign, dS(id).F);
        dS(id).targetAlignDFoverF = bsxfun(@rdivide, dS(id).targetAlignDF, dS(id).F);
        dS(id).targetAlignDFoverF_npsub = squeeze(dS(id).targetAlignDFoverF(:,:,1,:)-dS(id).targetAlignDFoverF(:,:,2,:));
        dS(id).releaseAlignDF = bsxfun(@minus, dS(id).releaseAlign, dS(id).F);
        dS(id).releaseAlignDFoverF = bsxfun(@rdivide, dS(id).releaseAlignDF, dS(id).F);
        dS(id).releaseAlignDFoverF_npsub = squeeze(dS(id).releaseAlignDFoverF(:,:,1,:)-dS(id).releaseAlignDFoverF(:,:,2,:));
        
        nROI = size(dS(id).F,2);
        
        for iR = 1:nROI
            roi(iR).F(:,id,:) = squeeze(mean(dS(id).F(1,iR,:,:),4));
            if length(dS(id).successIx>0)
                roi(iR).pressAlignSIx(:,id) = squeeze(mean(dS(id).pressAlignDFoverF_npsub(:,iR,dS(id).successIx),3));
                roi(iR).releaseAlignSIx(:,id) = squeeze(mean(dS(id).releaseAlignDFoverF_npsub(:,iR,dS(id).successIx),3));
                roi(iR).targetAlignSIx(:,id) = squeeze(mean(dS(id).targetAlignDFoverF_npsub(:,iR,dS(id).successIx),3));
            else
                roi(iR).pressAlignSIx(:,id) = nan(size(dS(id).pressAlign,1),1);
                roi(iR).releaseAlignSIx(:,id) = nan(size(dS(id).pressAlign,1),1);
                roi(iR).targetAlignSIx(:,id) = nan(size(dS(id).pressAlign,1),1);
            end
            if length(dS(id).failureIx>0)
                roi(iR).pressAlignFIx(:,id) = squeeze(mean(dS(id).pressAlignDFoverF_npsub(:,iR,dS(id).failureIx),3));
                roi(iR).releaseAlignFIx(:,id) = squeeze(mean(dS(id).releaseAlignDFoverF_npsub(:,iR,dS(id).failureIx),3));
            else
                roi(iR).pressAlignFIx(:,id) = nan(size(dS(id).pressAlign,1),1);
                roi(iR).releaseAlignFIx(:,id) = nan(size(dS(id).pressAlign,1),1);
            end
            if length(dS(id).missedIx>0)    
                roi(iR).pressAlignMIx(:,id) = squeeze(mean(dS(id).pressAlignDFoverF_npsub(:,iR,dS(id).missedIx),3));
                roi(iR).targetAlignMIx(:,id) = squeeze(mean(dS(id).targetAlignDFoverF_npsub(:,iR,dS(id).missedIx),3));
            else
                roi(iR).pressAlignMIx(:,id) = nan(size(dS(id).pressAlign,1),1);
                roi(iR).targetAlignMIx(:,id) = nan(size(dS(id).pressAlign,1),1);
            end
            if length(dS(id).goodSIx>0)
                roi(iR).targetAlignGoodSIx(:,id) = squeeze(mean(dS(id).targetAlignDFoverF_npsub(:,iR,dS(id).goodSIx),3));
            else
                roi(iR).targetAlignGoodSIx(:,id) = nan(size(dS(id).pressAlign,1),1);
            end
            if length(dS(id).earlySIx>0)
                roi(iR).targetAlignEarlySIx(:,id) = squeeze(mean(dS(id).targetAlignDFoverF_npsub(:,iR,dS(id).earlySIx),3));
            else
                roi(iR).targetAlignEarlySIx(:,id) = nan(size(dS(id).pressAlign,1),1);
            end
            if length(dS(id).lateSIx>0)
                roi(iR).targetAlignLateSIx(:,id) = squeeze(mean(dS(id).targetAlignDFoverF_npsub(:,iR,dS(id).lateSIx),3));
            else
                roi(iR).targetAlignLateSIx(:,id) = nan(size(dS(id).pressAlign,1),1);
            end
        end
        
        if doPlot
            figure
            for iR = 1:nROI
                subplot(3,3,iR)
                plot(nanmean(dS(id).pressAlignDFoverF_npsub(:,iR,dS(id).successIx),3), 'k')
                hold on
                plot(nanmean(dS(id).pressAlignDFoverF_npsub(:,iR,dS(id).failureIx),3), 'r')
                hold on
                plot(nanmean(dS(id).pressAlignDFoverF_npsub(:,iR,dS(id).missedIx),3), 'b')
                title(area_list(iR,:))
                vline(dS(id).eventBufferFrames)
                hline(0)
                xlim([0 dS(id).eventBufferFrames*2])
            end
            suptitle(['Align to press: ' mouse ' ' dS(id).date])
            if ~exist(fullfile(rc.structOutput, dS(id).date))
                mkdir(fullfile(rc.structOutput, dS(id).date))
            end
            print(fullfile(rc.structOutput, dS(id).date, [mouse '_' dS(id).date '_pressAlignROIs.pdf']), '-dpdf')

            figure
            for iR = 1:nROI
                subplot(3,3,iR)
                plot(nanmean(dS(id).releaseAlignDFoverF_npsub(:,iR,dS(id).successIx),3), 'k')
                hold on
                plot(nanmean(dS(id).releaseAlignDFoverF_npsub(:,iR,dS(id).failureIx),3), 'r')
                title(area_list(iR,:))
                vline(dS(id).eventBufferFrames)
                hline(0)
                xlim([0 dS(id).eventBufferFrames*2])
            end
            alignYaxes(gcf)
            suptitle(['Align to release: ' mouse ' ' dS(id).date])
            print(fullfile(rc.structOutput, dS(id).date, [mouse '_' dS(id).date '_releaseAlignROIs.pdf']), '-dpdf')

            figure
            for iR = 1:nROI
                subplot(3,3,iR)
                plot(nanmean(dS(id).targetAlignDFoverF_npsub(:,iR,dS(id).successIx),3), 'k')
                hold on
                plot(nanmean(dS(id).targetAlignDFoverF_npsub(:,iR,dS(id).missedIx),3), 'b')
                title(area_list(iR,:))
                vline(dS(id).eventBufferFrames)
                hline(0)
                xlim([0 dS(id).eventBufferFrames*2])
            end
            alignYaxes(gcf)
            suptitle(['Align to target: ' mouse ' ' dS(id).date])
            print(fullfile(rc.structOutput, dS(id).date, [mouse '_' dS(id).date '_targetAlignROIs.pdf']), '-dpdf')
        end
    end
    save(fullfile(rc.structOutput, [mouse '_dSStruct.mat']), 'dS')
    save(fullfile(rc.structOutput, [mouse '_roiStruct.mat']), 'roi')
    
    if doPlot
        %plot learning within area
        cmap = copper(length(dS));
        for iR = 1:nROI
            figure
            tt = (1-dS(1).eventBufferFrames:dS(1).eventBufferFrames)./(dS(1).frameRate./1000);
            for id = 1:length(dS)
                subplot(3,3,1)
                plot(tt, roi(iR).pressAlignSIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                title('Press Align Success')
                subplot(3,3,2)
                plot(tt, roi(iR).pressAlignFIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                title('Press Align Early')
                subplot(3,3,3)
                plot(tt, roi(iR).pressAlignMIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                title('Press Align Missed')
                subplot(3,3,4)
                plot(tt, roi(iR).targetAlignSIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                title('Target Align Success')
                subplot(3,3,5)
                plot(tt, roi(iR).targetAlignMIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                title('Target Align Missed')
                subplot(3,3,6)
                plot(tt, roi(iR).releaseAlignSIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                title('Release Align Success')
                subplot(3,3,7)
                plot(tt, roi(iR).releaseAlignFIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                title('Release Align Early')
            end
            subplot(3,3,8)
            plot(roi(iR).F')
            xlabel('Days')
            ylabel('F')
            title('F over time')
            suptitle(['Mouse ' mouse ' area ' area_list(iR,:)])
            print(fullfile(rc.structOutput, [mouse '_allAlignArea' area_list(iR,:) '.pdf']), '-dpdf')
        end        
                
        %plot learning by kind of success
        for iR = 1:nROI
            figure
            tt = (1-dS(1).eventBufferFrames:dS(1).eventBufferFrames)./(dS(1).frameRate./1000);
            for id = 1:length(dS)
                subplot(1,3,1)
                plot(tt, roi(iR).targetAlignEarlySIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                vline([100 250])
                title('Target Align Early Success')
                subplot(1,3,2)
                plot(tt, roi(iR).targetAlignGoodSIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                vline([250 800])
                title('Target Align Good Success')
                subplot(1,3,3)
                plot(tt, roi(iR).targetAlignLateSIx(:,id), 'Color', cmap(id,:))
                hold on 
                vline(dS(1).eventBufferFrames./(dS(1).frameRate./1000))
                xlabel('Time (ms)')
                ylabel('dF/F')
                vline(1500)
                title('Target Align Late Success')
            end
            suptitle(['Mouse ' mouse ' area ' area_list(iR,:)])
            print(fullfile(rc.structOutput, [mouse '_TargetAlignSIxByTimingArea' area_list(iR,:) '.pdf']), '-dpdf')
        end

        for iR = 1:nROI;
            figure;
            for id = 1:length(dS)
                subplot(4,4,id)
                plot(tt, roi(iR).targetAlignEarlySIx(:,id), 'r')
                hold on;
                plot(tt, roi(iR).targetAlignGoodSIx(:,id), 'k')
                hold on
                plot(tt, roi(iR).targetAlignLateSIx(:,id), 'b')
                hold on
                vline(0)
                ylabel('dF/F')
                xlabel('Time (ms)')
                title(dS(id).date)
            end
            suptitle([mouse ' area ' area_list(iR,:)])
            print(fullfile(rc.structOutput, [mouse '_TargetAlignSIxByTimingByDayArea' area_list(iR,:) '.pdf']), '-dpdf')
        end


        %learning plots
        figure;
        nTrials = zeros(1,length(dS));
        nFA = zeros(1,length(dS));
        nMiss = zeros(1,length(dS));
        nResp = zeros(1,length(dS));
        nGuess = zeros(1,length(dS));
        for id = 1:length(dS)
            subplot(2,2,1)
            h(id) = cdfplot(dS(id).reactTimesMs);
            set(h(id), 'Color', cmap(id,:))
            hold on
            xlabel('Time (ms)')
            title('Reaction over time')
            xlim([-1000 2000])
            subplot(2,2,2)
            nTrials(:,id) = length(dS(id).reactTimesMs);
            nMiss(:,id) = length(find(dS(id).reactTimesMs>800));
            subplot(2,2,3)
            nFA(:,id)= length(find(dS(id).reactTimesMs<250));
            subplot(2,2,4)
            nResp(:,id) = length(find(dS(id).reactTimesMs>250 & dS(id).reactTimesMs<550));
            nGuess(:,id) = length(find(dS(id).reactTimesMs>-250 & dS(id).reactTimesMs<0));
        end
        subplot(2,2,2)
        plot(nMiss./nTrials)
        title('Lapse rate')
        ylim([0 1])
        subplot(2,2,3)
        plot(nFA./nTrials)
        title('False alarm rate')
        ylim([0 1])
        subplot(2,2,4)
        plot(nGuess./nResp);
        title('FA/Correct Ratio')
        ylim([0 1])
        suptitle(['Mouse ' mouse 'reaction times during learning'])
        print(fullfile(rc.structOutput, [mouse '_behaviorProgress.pdf']), '-dpdf')

        
        %quantification
        baseSIx = zeros(nROI, length(dS));
        peakSIx = zeros(nROI, length(dS));
        frameSIx = zeros(nROI, length(dS));
        baseMIx = zeros(nROI, length(dS));
        peakMIx = zeros(nROI, length(dS));
        frameMIx = zeros(nROI, length(dS));
        figure;
        for iR = 1:nROI
            for id = 1:length(dS)
                baseSIx(iR,id) = nanmean(roi(iR).targetAlignSIx(dS(id).eventBufferFrames./2:dS(id).eventBufferFrames, id),1);
                [peakSIx(iR,id), frameSIx(iR,id)]= max(roi(iR).targetAlignSIx(dS(id).eventBufferFrames:2.*dS(id).eventBufferFrames, id),[],1);
                baseMIx(iR,id) = nanmean(roi(iR).targetAlignMIx(dS(id).eventBufferFrames./2:dS(id).eventBufferFrames, id),1);
                [peakMIx(iR,id), frameMIx(iR,id)]= max(roi(iR).targetAlignMIx(dS(id).eventBufferFrames:2.*dS(id).eventBufferFrames, id),[],1);
            end
            deltaPeakSIx(iR,:) = peakSIx(iR,:)'-baseSIx(iR,:)';
            subplot(3,3,iR)
            plot(baseSIx(iR,:), 'k')
            hold on;
            plot(deltaPeakSIx(iR,:), 'm')
            title(area_list(iR,:))
        end
        suptitle('Base (black) and Resp (magenta) over days')
        print(fullfile(rc.structOutput, [mouse '_ampChangeByArea.pdf']), '-dpdf')
        figure;
        subplot(1,3,1)
        plot(baseSIx')
        title('Base')
        subplot(1,3,2)
        plot(deltaPeakSIx')
        title('Peak Change')
        subplot(1,3,3)
        plot(frameSIx'./(dS(1).frameRate./1000))
        title('Peak latency (ms)')
        print(fullfile(rc.structOutput, [mouse '_ampChangeAllArea.pdf']), '-dpdf')

        %comparison of response across areas across days
        tt = (1-dS(1).eventBufferFrames:dS(1).eventBufferFrames)./(dS(1).frameRate./1000);
        for id = 1:length(dS)
            figure;
            start = 1;
            for iR1 = 1:nROI
                for iR2 = 1:nROI
                    if iR1 >= iR2
                        continue
                    end
                    subplot(6,5,start)
                    plot(tt, roi(iR1).targetAlignSIx(:,id), 'k')
                    hold on
                    plot(tt, roi(iR2).targetAlignSIx(:,id), 'b')
                    hold on
                    vline(0)
                    axis off
                    title([area_list(iR1,:) '; ' area_list(iR2,:)])
                    start = start+1;
                end
            end
            suptitle([mouse ' ' dS(id).date ' target align'])
            print(fullfile(rc.structOutput, dS(id).date, [mouse '_' dS(id).date '_allAreaTCs_targetAlignSIx.pdf']), '-dpdf')
        end
    
    
        % plot changes in dF/F for roi, np, and np subtraction
        cmap = copper(length(dS));
        nROI = size(dS(1).F,2);
        for iR = 1:nROI
            figure
            tt = (1-dS(1).eventBufferFrames:dS(1).eventBufferFrames)./(dS(1).frameRate./1000);
            for id = 1:length(dS)
                subplot(3,3,1)
                plot(tt, mean(dS(id).pressAlignDFoverF(:,iR,1,dS(id).successIx),4), 'Color', cmap(id,:))
                hold on
                title('press- ROI')
                ylim([-0.03 0.1])
                subplot(3,3,2)
                plot(tt, mean(dS(id).pressAlignDFoverF(:,iR,2,dS(id).successIx),4), 'Color', cmap(id,:))
                hold on
                title('press- surround')
                ylim([-0.03 0.1])
                subplot(3,3,3)
                plot(tt, mean(dS(id).pressAlignDFoverF_npsub(:,iR,dS(id).successIx),3), 'Color', cmap(id,:))
                hold on
                title('press- NP sub')
                ylim([-0.03 0.1])
                subplot(3,3,4)
                plot(tt, mean(dS(id).targetAlignDFoverF(:,iR,1,dS(id).successIx),4), 'Color', cmap(id,:))
                hold on
                title('target- ROI')
                ylim([-0.03 0.1])
                subplot(3,3,5)
                plot(tt, mean(dS(id).targetAlignDFoverF(:,iR,2,dS(id).successIx),4), 'Color', cmap(id,:))
                hold on
                title('target- surround')
                ylim([-0.03 0.1])
                subplot(3,3,6)
                plot(tt, mean(dS(id).targetAlignDFoverF_npsub(:,iR,dS(id).successIx),3), 'Color', cmap(id,:))
                hold on
                title('target- NP sub')
                ylim([-0.03 0.1])
                subplot(3,3,7)
                plot(tt, mean(dS(id).releaseAlignDFoverF(:,iR,1,dS(id).successIx),4), 'Color', cmap(id,:))
                hold on
                title('release- ROI')
                ylim([-0.03 0.1])
                subplot(3,3,8)
                plot(tt, mean(dS(id).releaseAlignDFoverF(:,iR,2,dS(id).successIx),4), 'Color', cmap(id,:))
                hold on
                title('release- surround')
                ylim([-0.03 0.1])
                subplot(3,3,9)
                plot(tt, mean(dS(id).releaseAlignDFoverF_npsub(:,iR,dS(id).successIx),3), 'Color', cmap(id,:))
                hold on
                title('release- NP sub')
                ylim([-0.03 0.1])
            end
            suptitle(['Area ' area_list(iR,:)])
            print(fullfile(rc.structOutput, [mouse '_allAreaTCs_allAlignSIx_byNP.pdf']), '-dpdf')
        end
        figure; 
        for iR = 1:nROI
            subplot(3,3,iR)
            plot(squeeze(roi(iR).F)')
            title(area_list(iR,:))
        end
        legend('roi','np')
        title('Change in F over days')
        print(fullfile(rc.structOutput, [mouse '_FoverDays.pdf']), '-dpdf')
    end

end

    

