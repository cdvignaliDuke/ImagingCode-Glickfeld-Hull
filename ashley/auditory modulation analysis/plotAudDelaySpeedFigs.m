%% speed
speedDataAllExpt = getSpeedDataAcrossExpt(ms,fn);

[alignVisStimEaExptFig, alignVisStimFig] = ...
    plotSpeedAlignVisStim(speedDataAllExpt,'speedAlignVisStim',...
    delayColors,exptNames,exptColors);
[alignStartEaExptFig, alignStartFig] = ...
    plotSpeedAlignVisStim(speedDataAllExpt(1:nDelay-2),'speedAlignStart',...
    delayColors,exptNames,exptColors);

speedRespFig = plotSpeedResp2VisStimEaExpt(speedDataAllExpt,exptColors);

%% functions for speed plotting
function [alignSpeedEaExptFig, alignSpeedFig] = plotSpeedAlignVisStim(...
    speedData,speedField2Plot,delayColors,exptNames,exptColors)
    ndelay = size(speedData,2);
    nexpt = size(exptColors,1);

    eaExptLegend = [];
    allExptLegend = [];
    alignSpeedEaExptFig = figure;
    alignSpeedFig = figure;
    [spRowsDelays,spColDelays] = optimizeSubplotDim(ndelay);
    for idelay = 1:ndelay
        timestamp = eval(['speedData(idelay).' speedField2Plot 'TimeStampS']);
        figure(alignSpeedEaExptFig)
        subplot(spRowsDelays,spColDelays,idelay)    
        for iexp = 1:nexpt
            y = eval(['speedData(idelay).' speedField2Plot 'EaExpt(:,iexp)']);
            yerr = eval(['speedData(idelay).' speedField2Plot 'ErrEaExpt(:,iexp)']);
            h = shadedErrorBar(timestamp,y,yerr,[],1);
            h.mainLine.Color = exptColors(iexp,:);
            hold on
            v = vline(0,'k--');
            if idelay == 1
                if ~isnan(mean(y))
                    eaExptLegend(iexp) = h.mainLine;
                else
                    eaExptLegend(iexp) = v;
                end 
            end
        end
        title(sprintf('%sms delay',speedData(idelay).trType));
        figXAxis([],'time(s)',[timestamp(1) timestamp(end)])
        figYAxis([],'speed diff from baseline(m/s)',[])
        ax = gca;
        ax.TickDir = 'out';
        ax.Box = 'off';
        
        y = nanmean(eval(['speedData(idelay).' speedField2Plot 'EaExpt']),2);
        yerr = ste(eval(['speedData(idelay).' speedField2Plot 'EaExpt']),2);
        figure(alignSpeedFig)
        hold on
        h = shadedErrorBar(timestamp,y,yerr,[],1);
        h.mainLine.Color = delayColors(idelay,:);
        h.mainLine.LineWidth = 2;
        allExptLegend(idelay) = h.mainLine;
        hold on
        vline(0,'k--')
    end
    figure(alignSpeedEaExptFig)
    legend(eaExptLegend,exptNames)
    figure(alignSpeedFig)
    legend(allExptLegend,{speedData.trType})
    figXAxis([],'time(s)',[timestamp(1) timestamp(end)])
    figYAxis([],'speed diff from baseline(m/s)',[])
    ax = gca;
    ax.TickDir = 'out';
    ax.Box = 'off';
end

function speedRespFig = plotSpeedResp2VisStimEaExpt(speedData,exptColors)
    ndelay = size(speedData,2);
    speedRespFig = figure;
    [spRowsDelays,spColDelays] = optimizeSubplotDim(ndelay);
    visOnlyInd = strcmp({speedData.trType},'vis only');
    visOnlySpeedResp = speedData(visOnlyInd).speedResp;
    visOnlySpeedErr = speedData(visOnlyInd).speedRespErr;
    suptitle('avg speed during vis stim response window')
    for idelay = 1:ndelay-1
        audDelaySpeedResp = speedData(idelay).speedResp;
        audDelaySpeedErr = speedData(idelay).speedRespErr;
        respLim = minMaxDFFRespAxLim([visOnlySpeedResp,audDelaySpeedResp],0.01);
        subplot(spRowsDelays,spColDelays,idelay)
    %     h = scatter(visOnlySpeedResp,audDelaySpeedResp,'ko');
        h = errorbarxy(visOnlySpeedResp,audDelaySpeedResp,visOnlySpeedErr,...
            audDelaySpeedErr,{'ko','k','k'});
        h.hMain.MarkerFaceColor = 'k';
        h.hMain.MarkerEdgeColor = [1,1,1];
        hold on
        h = scatter(visOnlySpeedResp,audDelaySpeedResp,[],exptColors,'filled');
        h.MarkerEdgeColor = [1,1,1];
        hold on
        plot(respLim,respLim)
        figXAxis([],'vis only',respLim)
        figYAxis([],speedData(idelay).trType,respLim)
        figAxForm([])
    end
end

function neuronCorr2SpeedFig = plotMeanNeuron2SpeedCorr(speedData)
    ndelay = size(speedData,2);
    audDelaySpeedCorr = zeros(1,ndelay);
    audDelaySpeedCorrErr = zeros(1,ndelay);
    for idelay = 1:ndelay
        audDelaySpeedCorr(idelay) = nanmean(speedData(idelay).neuronCorr2Speed);
        audDelaySpeedCorrErr(idelay) = ste(speedData(idelay).neuronCorr2Speed,2);
    end
    neuronCorr2SpeedFig = figure;
    h = errorbar(1:ndelay,audDelaySpeedCorr,audDelaySpeedCorrErr,'ko');
    h.MarkerFaceColor = 'k';
    h.MarkerEdgeColor = [1,1,1];
    figXAxis([],'aud delay',[0 nDelay+1],1:ndelay,{speedData.trType})
    figYAxis([],'mean corr (dff:speed)',[-0.2 0.2])
    figAxForm([])
end

