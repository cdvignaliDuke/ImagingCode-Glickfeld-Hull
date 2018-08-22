function plotEachCellTC(exptName,tc_mean,tc_err,nBaselineFr,cycLengthFr,frRateHz)
    [trLength,nc] = size(tc_mean);
    
    if nc > 64
        nRows = 8;
        nCols = 8;
        nFig = ceil(nc./64);
    else
        [nRows,nCols] = optimizeSubplotDim(nc);
        nFig = 1;
    end
    ttMs = (((1:trLength) - nBaselineFr - 2)./frRateHz).*1000;
    stimOnTimeMs = (((0:1).*cycLengthFr)./frRateHz).*1000;
    stimOffTimeMs = (((1:2).*cycLengthFr)./frRateHz).*1000-250;

    set(0,'defaultAxesFontSize',6)
    for ifig = 1:nFig
        if ifig == 1
            start = 1;
            if nc > 64
                cellEnd = 64;
                maxPlot = 64;
            else
                cellEnd = nc;
                maxPlot = nc;
            end
        else
            start = start+64;
            if ifig == nFig
                cellEnd = nc;
                maxPlot = nc - ((ifig-1)*64);
            else
                cellEnd = cellEnd+64;
                maxPlot = 64;
            end
        end
        cellInd = start:cellEnd;
        figure
        suptitle(exptName)
        for i = 1:maxPlot
            subplot(nRows,nCols,i)
            h = shadedErrorBar(ttMs, tc_mean(:,cellInd(i)), ...
                tc_err(:,cellInd(i)), 'k');
            hold on
            vline(0,'k:')
            vline(stimOnTimeMs,'g-')
            vline(stimOffTimeMs, 'b-')
            figXAxis([],'Time (ms)',[])
            figYAxis([],'dF/F',[])
            figAxForm([],0)
            title(cellInd(i))
        end
    end
end