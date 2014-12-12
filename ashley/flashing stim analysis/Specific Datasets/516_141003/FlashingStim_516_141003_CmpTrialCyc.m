mouse = '516';
date = '141003';
ImgFolder = '002+003+004';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructDFoverF.mat')

%get average response to vis stim for each cycle
nCells = dataStructDFoverF.nCells;

for icyc = dataStructDFoverF.Cycles
    for icell = 1:nCells
        thisL = dataStructDFoverF.cycTrialL{icyc};
        thisData = dataStructDFoverF.cycData{icyc};
        thisTC = dataStructDFoverF.cellTC{icyc};
        for itrial = 1:dataStructDFoverF.cycNTrials{icyc}
            thisMat(:,itrial,icell) = thisTC(1+(thisL.*(itrial-1)):thisL.*itrial,icell);
        end
        dataTrialsCells{icyc} = thisMat;
    end
    clear thisL thisData thisTC thisMat
end

    %avg all trials all cells
for icyc = dataStructDFoverF.Cycles
    thisData = dataTrialsCells{icyc};
    for itrial = 1:dataStructDFoverF.cycNTrials{icyc}
        thisDataAvgCells(:,itrial) = mean(thisData(:,itrial,:),3);
    end
    thisDataAvgTrials = mean(thisDataAvgCells,2);
    cycAvgTC{1,icyc} = thisDataAvgTrials;
    clear thisData thisDataAvgCells thisDataAvgTrials
end

    %FIG for traces
nCycles = max(dataStructDFoverF.Cycles);
figure;
color = colormap(jet(nCycles));
for icyc = dataStructDFoverF.Cycles
    cycColor = color(icyc,:);
    plot(cycAvgTC{icyc},'Color',cycColor);
    hold on
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    hold on
end


