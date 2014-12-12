mouse = '516';
date = '141003';
ImgFolder = '002+003+004';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructDFoverF.mat');
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

%% organize dataset
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


%% Auditory vs visual

%
for icyc = dataStructDFoverF.Cycles
    thisInd = find(dataStructDFoverF.cycBlock2ON{icyc} == 1);
    cycA_ind{1,icyc} = thisInd;
    clear thisInd
    thisInd = find(dataStructDFoverF.cycBlock2ON{icyc} == 0);
    cycV_ind{1,icyc} = thisInd;
    clear thisInd    
end

for icyc = dataStructDFoverF.Cycles
    thisData = dataTrialsCells{icyc};
    start = 1;
    for itrial = cycA_ind{icyc}
        thisDataAvgCells(:,start) = mean(thisData(:,itrial,:),3);
        start = start+1;
    end
    thisDataAvgTrials = mean(thisDataAvgCells,2);
    cycAvgTC_A{1,icyc} = thisDataAvgTrials;
    clear thisData thisDataAvgCells thisDataAvgTrials
end

for icyc = dataStructDFoverF.Cycles
    thisData = dataTrialsCells{icyc};
    start = 1;
    for itrial = cycV_ind{icyc}
        thisDataAvgCells(:,start) = mean(thisData(:,itrial,:),3);
        start = start+1;
    end
    thisDataAvgTrials = mean(thisDataAvgCells,2);
    cycAvgTC_V{1,icyc} = thisDataAvgTrials;
    clear thisData thisDataAvgCells thisDataAvgTrials
end


nCycles = max(dataStructDFoverF.Cycles);
figure;
colorA = colormap(copper(nCycles));
colorV = colormap(hsv(nCycles));
for icyc = dataStructDFoverF.Cycles
    cycColor = colorA(icyc,:);
    plot(cycAvgTC_A{icyc},'Color',cycColor);
    hold on
    cycColor = colorV(icyc,:);
    plot(cycAvgTC_V{icyc},'Color',cycColor);
    hold on
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    axis([0 120 -0.05 0.1]);
    hold on
end

