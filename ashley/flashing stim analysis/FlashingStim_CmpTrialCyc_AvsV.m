mouse = '516';
SubNum = '516';
date = '141003';
ImgFolder = '002+003+004';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructVar.mat');
load('dataTC.mat');
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

%% dirTuning dataset
%get average response to vis stim for each cycle
siz = size(dataStructVar.Cycles,2);
data = dataTC.dirTuning;
nCells = size(data{1},2);

for icyc = 1:siz
    for icell = 1:nCells
        thisL = dataStructVar.cycTrialL{icyc};
        thisData = data{icyc};
        for itrial = 1:dataStructVar.cycNTrials{icyc}
            thisMat(:,itrial,icell) = thisData(1+(thisL.*(itrial-1)):thisL.*itrial,icell);
        end
        dataTrialsCells{icyc} = thisMat;
    end
    clear thisL thisData thisMat
end

    %avg all trials all cells
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    for itrial = 1:dataStructVar.cycNTrials{icyc}
        thisDataAvgCells(:,itrial) = mean(thisData(:,itrial,:),3);
    end
    thisDataAvgTrials = mean(thisDataAvgCells,2);
    cycAvgTC{1,icyc} = thisDataAvgTrials;
    clear thisData thisDataAvgCells thisDataAvgTrials
end

% Auditory vs visual

for icyc = 1:siz
    thisInd = find(dataStructVar.cycBlock2ON{icyc} == 1);
    cycA_ind{1,icyc} = thisInd;
    clear thisInd
    thisInd = find(dataStructVar.cycBlock2ON{icyc} == 0);
    cycV_ind{1,icyc} = thisInd;
    clear thisInd    
end

for icyc = 1:siz
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

for icyc = 1:siz
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

for icyc = 1:siz
        thisData = dataTrialsCells{icyc};   
    for icell = 1:nCells
        thisAvg(:,icell) = squeeze(mean(thisData(:,:,icell),2));
        dataCells{icyc} = thisAvg;
        clear thisAvg
    end
    clear thisData
end

nCycles = max(dataStructVar.Cycles);
figure;
colorA = colormap(copper(nCycles));
colorV = colormap(hsv(nCycles));
for icyc = 1:siz
    cycColor = colorA(icyc,:);
    plot(cycAvgTC_A{icyc},'r');
    hold on
    cycColor = colorV(icyc,:);
    plot(cycAvgTC_V{icyc},'g');
    hold on
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    axis([0 125 -0.05 0.05]);
    hold on
end

figure;
for icyc = 1:siz
    subplot(2,4,icyc);
    plot(cycAvgTC_A{icyc},'r');
    hold on
    plot(cycAvgTC_V{icyc},'g');
    hold on
    axis([0 125 -0.05 0.05]);
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    hold on
end

%% FSfirstRsp dataset
data = dataTC.FSfirstRsp;
nCells = size(data{1},2);

for icyc = 1:siz
    for icell = 1:nCells
        thisL = dataStructVar.cycTrialL{icyc};
        thisData = data{icyc};
        for itrial = 1:dataStructVar.cycNTrials{icyc}
            thisMat(:,itrial,icell) = thisData(1+(thisL.*(itrial-1)):thisL.*itrial,icell);
        end
        dataTrialsCells{icyc} = thisMat;
    end
    clear thisL thisData thisMat
end

    %avg all trials all cells
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    for itrial = 1:dataStructVar.cycNTrials{icyc}
        thisDataAvgCells(:,itrial) = mean(thisData(:,itrial,:),3);
    end
    thisDataAvgTrials = mean(thisDataAvgCells,2);
    cycAvgTC{1,icyc} = thisDataAvgTrials;
    clear thisData thisDataAvgCells thisDataAvgTrials
end

% Auditory vs visual

for icyc = 1:siz
    thisInd = find(dataStructVar.cycBlock2ON{icyc} == 1);
    cycA_ind{1,icyc} = thisInd;
    clear thisInd
    thisInd = find(dataStructVar.cycBlock2ON{icyc} == 0);
    cycV_ind{1,icyc} = thisInd;
    clear thisInd    
end

for icyc = 1:siz
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

for icyc = 1:siz
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


nCycles = max(dataStructVar.Cycles);
figure;
for icyc = 1:siz
    cycColor = colorA(icyc,:);
    plot(cycAvgTC_A{icyc},'r');
    hold on
    cycColor = colorV(icyc,:);
    plot(cycAvgTC_V{icyc},'g');
    hold on
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    axis([0 125 -0.05 0.05]);
    hold on
end

figure;
for icyc = 1:siz
    subplot(2,4,icyc);
    plot(cycAvgTC_A{icyc},'r');
    hold on
    plot(cycAvgTC_V{icyc},'g');
    hold on
    axis([0 125 -0.05 0.05]);
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    hold on
end
