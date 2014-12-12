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



    %response to target by cell by trial
for icyc = 1:siz
    thisRL = 40;
    thisTarget = dataStructVar.cycTargetOn{icyc};
    thisData = dataTrialsCells{icyc};
    for itrial = dataStructVar.cycNTrials{icyc}
        for icell = 1:nCells
        if isnan(thisTarget(itrial))
            thisTarRsp(:,itrial,icell) = NaN(thisRL,1,1);
        else
            thisInd = thisTarget-19:thisTarget+20;
            thisTarRsp(:,itrial,icell) = thisData(thisInd,itrial,icell);
        end
        end
    end
    dataTrialsCells_TarRsp{1,icyc} = thisTarRsp;
    clear thisRL thisTarget thisData thisInd thisTarRsp
end

    % target presence index
for icyc = 1:siz
    thisTarget = dataStructVar.cycTargetOn{icyc};
    thisIndNaN = find(isnan(thisTarget));
    thisIndTar = find(isnan(thisTarget)==0);
    cycNaN_ind{1,icyc} = thisIndNaN;
    cycTar_ind{1,icyc} = thisIndTar;
    clear thisTarget thisIndNaN thisIndTar
end
    





% Auditory vs visual
    %auditory and visual indexes
for icyc = 1:siz
    thisInd = find(dataStructVar.cycBlock2ON{icyc} == 1);
    cycA_ind{1,icyc} = thisInd;
    clear thisInd
    thisInd = find(dataStructVar.cycBlock2ON{icyc} == 0);
    cycV_ind{1,icyc} = thisInd;
    clear thisInd    
end  

% for icyc = 1:siz
%     thisData = dataTrialsCells_TarRsp{icyc};
%     start = 1;
%     thisIndTar = cycTar_ind{icyc};
%     thisIndA = cycA_ind{icyc};
%     thisInd = intersect(thisIndTar,thisIndA);
%     thisDataAvgCells = [];
%     for itrial = thisInd
%         thisDataAvgCells(:,start) = mean(thisData(:,itrial,:),3);
%         start = start+1;
%     end
%     if isempty(thisDataAvgCells)
%         cycAvgTC_A{1,icyc} = [];
%     else        
%         thisDataAvgTrials = mean(thisDataAvgCells,2);
%         cycAvgTC_A{1,icyc} = thisDataAvgTrials;
%     end
%     clear thisData thisDataAvgCells thisDataAvgTrials thisIndTar thisIndA thisInd
% end
% 
% for icyc = 1:siz
%     thisData = dataTrialsCells_TarRsp{icyc};
%     start = 1;
%     thisIndTar = cycTar_ind{icyc};
%     thisIndA = cycV_ind{icyc};
%     thisInd = intersect(thisIndTar,thisIndA);
%     thisDataAvgCells = [];
%     for itrial = thisInd
%         thisDataAvgCells(:,start) = mean(thisData(:,itrial,:),3);
%         start = start+1;
%     end
%     if isempty(thisDataAvgCells)
%         cycAvgTC_V{1,icyc} = [];
%     else        
%         thisDataAvgTrials = mean(thisDataAvgCells,2);
%         cycAvgTC_V{1,icyc} = thisDataAvgTrials;
%     end
%     clear thisData thisDataAvgCells thisDataAvgTrials thisIndTar thisIndA thisInd
% end

for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    start = 1;
    thisIndTar = cycTar_ind{icyc};
    thisIndA = cycA_ind{icyc};
    thisInd = intersect(thisIndTar,thisIndA);
    thisDataAvgCells = [];
    for itrial = thisInd
        thisDataAvgCells(:,start) = mean(thisData(:,itrial,:),3);
        start = start+1;
    end
    if isempty(thisDataAvgCells)
        cycAvgTC_A{1,icyc} = [];
    else        
        thisDataAvgTrials = mean(thisDataAvgCells,2);
        cycAvgTC_A{1,icyc} = thisDataAvgTrials;
    end
    clear thisData thisDataAvgCells thisDataAvgTrials thisIndTar thisIndA thisInd
end

for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    start = 1;
    thisIndTar = cycTar_ind{icyc};
    thisIndA = cycV_ind{icyc};
    thisInd = intersect(thisIndTar,thisIndA);
    thisDataAvgCells = [];
    for itrial = thisInd
        thisDataAvgCells(:,start) = mean(thisData(:,itrial,:),3);
        start = start+1;
    end
    if isempty(thisDataAvgCells)
        cycAvgTC_V{1,icyc} = [];
    else        
        thisDataAvgTrials = mean(thisDataAvgCells,2);
        cycAvgTC_V{1,icyc} = thisDataAvgTrials;
    end
    clear thisData thisDataAvgCells thisDataAvgTrials thisIndTar thisIndA thisInd
end

nCycles = max(dataStructVar.Cycles);
figure;
for icyc = 1:siz
    plot(cycAvgTC_A{icyc},'r');
    hold on
    plot(cycAvgTC_V{icyc},'g');
    hold on
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    axis([0 125 -0.05 0.2]);
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

    %average response to target by cell by trial
for icyc = 1:siz
    thisRL = 40;
    thisTarget = dataStructVar.cycTargetOn{icyc};
    thisData = dataTrialsCells{icyc};
    for itrial = dataStructVar.cycNTrials{icyc}
        for icell = 1:nCells
        if isnan(thisTarget(itrial))
            thisTarRsp(:,itrial,icell) = NaN(thisRL,1,1);
        else
            thisInd = thisTarget-19:thisTarget+20;
            thisTarRsp(:,itrial,icell) = thisData(thisInd,itrial,icell);
        end
        end
    end
    dataTrialsCells_TarRsp{1,icyc} = thisTarRsp;
    clear thisRL thisTarget thisData thisInd thisTarRsp
end

% Auditory vs visual
    %auditory and visual indexes
for icyc = 1:siz
    thisInd = find(dataStructVar.cycBlock2ON{icyc} == 1);
    cycA_ind{1,icyc} = thisInd;
    clear thisInd
    thisInd = find(dataStructVar.cycBlock2ON{icyc} == 0);
    cycV_ind{1,icyc} = thisInd;
    clear thisInd    
end

    %avg all trials all cells
for icyc = 1:siz
    thisData = dataTrialsCells_TarRsp{icyc};
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
    thisData = dataTrialsCells_TarRsp{icyc};
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
    plot(cycAvgTC_A{icyc},'r');
    hold on
    plot(cycAvgTC_V{icyc},'g');
    hold on
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    axis([0 125 -0.05 0.05]);
    hold on
end

    %avg all auditory and visual trials
for icyc = 1:siz
        thisData = dataTrialsCells_TarRsp{icyc};
        thisInd = cycA_ind{icyc};
    for icell = 1:nCells
        thisAvg(:,icell) = squeeze(mean(thisData(:,thisInd,icell),2));
        dataCells_A{icyc} = thisAvg;
        clear thisAvg
    end
    clear thisData thisInd
end

for icyc = 1:siz
        thisData = dataTrialsCells_TarRsp{icyc};
        thisInd = cycV_ind{icyc};
    for icell = 1:nCells
        thisAvg(:,icell) = squeeze(mean(thisData(:,thisInd,icell),2));
        dataCells_V{icyc} = thisAvg;
        clear thisAvg
    end
    clear thisData thisInd
end

