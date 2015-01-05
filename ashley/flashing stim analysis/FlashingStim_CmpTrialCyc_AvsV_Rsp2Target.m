mouse = 'AW07';
SubNum = '607';
date = '141215';
ImgFolder = '004';
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
    thisIndV = cycV_ind{icyc};
    thisInd = intersect(thisIndTar,thisIndV);
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
    v = 5+ ((i-1)*(dataStructVar.ONfr+dataStructVar.OFFfr));
    vline(v,'k');
    end
    axis([0 125 -0.05 0.1]);
    hold on
end

figure;
for icyc = 1:siz
    subplot(2,4,icyc);
    plot(cycAvgTC_A{icyc},'r');
    hold on
    plot(cycAvgTC_V{icyc},'g');
    hold on
    axis([0 125 -0.02 0.06]);
    for i = 1:nCycles+1
    v = 5+ (i-1)*(dataStructVar.ONfr+dataStructVar.OFFfr);
    vline(v,'k');
    end
    t = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
    vline(t,'c');
    hold on
end

%% frames around target only

%     %avg response for each cycle type
% for icyc = 1:siz
%     thisData = cycAvgTC_A{icyc};
%     thisTar = ceil(((dataStructVar.ONms+dataStructVar.OFFms).*dataStructVar.RateFRperMS).*(icyc+1));
%     thisInd = thisTar-10:thisTar+19;
%     if isempty(thisData)
%         thisTrace = [];
%     else
%         thisTrace = thisData(thisInd,:);
%     end
%     cycAvgTC_ATar{1,icyc} = thisTrace;
%     clear thisData thisTar thisInd thisTrace
% end
% 
% for icyc = 1:siz
%     thisData = cycAvgTC_V{icyc};
%     thisTar = ceil(((dataStructVar.ONms+dataStructVar.OFFms).*dataStructVar.RateFRperMS).*(icyc+1));
%     thisInd = thisTar-10:thisTar+19;
%     if isempty(thisData)
%         thisTrace = [];
%     else
%         thisTrace = thisData(thisInd,:);
%     end
%     cycAvgTC_VTar{1,icyc} = thisTrace;
%     clear thisData thisTar thisInd thisTrace
% end
% 
%     %avg response for all targets (not weighted)
% 
% for icyc = 1:siz
%     thisData = dataTrialsCells{icyc};
%     thisIndA = cycA_ind{icyc};
%     thisIndTar = cycTar_ind{icyc};
%     thisInd = intersect(thisIndTar,thisIndA);
%     thisTar = ceil(((dataStructVar.ONms+dataStructVar.OFFms).*dataStructVar.RateFRperMS).*(icyc+1));
%     thisTarL = thisTar-10:thisTar+19;
%         if isempty(thisInd)
%             thisDataTrials_A = [];
%         else
%         thisDataTrials_A(:,thisInd) = mean(thisData(thisTarL,thisInd,:),3);
%         end
%     dataTrials_A{icyc} = thisDataTrials_A;
%     clear thisData thisInd thisTar thisTarInd thisDataTrials_A
% end
% 
% meanDataTrials_A = mean((cell2mat(dataTrials_A)),2);
% 
% for icyc = 1:siz
%     thisData = dataTrialsCells{icyc};
%     thisIndV = cycV_ind{icyc};
%     thisIndTar = cycTar_ind{icyc};
%     thisInd = intersect(thisIndTar,thisIndV);
%     thisTar = ceil(((dataStructVar.ONms+dataStructVar.OFFms).*dataStructVar.RateFRperMS).*(icyc+1));
%     thisTarL = thisTar-10:thisTar+19;
%         if isempty(thisInd)
%             thisDataTrials_V = [];
%         else
%         thisDataTrials_V(:,thisInd) = mean(thisData(thisTarL,thisInd,:),3);
%         end
%     dataTrials_V{icyc} = thisDataTrials_V;
%     clear thisData thisInd thisTar thisTarInd thisDataTrials_A
% end
% 
% meanDataTrials_V = mean((cell2mat(dataTrials_V)),2);

% for datasets that used FlashingStim_2P_Frames
%     %avg response for each cycle type
for icyc = 1:siz
    thisData = cycAvgTC_A{icyc};
    thisTar = (dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
    thisInd = thisTar-10:thisTar+19;
    if isempty(thisData)
        thisTrace = [];
    else
        thisTrace = thisData(thisInd,:);
    end
    cycAvgTC_ATar{1,icyc} = thisTrace;
    clear thisData thisTar thisInd thisTrace
end

for icyc = 1:siz
    thisData = cycAvgTC_V{icyc};
    thisTar = (dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
    thisInd = thisTar-10:thisTar+19;
    if isempty(thisData)
        thisTrace = [];
    else
        thisTrace = thisData(thisInd,:);
    end
    cycAvgTC_VTar{1,icyc} = thisTrace;
    clear thisData thisTar thisInd thisTrace
end

%     %avg response for all targets (not weighted)
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    thisIndA = cycA_ind{icyc};
    thisIndTar = cycTar_ind{icyc};
    thisInd = intersect(thisIndTar,thisIndA);
    thisTar = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
    thisTarL = thisTar-10:thisTar+19;
        if isempty(thisInd)
            thisDataTrials_A = [];
        else
        thisDataTrials_A(:,thisInd) = mean(thisData(thisTarL,thisInd,:),3);
        end
    dataTrials_A{icyc} = thisDataTrials_A;
    clear thisData thisInd thisTar thisTarInd thisDataTrials_A
end

meanDataTrials_A = mean((cell2mat(dataTrials_A)),2);

for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    thisIndV = cycV_ind{icyc};
    thisIndTar = cycTar_ind{icyc};
    thisInd = intersect(thisIndTar,thisIndV);
    thisTar = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
    thisTarL = thisTar-10:thisTar+19;
        if isempty(thisInd)
            thisDataTrials_V = [];
        else
        thisDataTrials_V(:,thisInd) = mean(thisData(thisTarL,thisInd,:),3);
        end
    dataTrials_V{icyc} = thisDataTrials_V;
    clear thisData thisInd thisTar thisTarInd thisDataTrials_A
end

meanDataTrials_V = mean((cell2mat(dataTrials_V)),2);



    %plots
    
figure;
for icyc = 1:siz
    plot(cycAvgTC_ATar{icyc},'r');
    hold on
    plot(cycAvgTC_VTar{icyc},'g');
    hold on
    vline(10,'k');
    axis([0 35 -0.05 0.2]);
    hold on
end

figure;
for icyc = 1:siz
    plot(meanDataTrials_A,'r');
    hold on
    plot(meanDataTrials_V,'g');
    hold on
    vline(10,'k');
    axis([0 35 -0.05 0.05]);
    hold on
end