% auditory vs. visual trace for average of all trials all cells up to, but
% not including the target, this way there are many trials at least for low
% numbers of cycles - this is an attempt to smooth the curves and reduce
% noise


mouse = 'AW07';
SubNum = '607';
date = '150121';
ImgFolder = '003';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructVar.mat');
load('dataTC.mat');
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

%% dirTuning dataset
%get average response to vis stim for each cycle
siz = size(dataStructVar.Cycles,2);
data = dataTC.retTuning;
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
    if isempty(cycV_ind{icyc})
        thisDataAvgCells = [];
    else
        for itrial = cycV_ind{icyc}
            thisDataAvgCells(:,start) = mean(thisData(:,itrial,:),3);
            start = start+1;
        end
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

% average all cells all trials, no target, not sorted by cyc type for that
% trial, but cumulatively
start = 1;
for ilength = dataStructVar.Cycles
    L(start) = ceil(5+ (ilength-1)*11);
    start = start+1;
end

for ilength = 1:siz
    L_ind{ilength} = ilength:siz;
end

for icyc = 1:siz
    thisL = L(icyc);
    x = 1;
    for i = L_ind{icyc}
        thisInd = cycA_ind{i};
        y = x + (size(thisInd,2))-1;
        thisData = dataTrialsCells{i};
        thisData_A(:,x:y,:) = thisData(1:thisL,thisInd,:);
        x = y+1;
        clear thisData
    end
    cycNoTarget_A{1,icyc} = thisData_A;
    thisData_meanA = mean(thisData_A,3);
    cycNoTargetAvgTC_A{1,icyc} = mean(thisData_meanA,2);
    clear thisData_A thisData_meanA
end

for icyc = 1:siz
    thisL = L(icyc);
    x = 1;
    for i = L_ind{icyc}
        thisInd = cycV_ind{i};
        y = x + (size(thisInd,2))-1;
        thisData = dataTrialsCells{i};
        thisData_V(:,x:y,:) = thisData(1:thisL,thisInd,:);
        x = y+1;
        clear thisData
    end
    cycNoTarget_V{1,icyc} = thisData_V;
    thisData_meanV = mean(thisData_V,3);
    cycNoTargetAvgTC_V{1,icyc} = mean(thisData_meanV,2);
    clear thisData_V thisData_meanV
end


% nCycles = max(dataStructVar.Cycles);
% figure;
% for icyc = 1:siz
%     plot(cycNoTargetAvgTC_A{icyc},'r');
%     hold on
%     plot(cycNoTargetAvgTC_V{icyc},'g');
%     hold on
%     for i = 1:nCycles+1
%     v = 5+ (i-1)*10.5;
%     vline(v,'k');
%     end
%     axis([0 125 -0.05 0.08]);
%     hold on
% end

figure;
for icyc = 1:siz
    subplot(2,4,icyc);
    plot(cycNoTargetAvgTC_A{icyc},'r');
    hold on
    plot(cycNoTargetAvgTC_V{icyc},'g');
    hold on
    C = dataStructVar.Cycles(icyc);
    for i = 1:C
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    axis([0 v -0.05 0.05]);
    hold on
end
    