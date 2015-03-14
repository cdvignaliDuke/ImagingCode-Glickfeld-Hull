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


%     %response to target by cell by trial
% for icyc = 1:siz
%     thisRL = 60;
%     thisTarget = dataStructVar.cycTargetOn{icyc}-dataStructVar.cycLeverDown{icyc};
%     thisData = dataTrialsCells{icyc};
%     for itrial = dataStructVar.cycNTrials{icyc}
%         for icell = 1:nCells
%         if isnan(thisTarget(itrial))
%             thisTarRsp(:,itrial,icell) = NaN(thisRL,1,1);
%         else
%             thisInd = thisTarget-29:thisTarget+30;
%             thisTarRsp(:,itrial,icell) = thisData(thisInd,itrial,icell);
%         end
%         end
%     end
%     dataTrialsCells_TarRsp{1,icyc} = thisTarRsp;
%     clear thisRL thisTarget thisData thisInd thisTarRsp
% end

    % target presence index
for icyc = 1:siz
    thisTarget = dataStructVar.cycTargetOn{icyc};
    thisIndNaN = find(isnan(thisTarget));
    thisIndTar = find(isnan(thisTarget)==0);
    cycNaN_ind{1,icyc} = thisIndNaN;
    cycTar_ind{1,icyc} = thisIndTar;
    clear thisTarget thisIndNaN thisIndTar
end
    
    % frames around target (only trials with actual target present) - cycles unsorted
% start = 1;
% for icyc = 1:siz
%     thisInd = cycTar_ind{icyc};
%     if isempty(thisInd)
%         x = 1;
%     else
%     thisData = dataTrialsCells_TarRsp{icyc};
%     thisDataTar = thisData(:,thisInd,:);
%     tTrials = size(thisDataTar,2);
%     targetStimResp_data(:,start:start+(tTrials-1),1:nCells) = thisDataTar;
%     start = start+tTrials;
%     end
%     clear thisData tTrials thisInd thisDataTar
% end
% 
% figure;
% for icyc = 1:siz
%     plot(mean(mean(dataTrialsCells_TarRsp{icyc},2),3),'r');
%     hold on
%     vline(30,'k');
%     axis([0 60 -0.02 0.02]);
%     hold on
% end
%     



%% Auditory vs visual
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
    v = 5+ (i-1)*(dataStructVar.ONfr+dataStructVar.OFFfr)
    vline(v,'k');
    end
    t = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
    vline(t,'c');
    hold on
end

% nCycles = max(dataStructVar.Cycles);
% figure;
% for icyc = 1:siz
%     plot(cycAvgTC_A{icyc},'r');
%     hold on
%     plot(cycAvgTC_V{icyc},'g');
%     hold on
%     v = 5+ ((icyc-1)*((dataStructVar.ONms+dataStructVar.OFFms)*dataStructVar.RateFRperMS));
%     vline(v,'k');
%     hold on
%     axis([0 125 -0.05 0.2]);
%     hold on
% end
% 
% figure;
% for icyc = 1:siz
%     subplot(2,4,icyc);
%     plot(cycAvgTC_A{icyc},'r');
%     hold on
%     plot(cycAvgTC_V{icyc},'g');
%     hold on
%     axis([0 125 -0.02 0.2]);
%     for i = 1:nCycles
%     v = 5+ (i-1)*((dataStructVar.ONms+dataStructVar.OFFms).*dataStructVar.RateFRperMS);
%     vline(v,'k');
%     end
%     t = 5+((dataStructVar.ONms+dataStructVar.OFFms).*dataStructVar.RateFRperMS).*(icyc);
%     vline(t,'c');
%     hold on
% end

%% frames around target only

    %avg response for each cycle type
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

    %avg response for all targets (not weighted)

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

start = 1;
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    thisIndA = cycA_ind{icyc};
    thisIndTar = cycTar_ind{icyc};
    thisInd = intersect(thisIndTar,thisIndA);
    thisTarMat = dataStructVar.cycTargetOn{icyc}-dataStructVar.cycLeverDown{icyc};
    if isempty(thisInd)
        x = []
    else
        thisTrialN = size(thisInd,2);
        thisTarMat = thisTarMat(thisInd);
        thisTarL_start = thisTarMat-10';
        thisTarL_finish = thisTarMat+19';
        thisDataA = thisData(:,thisInd,:);
        for i = 1:thisTrialN
            dataTrialsCells_A(:,start,:) = thisDataA(thisTarL_start(i):thisTarL_finish(i),i,:);
        start = start+1;
        end
    clear thisData thisIndA thisIndTar thisInd thisTarMat x thisTrialN thisTarL_start thisTarL_finish clear thisDataA
    end
end

meanDataTrials_A = mean(mean(dataTrialsCells_A,3),2);
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

start = 1;
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    thisIndV = cycV_ind{icyc};
    thisIndTar = cycTar_ind{icyc};
    thisInd = intersect(thisIndTar,thisIndV);
    thisTarMat = dataStructVar.cycTargetOn{icyc}-dataStructVar.cycLeverDown{icyc};
    if isempty(thisInd)
        x = []
    else
        thisTrialN = size(thisInd,2);
        thisTarMat = thisTarMat(thisInd);
        thisTarL_start = thisTarMat-10';
        thisTarL_finish = thisTarMat+19';
        thisDataA = thisData(:,thisInd,:);
        for i = 1:thisTrialN
            dataTrialsCells_V(:,start,:) = thisDataA(thisTarL_start(i):thisTarL_finish(i),i,:);
        start = start+1;
        end
    clear thisData thisIndA thisIndTar thisInd thisTarMat x thisTrialN thisTarL_start thisTarL_finish clear thisDataA
    end
end

meanDataTrials_V = mean(mean(dataTrialsCells_V,3),2);

% % for datasets that used FlashingStim_2P_Frames
% %     %avg response for each cycle type
% for icyc = 1:siz
%     thisData = cycAvgTC_A{icyc};
%     thisTar = (dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
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
%     thisTar = (dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
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
% %     %avg response for all targets (not weighted)
% for icyc = 1:siz
%     thisData = dataTrialsCells{icyc};
%     thisIndA = cycA_ind{icyc};
%     thisIndTar = cycTar_ind{icyc};
%     thisInd = intersect(thisIndTar,thisIndA);
%     thisTar = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
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
%     thisTar = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
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



    %plots
    
% figure;
% for icyc = 1:siz
%     plot(cycAvgTC_ATar{icyc},'r');
%     hold on
%     plot(cycAvgTC_VTar{icyc},'g');
%     hold on
%     vline(10,'k');
%     axis([0 35 -0.05 0.2]);
%     hold on
% end

figure;
for icyc = 1:siz
    plot(meanDataTrials_A,'r');
    hold on
    plot(meanDataTrials_V,'g');
    hold on
    vline(10,'k');
    axis([0 35 0 0.05]);
    hold on
end

%% find responsive cells (to target)

    %auditory
nTrials = size(dataTrialsCells_A,2);
for itrial = 1:nTrials
    for icell = 1:nCells
        preStimResp(itrial,icell) = mean(dataTrialsCells_A(1:10,itrial,icell),1);
    end
end

for itrial = 1:nTrials
    for icell = 1:nCells
        visStimResp(itrial,icell) = mean(dataTrialsCells_A(11:end,itrial,icell),1);
    end
end

visResp_ttest = ttest(preStimResp,visStimResp);
visResp_ind_A = find(visResp_ttest == 1);

clear nTrials preStimResp visStimResp visResp_ttest
    %visual
nTrials = size(dataTrialsCells_V,2);
for itrial = 1:nTrials
    for icell = 1:nCells
        preStimResp(itrial,icell) = mean(dataTrialsCells_V(1:10,itrial,icell),1);
    end
end

for itrial = 1:nTrials
    for icell = 1:nCells
        visStimResp(itrial,icell) = mean(dataTrialsCells_V(11:end,itrial,icell),1);
    end
end

visResp_ttest = ttest(preStimResp,visStimResp);
visResp_ind_V = find(visResp_ttest == 1);

visResp_ind_AandV = intersect(visResp_ind_A,visResp_ind_V);

%% Plot cells

%auditory responsive cells
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    dataTrialsCells_AResp{icyc} = thisData(:,:,visResp_ind_A);
end
%visually responsive cells
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    dataTrialsCells_VResp{icyc} = thisData(:,:,visResp_ind_V);
end

%both
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    dataTrialsCells_AandVResp{icyc} = thisData(:,:,visResp_ind_AandV);
end

%plot averages
figure;
for icyc = 1:siz
    subplot(2,4,icyc);
    plot(mean(mean(dataTrialsCells_AResp{icyc},3),2),'r');
    hold on
    plot(mean(mean(dataTrialsCells_VResp{icyc},3),2),'g');
    hold on
    axis([0 125 -0.02 0.06]);
    for i = 1:nCycles+1
    v = 5+ (i-1)*(dataStructVar.ONfr+dataStructVar.OFFfr)
    vline(v,'k');
    end
    t = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(icyc+1);
    vline(t,'c');
    hold on
end

%plot cells that respond to the auditory target
figure;
nACells = size(dataTrialsCells_AResp{1},3);

for icell = 1:nACells
    subplot(4,4,icell);
    thisIndA = cycA_ind{4};
    thisData = dataTrialsCells_AResp{4};
    plot(squeeze(mean(thisData(:,thisIndA,icell),2)),'r');
    hold on
    thisIndV = cycV_ind{icyc};
    plot(squeeze(mean(thisData(:,thisIndV,icell),2)),'g');
    hold on
    axis([0 100 -0.02 0.15]);
    for i = 1:4
    v = 5+ (i-1)*(dataStructVar.ONfr+dataStructVar.OFFfr);
    vline(v,'k');
    end
    t = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(i);
    vline(t,'c');
    hold on
    clear thisIndA thisDataA thisIndV thisDataV
end

%respond to visual target
figure;
nVCells = size(dataTrialsCells_VResp{1},3);

for icell = 1:nVCells
    subplot(4,4,icell);
    thisIndA = cycA_ind{4};
    thisData = dataTrialsCells_VResp{4};
    plot(squeeze(mean(thisData(:,thisIndA,icell),2)),'r');
    hold on
    thisIndV = cycV_ind{icyc};
    plot(squeeze(mean(thisData(:,thisIndV,icell),2)),'g');
    hold on
    axis([0 100 -0.02 0.15]);
    for i = 1:4
    v = 5+ (i-1)*(dataStructVar.ONfr+dataStructVar.OFFfr);
    vline(v,'k');
    end
    t = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(i);
    vline(t,'c');
    hold on
    clear thisIndA thisDataA thisIndV thisDataV
end

%respond to both
figure;
nAVCells = size(dataTrialsCells_AandVResp{1},3);

for icell = 1:nAVCells
    subplot(4,4,icell);
    thisIndA = cycA_ind{4};
    thisData = dataTrialsCells_AandVResp{4};
    plot(squeeze(mean(thisData(:,thisIndA,icell),2)),'r');
    hold on
    thisIndV = cycV_ind{icyc};
    plot(squeeze(mean(thisData(:,thisIndV,icell),2)),'g');
    hold on
    axis([0 100 -0.02 0.15]);
    for i = 1:4
    v = 5+ (i-1)*(dataStructVar.ONfr+dataStructVar.OFFfr);
    vline(v,'k');
    end
    t = 5+(dataStructVar.ONfr+dataStructVar.OFFfr).*(i);
    vline(t,'c');
    hold on
    clear thisIndA thisDataA thisIndV thisDataV
end
    
    
    

