%% using dataStruct, plot average trace for all trials, sucess, early, and miss 
date = '140923';
mouse = 'AW04';
ImgFolder = '005+006+007';       
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
%% dF/F
siz = size(dataStruct.Cycles,2);
for icyc = 1:siz
    for itrial = 1:dataStruct.cycNTrials{icyc}
        thisCycData(:,:,itrial) = mean(dataStruct.cycData{icyc}(:,:,(dataStruct.cycTrialL{icyc}.*(itrial-1)+1):(dataStruct.cycTrialL{icyc}.*(itrial-1)+1)+5),3);
    end
    cycF{1,icyc} = thisCycData;
    clear thisCycData
end

for icyc = 1:siz
    for itrial = 1:dataStruct.cycNTrials{icyc}
        start1 = dataStruct.cycTrialL{icyc}.*(itrial-1)+1;
        start2 = dataStruct.cycTrialL{icyc}.*(itrial);
        thisData = dataStruct.cycData{icyc};
        thisDataF = cycF{icyc};
        thisCycData(:,:,start1:start2) = bsxfun(@minus,thisData(:,:,start1:start2),thisDataF(:,:,itrial));
    end
    cycDF{1,icyc} = thisCycData;
    clear thisCycData thisData thisDataF
end

for icyc = 1:siz
    for itrial = 1:dataStruct.cycNTrials{icyc}
        start1 = dataStruct.cycTrialL{icyc}.*(itrial-1)+1;
        start2 = dataStruct.cycTrialL{icyc}.*(itrial);
        thisDataDF = cycDF{icyc};
        thisDataF = cycF{icyc};
        thisCycData(:,:,start1:start2) = bsxfun(@rdivide,thisDataDF(:,:,start1:start2),thisDataF(:,:,itrial));
    end
    cycDFoverF{1,icyc} = thisCycData;
    clear thisCycData thisDataDF thisDataF
end

dataStructDFoverF.cycData = cycDFoverF;
dataStructDFoverF.cycTrialL = dataStruct.cycTrialL; 
dataStructDFoverF.cycNTrials = dataStruct.cycNTrials;
dataStructDFoverF.cycLeverDown = dataStruct.cycLeverDown;
dataStructDFoverF.cycTargetOn = dataStruct.cycTargetOn;
dataStructDFoverF.cycLeverOn = dataStruct.cycLeverOn;
dataStructDFoverF.cycBlock2ON = dataStruct.cycBlock2ON;
dataStructDFoverF.cycTrialOutcome = dataStruct.cycTrialOutcome;
dataStructDFoverF.ONms = dataStruct.ONms;
dataStructDFoverF.OFFms = dataStruct.OFFms;
dataStructDFoverF.RateFRperMS = dataStruct.RateFRperMS;
dataStructDFoverF.Cycles = dataStruct.Cycles;
dataStructDFoverF.minCyclesOn = dataStruct.minCyclesOn;
dataStructDFoverF.maxCyclesOn = dataStruct.maxCyclesOn;
dataStructDFoverF.mouse = dataStruct.mouse;
dataStructDFoverF.date = dataStruct.date;
dataStructDFoverF.ImgFolder = dataStruct.ImgFolder;


%max image for each cycle

for icyc = 1:siz
    thisCycleData = dataStructDFoverF.cycData{icyc};
    thisMax = max(thisCycleData,[],3);
    maxDFoverF{icyc} = thisMax;
    clear thisMax thisCycleData
end

for icyc = 1:siz
    FIG = maxDFoverF{icyc};
    figure; imagesq(FIG); colormap(gray)
end

% bwout = imCellEditInteractive(max_dF_VisRspMin);
% mask_cellMin = bwlabel(bwout);


%get data_TC for each cycle using direction tuning mask

CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\003 - Direction Selectivity'];
cd(CD);
load('mask_cell_DirTuning.mat');
mask_cellDir = mask_cell;
clear mask_cell

for icyc = 1:siz
    thisCycleData = dataStructDFoverF.cycData{icyc};
    thisTC = stackGetTimeCourses(thisCycleData, mask_cellDir);
    cellTC{icyc} = thisTC;
    nCells  = size(thisTC,2);
end

for icyc = 1:siz
    FIG = cellTC{icyc};
    figure; tcOffsetPlot(FIG); colormap(gray)
end

dataStructDFoverF.maxDFoverF = maxDFoverF;
dataStructDFoverF.cellTC = cellTC;
dataStructDFoverF.nCells = nCells;
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
save('dataStructDFoverF.mat','dataStructDFoverF');

