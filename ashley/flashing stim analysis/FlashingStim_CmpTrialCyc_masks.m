mouse = 'AW04';
SubNum = '604';
date = '140923';
ImgFolder = '005+006+007';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructDFoverF.mat')

%%
% create structure dataMask for max dF/F for datasets: direction tuning,
% flashing stim (first flash rsp), flashing stim rsp to target (for each
% cycle)
% parameters for imCellEditInteractive are 0.9 and 3 pixels
siz = size(dataStructDFoverF.Cycles,2);

% direction tuning dataset
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\003 - Direction Selectivity'];
cd(CD);
load('mask_cell_DirTuning.mat');
dirTuning = mask_cell;
clear mask_cell

dataMasks.dirTuning = dirTuning;

% flashing stim datasets
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

for icyc = 1:siz
    thisMax = dataStructDFoverF.maxDFoverF{icyc};
    bwout = imCellEditInteractive(thisMax);
    thisMask = bwlabel(bwout);
    FScycMaxDF{icyc} = thisMask;
    clear thisMask thisMax bwout
end

dataMasks.FScycMaxDF = FScycMaxDF;

%flashing stim - average first response
for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    thisL = dataStructDFoverF.cycTrialL{icyc};
    start = 1;
    for itrial = 1:dataStructDFoverF.cycNTrials{icyc}
        thisDatabyTrial = thisData(:,:,(1+(thisL.*(itrial-1)):thisL.*itrial));
        firstsData(:,:,start:start+29) = thisDatabyTrial(:,:,1:30);
        start = start+30;
    end
    maxDFoverF(:,:,icyc) = max(firstsData,[],3);
    clear thisData thisL thisDatabyTrial firstsData
end

FSfirstRspMax = max(maxDFoverF,[],3);
clear maxDFoverF
bwout = imCellEditInteractive(FSfirstRspMax);
FSfirstRsp = bwlabel(bwout);

dataMasks.FSfirstRsp = FSfirstRsp;

%flashing stim - target response
for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    thisL = dataStructDFoverF.cycTrialL{icyc};
    start = 1;
    for itrial = 1:dataStructDFoverF.cycNTrials{icyc}
        thisDatabyTrial = thisData(:,:,(1+(thisL.*(itrial-1)):thisL.*itrial));
        targetData(:,:,start:start+49) = thisDatabyTrial(:,:,thisL-49:thisL);
        start = start+50;
    end
    cycTargetMax{icyc} = max(targetData,[],3);
    clear thisData thisL thisDatabyTrial targetData
end

for icyc = 1:siz
    thisMax = cycTargetMax{icyc};
    bwout = imCellEditInteractive(thisMax);
    thisMask = bwlabel(bwout);
    FScycTargetRsp{icyc} = thisMask;
    clear thisMask thisMax bwout
end

dataMasks.FScycTargetRsp = FScycTargetRsp;

CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
save('dataMasks.mat','dataMasks')