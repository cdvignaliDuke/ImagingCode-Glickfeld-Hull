%% Load registered dataset, retinotopy, DSI, OSI, etc.
SubNum = '001';
date = '140828';
time = '009';
ImgFolder = '009';
mouse = 'AW01';
fName = '009_000_000';
experiment = 'Flashing Stim';

%load retinotopy
retfile = ['Z:\2P imaging\Analysis\' mouse '\' experiment '\' date '\Retinotopy_V1\analysis']; 
load (retfile);
clear data_reg

%load registered data, cell timecourses, and mworks variables
FSAfile = ['Z:\2P imaging\Analysis\' mouse '\' experiment '\' date '\FlashingStimAnalysis\analysis'];
load(FSAfile);

%set directory for saving 
Save = ['Z:\2P imaging\Analysis\' mouse '\' experiment '\' date '\FlashingStimAnalysis\Rsp2OriChange'];
cd(Save)
%% call some mworks variables

nTrials = input.trialSinceReset;
cLeverDown = double(cell2mat(input.cLeverDown));
cTargetOn = double(cell2mat(input.cTargetOn));
cLeverUp = double(cell2mat(input.cLeverUp));
tCyclesOn = double(cell2mat(input.tCyclesOn));
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
Cycles = unique(tCyclesOn);
ChangeMag = double(cell2mat(input.tGratingMaxDirectionStepDeg));
ChangeMag = ChangeMag(:,1:end-1);

%% Find off and on indices for each trial

start = 1;
start1 = 1;
for itrial = 1:nTrials
    start2 = cLeverDown(itrial)-start +1;
    nOFF_ind (1,start1:start1+start2-1) =  start:cLeverDown(itrial);
    start1 = start1+start2;
    start = cLeverUp(itrial);
    
end
siz = size(nOFF_ind,2);
nOFF_ind = nOFF_ind(:,1:siz-2);

%F1=average fluorescence for all iti (stim OFF)
F1 = mean(data_reg(:,:,nOFF_ind),3);

%% Find dFoverF for each trial for F1 iti

dF1_data = bsxfun(@minus,data_reg, F1);
dF1overF1_data = bsxfun(@rdivide, dF1_data, F1);
max_dF1 = max(dF1overF1_data,[],3);
figure; imagesq(max_dF1); colormap(gray)   

%% Get average response for each TargetON 100ms timeperiod

%per full field
for itrial = 1:nTrials-1
    x = cTargetOn(itrial):cTargetOn(itrial) + ceil(100*RateFRperMS);
    siz = size(x,2);
    TargetRspAll(1:siz,itrial) = dF1overF1all_mean(x);
end
TargetRspAll_mean = mean(TargetRspAll,1);

%per cell
for itrial = 1:nTrials-1
    x = cTargetOn(itrial):cTargetOn(itrial) + ceil(100*RateFRperMS);
    siz = size(x,2);
    for icell = 1:nCell
        TargetRspCell(1:siz,icell,itrial) = data_TC_FS(x,icell);
    end
end
TargetRspCell_mean = squeeze(mean(TargetRspCell,1));

%% Plot target response by degree of orientation change

nChangeMag = unique(ChangeMag);
Changes = size(nChangeMag,2);


%per full field
for ideg = 1:Changes
    ChangeMag_ind = find(ChangeMag==nChangeMag(ideg));
    dF1overF1all_byChangeMag(1,ideg) = mean(TargetRspAll_mean(ChangeMag_ind));
end
Rsp2TargetbyChangeMag = figure;
hold on
title('Full field response to target stimulus by degree of orientation change');
hold on
xlabel('Ori change');
hold on
ylabel('dF/F1');
hold on
plot(dF1overF1all_byChangeMag);
saveas(Rsp2TargetbyChangeMag, 'Rsp2TargetbyChangeMag.fig');

%per cell
for icell = 1:nCell
    for ideg = 1:Changes
        ChangeMag_ind = find(ChangeMag==nChangeMag(ideg));
        dF1overF1cell_byChangeMag(icell,ideg) = mean(TargetRspCell_mean(icell,ChangeMag_ind),2);
    end
end
%sort cells by amount of dF/F for lowest change in orientation
Cell_Rsp2TargetbyChangeMag = sortrows(dF1overF1cell_byChangeMag,1);

Cell_Rsp2TargetbyChangeMag_fig = figure;
colormap('hot');
imagesc(Cell_Rsp2TargetbyChangeMag);
c = colorbar;
hold on
title('Cell response to target stimulus by degree of orientation change');
hold on
xlabel('Ori change');
hold on
ylabel('Cell');
hold on
ylabel(c,'dF/F1');
saveas(Cell_Rsp2TargetbyChangeMag_fig,'Cell_Rsp2TargetbyChangeMag_fig.fig');