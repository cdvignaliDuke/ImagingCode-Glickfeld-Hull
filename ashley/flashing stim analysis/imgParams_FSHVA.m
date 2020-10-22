%% common indices
visualTrials = 1;
auditoryTrials = 2;
allTrialsInd = 3;
alignStart = 1;
alignFA = 2;
alignCR = 3;
alignTarget = 4;
hitTrials = 1;
missTrials = 2;
oriBins = [0 1 32 90];
ampBins = [0 0.0001 0.05 1];
nStimBins = 3;
trOutNames = {'h','m','fa','cr'};
%% common variables
frameRateHz = 30;
nVisDelayFr = 3;
nVisDelayFr_target = 2;

%% testing
eyeAlpha = 0.05;
cellGroupsAlpha = 0.01;
tuningReliabilityThresh = 30;
tuningReliabilityThresh_decode = 11;
minRespThreshold = 0.002;
minRespThreshold_decode = 0.005;
minTrN = 5;
maxCellN = 15;
minCellN = 7;
minTrN_mdl = 20;
%% colors
cueColor = {[0 0 0];[.5 .5 1]};
AVColor = {[0 0 0];[0.5 0.5 0.5]};
hiLoColor = {[0.5 0.5 0.5];[0 0 0]};
aurocColor = {[0 0 0.75];[0.75 0 0];[0 0 0]};
taskTuneColor = {[0 0 0.75];[0.94 0.23 0.17];[0.6 0 0.05];[0 0 0]};
trOutColor = {[0 0 0];[0.75 0 0];[0.5 0.5 0.5];[0.75 0 0.75]};
%% example cells
if strcmp(ds,'FSAV_attentionV1')
    exampleCell_1 = 429; %418; %738 % first-stim responsive
    exampleCell_2 = 1297;%1269;%1269 % late responsive
    exampleCell_3 = 127;%386; % late suppressed
    
%     attnExCell_1 = 366; %first-stim responsive and modulated by attention(+V)
% elseif strcmp(ds,'FSAV_V1_100ms_naive')    
%     attnExCell_1 = 506;%603; %first-stim responsive 

    
end

%% eye tracking
eyeMmPerPix = 1/26.6; %mm/pix

preEventMs_eye = 1000; %ms
postEventMs_eye = 4500; %ms

shortTrialTimeS_eye = 1; %s
longTrialTimeS_eye = 2.9; %s
targetTimeS_eye = 1; %s

eyeBLFr = 16:30;
eyeRespWinFr = 48:54; % aligned with preEventMs_eye
eyeLateRespWinFr = 74:117; % aligned with preEventMs_eye

%% imaging align
preEventMs = 1000; %ms
postEventMs = 4500; %ms
longTrialLengthFr = 88;

area_list = {'LM','AL','PM'};
basewin = 1:34;
%LM
basewin_0{1} = 32:34;
basewin_0_target{1} = 30:32; 
respwin{1} = 36:38;
respwin_target{1} = 35:37;

%AL
basewin_0{2} = 32:34;
basewin_0_target{2} = 30:32; 
respwin{2} = 36:38;
respwin_target{2} = 35:37;

%PM
basewin_0{3} = 33:35;
basewin_0_target{3} = 31:33; 
respwin{3} = 37:39;
respwin_target{3} = 36:38;
