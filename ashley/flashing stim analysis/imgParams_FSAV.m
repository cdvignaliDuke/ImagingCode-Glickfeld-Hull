%% common indices
visualTrials = 1;
auditoryTrials = 2;
alignStart = 1;
alignTarget = 2;
%% common variables
frameRateHz = 30;
nVisDelayFr = 3;

%% testing
eyeAlpha = 0.05;
%% colors
cueColor = {[0 0 0];[.5 .5 1]};
AVColor = {[0 0 0];[0.5 0.5 0.5]};
hiLoColor = {[0.5 0.5 0.5];[0 0 0]};

%% example cells
if strcmp(ds,'FSAV_attentionV1')
    exampleCell_1 = 418; %738 % first-stim responsive
    exampleCell_2 = 1223;%1269;%1269 % late responsive
    exampleCell_3 = 543;%386; % late suppressed
    
    attnExCell_1 = 366; %first-stim responsive and modulated by attention(+V)
elseif strcmp(ds,'FSAV_V1_100ms_naive')    
    attnExCell_1 = 506;%603; %first-stim responsive 
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