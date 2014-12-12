%set current directory for saving
%create folders
ImgFolder = '008+009'
CD = ['Z:\2P imaging\Analysis\' mouse '\' date ];
cd(CD);
mkdir(ImgFolder);
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder];
cd(CD);
mkdir('FlashingStimAnalysis','Rsp2VisStim');
%save path
Save = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2VisStim'];
cd(Save)

%% Concatinate preLeverDownplus3visstim
SubNum = '516';
date = '141003';

ImgFolder = '002';
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2VisStim'];
cd(CD);
preLDP3VS002 = readtiff('Rsp2VisStim.tif');

ImgFolder = '003';
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2VisStim'];
cd(CD);
preLDP3VS03 = readtiff('Rsp2VisStim.tif');

ImgFolder = '004';
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2VisStim'];
cd(CD);
preLDP3VS04 = readtiff('Rsp2VisStim.tif');

preLeverDownplus3visstim = cat(3,preLDP3VS002,preLDP3VS03,preLDP3VS04);
clear preLDP3VS002 preLDP3VS03 preLDP3VS04
%% Concatenate mworks variables

%first dataset
SubNum = '604';
date = '141010';
time = '1201';
% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

cLeverDown = double(cell2mat(input.cLeverDown));
nTrials = (input.trialSinceReset)-1;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
Block2ON = Block2ON(:,end-1);
cLeverDown = double(cell2mat(input.cLeverDown));
cTargetOn = input.cTargetOn;
for itrial = 1:nTrials
    if isempty(cTargetOn{itrial})
        cTargetOn{itrial} = NaN;
    end
end
cTargetOn = (double(cell2mat_padded(cTargetOn)))'; %For now NaNs == 0, may need to change...
cLeverUp = double(cell2mat(input.cLeverUp));
tCyclesOn = double(cell2mat(input.tCyclesOn));
TrialOutcome = input.trialOutcomeCell;


RateFRperMS = 30./1000;
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;


%second dataset
SubNum = '604';
date = '141010';
time = '1207';
% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

cLeverDown2 = double(cell2mat(input.cLeverDown));
nTrials2 = (input.trialSinceReset)-1;
Block2ON2 = double(cell2mat(input.tBlock2TrialNumber));
Block2ON2 = Block2ON2(:,end-1);
cLeverDown2 = double(cell2mat(input.cLeverDown));
cTargetOn2 = input.cTargetOn;
for itrial = 1:nTrials2
    if isempty(cTargetOn2{itrial})
        cTargetOn2{itrial} = NaN;
    end
end
cTargetOn2 = (double(cell2mat_padded(cTargetOn2)))'; %For now NaNs == 0, may need to change...
cLeverUp2 = double(cell2mat(input.cLeverUp));
tCyclesOn2 = double(cell2mat(input.tCyclesOn));
TrialOutcome2 = input.trialOutcomeCell;

%third dataset
SubNum = '516';
date = '141003';
time = '1203';
% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

cLeverDown3 = double(cell2mat(input.cLeverDown));
nTrials3 = (input.trialSinceReset);
Block2ON3 = double(cell2mat(input.tBlock2TrialNumber));
cLeverDown3 = double(cell2mat(input.cLeverDown));
cTargetOn3 = input.cTargetOn;
for itrial = 1:nTrials3
    if isempty(cTargetOn3{itrial})
        cTargetOn3{itrial} = NaN;
    end
end
cTargetOn3 = (double(cell2mat_padded(cTargetOn3)))'; %For now NaNs == 0, may need to change...
cLeverUp3 = double(cell2mat(input.cLeverUp));
tCyclesOn3 = double(cell2mat(input.tCyclesOn));
TrialOutcome3 = input.trialOutcomeCell;

%concatinate
cLeverDown = cat(2,cLeverDown,cLeverDown2);
nTrials = nTrials+nTrials2;
Block2ON = cat(2,Block2ON,Block2ON2);
cTargetOn = cat(2,cTargetOn,cTargetOn2);
cLeverUp = cat(2,cLeverUp,cLeverUp2);
tCyclesOn = cat(2,tCyclesOn,tCyclesOn2);
TrialOutcome = cat(2,TrialOutcome,TrialOutcome2);
Cycles = unique(tCyclesOn);
clear cLeverDown2 cLeverDown3 nTrials2 nTrials3 Block2ON2 Block2ON3 cTargetOn2 cTargetOn3 cLeverUp2 cLeverUp3 tCyclesOn2 tCyclesOn3 TrialOutcome2 TrialOutcome3 

SuccessTrials_log = strcmp(TrialOutcome,'success');
SuccessTrials_ind = find(SuccessTrials_log == 1);
IgnoreTrials_log = strcmp(TrialOutcome,'ignore');
IgnoreTrials_ind = find(IgnoreTrials_log ==1);
SuccessANDIgnoreTrials_ind = sort([SuccessTrials_ind IgnoreTrials_ind]);
