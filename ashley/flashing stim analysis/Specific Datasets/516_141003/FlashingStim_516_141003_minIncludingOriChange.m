
%% 002
SubNum = '516';
date = '141003';
time = '1156';
ImgFolder = '002';
mouse = '516';
fName = '002_000_000';

% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

CD = ('Z:\2P imaging\Analysis\516\141003\002\FlashingStimAnalysis');
cd(CD);
data = readtiff('FSall.tif');

% CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'];
% cd(CD);
% mkdir('CmprTrialLength');
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

nTrials = input.trialSinceReset;
nTrials_mat = 1:nTrials;
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
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
Cycles = unique(tCyclesOn);
minCyclesOn = input.minCyclesOn;
maxCyclesOn = input.maxCyclesOn;

% save dataset of just min cycles on and max cycles on (success only)
success_log = strcmp(TrialOutcome,'success');
success_ind = find(success_log == 1);
%min dataset
minCycles_ind = find(tCyclesOn == minCyclesOn);
minCycles_trials = intersect(success_ind,minCycles_ind);
cLeverDown_minCycles = cLeverDown(minCycles_trials);
cLeverUp_minCycles = cLeverUp(minCycles_trials);
cTargetOn_minCycles = cTargetOn(minCycles_trials);

tTrials = size(minCycles_trials,2);
LminPlus = 10+min(cLeverUp_minCycles-cLeverDown_minCycles)+10;

minCyclesTrialsPlus = zeros(size(data,1),size(data,2),tTrials*LminPlus);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_minCycles(itrial)-10:cLeverDown_minCycles(itrial)+LminPlus-11;
    minCyclesTrialsPlus(:,:,start:start+LminPlus-1) = data(:,:,ind);
    start = start+LminPlus;
end

writetiff(minCyclesTrialsPlus,'minCyclesTrialsPlus.tif');

CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
cd(CD);
save ('LminPlus.mat','LminPlus');
tTrialsMinAll = tTrials;
save ('tTrialsMinAll.mat','tTrialsMinAll');
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

clear all

%% 003

SubNum = '516';
date = '141003';
time = '1159';
ImgFolder = '003';
mouse = '516';
fName = '003_000_000';

% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

CD = ('Z:\2P imaging\Analysis\516\141003\003\FlashingStimAnalysis');
cd(CD);
data = readtiff('FSall_003.tif');

% CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'];
% cd(CD);
% mkdir('CmprTrialLength');
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);
% name mworks variables
nTrials = input.trialSinceReset;
nTrials_mat = 1:nTrials;
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
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
Cycles = unique(tCyclesOn);
minCyclesOn = input.minCyclesOn;
maxCyclesOn = input.maxCyclesOn;

success_log = strcmp(TrialOutcome,'success');
success_ind = find(success_log == 1);
%min dataset
minCycles_ind = find(tCyclesOn == minCyclesOn);
minCycles_trials = intersect(success_ind,minCycles_ind);
cLeverDown_minCycles = cLeverDown(minCycles_trials);
cLeverUp_minCycles = cLeverUp(minCycles_trials);
cTargetOn_minCycles = cTargetOn(minCycles_trials);

tTrials = size(minCycles_trials,2);
CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
cd(CD);
load ('LminPlus.mat');

minCyclesTrialsPlus = zeros(size(data,1),size(data,2),tTrials*LminPlus);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_minCycles(itrial)-10:cLeverDown_minCycles(itrial)+LminPlus-11;
    minCyclesTrialsPlus(:,:,start:start+LminPlus-1) = data(:,:,ind);
    start = start+LminPlus;
end

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\002\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);
writetiff(minCyclesTrialsPlus,'minCyclesTrialsPlus.tif');

clear all

%% 004

SubNum = '516';
date = '141003';
time = '1203';
ImgFolder = '004';
mouse = '516';
fName = '004_000_000';

% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

CD = ('Z:\2P imaging\Analysis\516\141003\004\FlashingStimAnalysis');
cd(CD);
data = readtiff('FSall_004.tif');

% CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'];
% cd(CD);
% mkdir('CmprTrialLength');
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);
% name mworks variables
nTrials = input.trialSinceReset;
nTrials_mat = 1:nTrials;
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
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
Cycles = unique(tCyclesOn);
minCyclesOn = input.minCyclesOn;
maxCyclesOn = input.maxCyclesOn;

success_log = strcmp(TrialOutcome,'success');
success_ind = find(success_log == 1);
%min dataset
minCycles_ind = find(tCyclesOn == minCyclesOn);
minCycles_trials = intersect(success_ind,minCycles_ind);
cLeverDown_minCycles = cLeverDown(minCycles_trials);
cLeverUp_minCycles = cLeverUp(minCycles_trials);
cTargetOn_minCycles = cTargetOn(minCycles_trials);

tTrials = size(minCycles_trials,2);
CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
cd(CD);
load ('LminPlus.mat');

minCyclesTrialsPlus = zeros(size(data,1),size(data,2),tTrials*LminPlus);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_minCycles(itrial)-10:cLeverDown_minCycles(itrial)+LminPlus-11;
    minCyclesTrialsPlus(:,:,start:start+LminPlus-1) = data(:,:,ind);
    start = start+LminPlus;
end

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);
writetiff(minCyclesTrialsPlus,'minCyclesTrialsPlus.tif');

clear all
%% cat

%%002
SubNum = '516';
date = '141003';
time = '1156';
ImgFolder = '002';
mouse = '516';
fName = '002_000_000';

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);
minCyclesTrials002 = readtiff('minCyclesTrialsPlus.tif');

%003
SubNum = '516';
date = '141003';
time = '1159';
ImgFolder = '003';
mouse = '516';
fName = '003_000_000';

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);
minCyclesTrials003 = readtiff('minCyclesTrialsPlus.tif');

%004
SubNum = '516';
date = '141003';
time = '1203';
ImgFolder = '004';
mouse = '516';
fName = '004_000_000';

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);
minCyclesTrials004 = readtiff('minCyclesTrialsPlus.tif');

%cat
SubNum = '516';
date = '141003';
mouse = '516';
ImgFolder = '002+003+004';
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);

minCyclesTrialsAll = cat(3,minCyclesTrials002,minCyclesTrials003,minCyclesTrials004);
writetiff(minCyclesTrialsAll,'minCyclesTrialsPlusAll.tif');

clear all

