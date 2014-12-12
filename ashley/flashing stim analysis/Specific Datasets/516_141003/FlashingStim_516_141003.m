%% **002**

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
Lmin = 10+min(cLeverUp_minCycles-cLeverDown_minCycles);

minCyclesTrials = zeros(size(data,1),size(data,2),tTrials*Lmin);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_minCycles(itrial)-10:cLeverDown_minCycles(itrial)+Lmin-11;
    minCyclesTrials(:,:,start:start+Lmin-1) = data(:,:,ind);
    start = start+Lmin;
end

writetiff(minCyclesTrials,'minCyclesTrials.tif');

CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
cd(CD);
save ('Lmin.mat','Lmin');
tTrialsMinAll = tTrials;
save ('tTrialsMinAll.mat','tTrialsMinAll');
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

clear tTrials 

%max dataset
maxCycles_ind = find(tCyclesOn == maxCyclesOn);
maxCycles_trials = intersect(success_ind,maxCycles_ind);
cLeverDown_maxCycles = cLeverDown(maxCycles_trials);
cLeverUp_maxCycles = cLeverUp(maxCycles_trials);
cTargetOn_maxCycles = cTargetOn(maxCycles_trials);
if isempty(maxCycles_trials)
    maxCycles_ind = find(tCyclesOn == (maxCyclesOn-1));
    maxCycles_trials = intersect(success_ind,maxCycles_ind);
    cLeverDown_maxCycles = cLeverDown(maxCycles_trials);
    cLeverUp_maxCycles = cLeverUp(maxCycles_trials);
    cTargetOn_maxCycles = cTargetOn(maxCycles_trials);
    trueMax = 0;
    Lmax = 10+min(cLeverUp_maxCycles-cLeverDown_maxCycles);
else
    x = 'max is true max';
    trueMax = 1;
    disp(x);
    Lmax = 10+min(cLeverUp_maxCycles-cLeverDown_maxCycles);
    CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
    cd(CD);
    save ('Lmax.mat','Lmax');
    CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
    cd(CD);
end

tTrials = size(maxCycles_trials,2);
if trueMax == 1
    CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
    cd(CD);
    tTrialsMaxAll = tTrials;
    save ('tTrialsMaxAll.mat','tTrialsMaxAll');
    CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
    cd(CD);
else trueMax == 0
    CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
    cd(CD);
    tTrialsMaxAll = 0;
    save ('tTrialsMaxAll.mat','tTrialsMaxAll');
    CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
    cd(CD);   
end

maxCyclesTrials = zeros(size(data,1),size(data,2),tTrials*Lmax);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_maxCycles(itrial)-10:cLeverDown_maxCycles(itrial)+Lmax-11;
    maxCyclesTrials(:,:,start:start+Lmax-1) = data(:,:,ind);
    start = start+Lmax;
end

if trueMax == 1
    writetiff(maxCyclesTrials,'maxCyclesTrials.tif');
elseif trueMax == 0
    writetiff(maxCyclesTrials,'maxCyclesTrials_2ndHighest.tif');
end

clear all

%% **003**
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
load ('Lmin.mat');
load ('tTrialsMinAll.mat');
tTrialsMinAll = tTrials+tTrialsMinAll;
save ('tTrialsMinAll.mat','tTrialsMinAll');
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

minCyclesTrials = zeros(size(data,1),size(data,2),tTrials*Lmin);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_minCycles(itrial)-10:cLeverDown_minCycles(itrial)+Lmin-11;
    minCyclesTrials(:,:,start:start+Lmin-1) = data(:,:,ind);
    start = start+Lmin;
end

writetiff(minCyclesTrials,'minCyclesTrials.tif');

clear tTrials L

%max dataset
maxCycles_ind = find(tCyclesOn == maxCyclesOn);
maxCycles_trials = intersect(success_ind,maxCycles_ind);
cLeverDown_maxCycles = cLeverDown(maxCycles_trials);
cLeverUp_maxCycles = cLeverUp(maxCycles_trials);
cTargetOn_maxCycles = cTargetOn(maxCycles_trials);


CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
cd(CD);
try 
    load('Lmax.mat')
    if isempty(maxCycles_trials)
        maxCycles_ind = find(tCyclesOn == (maxCyclesOn-1));
        maxCycles_trials = intersect(success_ind,maxCycles_ind);
        cLeverDown_maxCycles = cLeverDown(maxCycles_trials);
        cLeverUp_maxCycles = cLeverUp(maxCycles_trials);
        cTargetOn_maxCycles = cTargetOn(maxCycles_trials);
        trueMax = 0;
    else
        x = 'max is true max';
        trueMax = 1;
        disp(x);
    end
catch 
    if isempty(maxCycles_trials)
        maxCycles_ind = find(tCyclesOn == (maxCyclesOn-1));
        maxCycles_trials = intersect(success_ind,maxCycles_ind);
        cLeverDown_maxCycles = cLeverDown(maxCycles_trials);
        cLeverUp_maxCycles = cLeverUp(maxCycles_trials);
        cTargetOn_maxCycles = cTargetOn(maxCycles_trials);
        trueMax = 0;
        Lmax = 10+min(cLeverUp_maxCycles-cLeverDown_maxCycles);
    else
        x = 'max is true max';
        trueMax = 1;
        disp(x);
        Lmax = 10+min(cLeverUp_maxCycles-cLeverDown_maxCycles);
        save ('Lmax.mat','Lmax');
    end
end

tTrials = size(maxCycles_trials,2);
load('tTrialsMaxAll.mat');
if trueMax == 1
    tTrialsMaxAll = tTrials+tTrialsMaxAll;
    save ('tTrialsMaxAll.mat','tTrialsMaxAll');
else trueMax == 0
    x = 'no new trials added to max';
    disp(x)
end

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

maxCyclesTrials = zeros(size(data,1),size(data,2),tTrials*Lmax);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_maxCycles(itrial)-10:cLeverDown_maxCycles(itrial)+Lmax-11;
    maxCyclesTrials(:,:,start:start+Lmax-1) = data(:,:,ind);
    start = start+Lmax;
end

if trueMax == 1
    writetiff(maxCyclesTrials,'maxCyclesTrials.tif');
elseif trueMax == 0
    writetiff(maxCyclesTrials,'maxCyclesTrials_2ndHighest.tif');
end

clear all

%% **004**
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
load ('Lmin.mat');
load ('tTrialsMinAll.mat');
tTrialsMinAll = tTrials+tTrialsMinAll;
save ('tTrialsMinAll.mat','tTrialsMinAll');
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

minCyclesTrials = zeros(size(data,1),size(data,2),tTrials*Lmin);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_minCycles(itrial)-10:cLeverDown_minCycles(itrial)+Lmin-11;
    minCyclesTrials(:,:,start:start+Lmin-1) = data(:,:,ind);
    start = start+Lmin;
end

writetiff(minCyclesTrials,'minCyclesTrials.tif');

clear tTrials

%max dataset
maxCycles_ind = find(tCyclesOn == maxCyclesOn);
maxCycles_trials = intersect(success_ind,maxCycles_ind);
cLeverDown_maxCycles = cLeverDown(maxCycles_trials);
cLeverUp_maxCycles = cLeverUp(maxCycles_trials);
cTargetOn_maxCycles = cTargetOn(maxCycles_trials);


CD = ('Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis');
cd(CD);

try 
    load('Lmax.mat')
catch 
    if isempty(maxCycles_trials)
        maxCycles_ind = find(tCyclesOn == (maxCyclesOn-1));
        maxCycles_trials = intersect(success_ind,maxCycles_ind);
        cLeverDown_maxCycles = cLeverDown(maxCycles_trials);
        cLeverUp_maxCycles = cLeverUp(maxCycles_trials);
        cTargetOn_maxCycles = cTargetOn(maxCycles_trials);
        trueMax = 0;
        Lmax = 10+min(cLeverUp_maxCycles-cLeverDown_maxCycles);
    else
        x = 'max is true max';
        trueMax = 1;
        disp(x);
        Lmax = 10+min(cLeverUp_maxCycles-cLeverDown_maxCycles);
        save ('Lmax.mat','Lmax');
    end
end
tTrials = size(maxCycles_trials,2);
load('tTrialsMaxAll.mat');
if trueMax == 1
    tTrialsMaxAll = tTrials+tTrialsMaxAll;
    save ('tTrialsMaxAll.mat','tTrialsMaxAll');
else trueMax == 0
    x = 'no new trials added to max'
    disp(x)
end
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

maxCyclesTrials = zeros(size(data,1),size(data,2),tTrials*Lmax);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_maxCycles(itrial)-10:cLeverDown_maxCycles(itrial)+Lmax-11;
    maxCyclesTrials(:,:,start:start+Lmax-1) = data(:,:,ind);
    start = start+Lmax;
end

if trueMax == 1
    writetiff(maxCyclesTrials,'maxCyclesTrials.tif');
elseif trueMax == 0
    writetiff(maxCyclesTrials,'maxCyclesTrials_2ndHighest.tif');
end

clear all

%% **cat datasets**
%%002
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

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

input002 = input; clear input
minCyclesTrials002 = readtiff('minCyclesTrials.tif');
maxCyclesTrials002 = NaN;

%003
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

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

input003 = input; clear input
minCyclesTrials003 = readtiff('minCyclesTrials.tif');
maxCyclesTrials003 = readtiff('maxCyclesTrials.tif');

%004
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

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
cd(CD);

input004 = input; clear input
minCyclesTrials004 = readtiff('minCyclesTrials.tif');
maxCyclesTrials004 = readtiff('maxCyclesTrials.tif');

%cat
SubNum = '516';
date = '141003';
mouse = '516';
ImgFolder = '002+003+004';
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);
% mkdir('CmprTrialLength');
% CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\CmprTrialLength'];
% cd(CD);

% minCyclesTrialsAll = cat(3,minCyclesTrials002,minCyclesTrials003,minCyclesTrials004);
% writetiff(minCyclesTrialsAll,'minCyclesTrialsAll.tif');
% maxCyclesTrialsAll = cat(3, maxCyclesTrials003,maxCyclesTrials004);
% writetiff(maxCyclesTrialsAll,'maxCyclesTrialsAll.tif');
% 
% inputAll = [input002 input003 input004];
% save('mWorksAll.mat', 'inputAll');
minCyclesTrialsAll = readtiff('minCyclesTrialsAll.tif');
maxCyclesTrialsAll = readtiff('maxCyclesTrialsAll.tif');
load('mWorksAll.mat');
load('tTrialsMinAll.mat');
load('tTrialsMaxAll.mat');
load('mask&TCMin.mat');
load('mask&TCMax.mat');
load('Lmax.mat');
load('Lmin.mat');

%% **put all datasets together - just successful trials**
% **002**

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
successCycles = minCyclesOn:maxCyclesOn;
success_log = strcmp(TrialOutcome,'success');
success_ind = find(success_log == 1);

%5 cylcles dataset
Cycles5_ind = find(tCyclesOn == minCyclesOn);
Cycles5_trials = intersect(success_ind,Cycles5_ind);
cLeverDown_Cycles5 = cLeverDown(Cycles5_trials);
cLeverUp_Cycles5 = cLeverUp(Cycles5_trials);
cTargetOn_Cycles5 = cTargetOn(Cycles5_trials);

tTrials = size(Cycles5_trials,2);
L5 = 10+min(cLeverUp_Cycles5-cLeverDown_Cycles5);

Cycles5Trials = zeros(size(data,1),size(data,2),tTrials*L5);
start = 1;
for itrial = 1:tTrials
    ind = cLeverDown_Cycles5(itrial)-10:cLeverDown_Cycles5(itrial)+L5-11;
    Cycles5Trials(:,:,start:start+L5-1) = data(:,:,ind);
    start = start+L5;
end
