awFSAVdatasets
for iexp = 1:size(expt,2)
SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
date = expt(iexp).date;
runs = expt(iexp).runs;
nrun = size(runs,1);
expFolder = expt(iexp).folder;
time_mat = expt(iexp).time_mat;

%%
%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
runstr = runs(1,:);
if nrun>1
    for irun = 2:nrun
        runstr = [runstr '-' runs(irun,:)];
    end
end
% fnout = fullfile('Z:\home\lindsey\Analysis\2P', mouse, date, [date '_' mouse '_' runstr '_']);
fnout = ['Z:\Analysis\' mouse '\' expFolder '\' date];

%% load and combine mworks data and timecourses
input = [];
for irun = 1:nrun
    time = time_mat(irun,:);
    fn_mworks = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' SubNum '-' date '-' time '.mat'];
    if irun == 1
        input = mwLoadData(fn_mworks, [], []);
    else
        input = [input mwLoadData(fn_mworks, [], [])];
    end
end
input = concatenateDataBlocks(input);

%%
run_trials = input.trialsSinceReset;
cLeverDown = cell2mat(input.cLeverDown);
cLeverUp = cell2mat(input.cLeverUp);
cTargetOn = celleqel2mat_padded(input.cTargetOn);
cStimOn = celleqel2mat_padded(input.cStimOn);
cItiStart = cell2mat(input.cItiStart);
if expt(iexp).catch == 1;
    cCatchOn = celleqel2mat_padded(input.cCatchOn);
    isFA = celleqel2mat_padded(input.tFalseAlarm);
end

dataTC = [];
offset = 0;
for irun = 1:nrun
    ImgFolder = runs(irun,:);
    fnTC = fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fnTC);
    load('Timecourses.mat')
    dataTC = cat(1, dataTC, dataTimecourse.dataTCsub);
    offset = offset+size(dataTimecourse.dataTCsub,1);
    if irun < nrun
        startTrial = sum(run_trials(1, 1:irun),2)+1;
        endTrial = sum(run_trials(1,1:irun+1),2);
        cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
        cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
        cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
        cStimOn(1,startTrial:endTrial) = cStimOn(1,startTrial:endTrial)+offset;
        cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
        if expt(iexp).catch == 1;
            cCatchOn(1,startTrial:endTrial) = cCatchOn(1,startTrial:endTrial)+offset;
            isFA(1,startTrial:endTrial) = cCatchOn(1,startTrial:endTrial)+offset;
        end
    end
end

ntrials = length(input.trialOutcomeCell);
trialOutcome = input.trialOutcomeCell;
tCyclesOn = cell2mat(input.tCyclesOn);
nCyclesOn = cell2mat(input.nCyclesOn);

cycles = unique(tCyclesOn);
block2 = cell2mat(input.tBlock2TrialNumber);
V_ind = find(cell2mat(input.tBlock2TrialNumber) == 0);
AV_ind = find(cell2mat(input.tBlock2TrialNumber) == 1);
cycTime = input.nFramesOn+input.nFramesOff;
tooFastTime = input.nFramesTooFast;
maxReactTime = input.nFramesReact;

DirectionDeg = cell2mat(input.tGratingDirectionDeg);
Dirs = unique(DirectionDeg);
tGratingDirectionDeg = chop(cell2mat(input.tGratingDirectionDeg),4);

if expt(iexp).catch == 1;
    catchCycle = cell2mat(input.catchCyclesOn);
    catchDirectionDeg = cell2mat_padded(input.tCatchGratingDirectionDeg);
    catchDirs = unique(catchDirectionDeg);
    isCatchTrial = catchDirectionDeg > 0;

    catchTrialOutcome = num2cell(NaN(length(nCyclesOn),1));
    catchIndex = find(isCatchTrial == 1);
    for i = 1:sum(isCatchTrial)
        if isFA(catchIndex(i)) == 1
            catchTrialOutcome{catchIndex(i),1} = 'FA';
        elseif cCatchOn(catchIndex(i)) == 0
            catchTrialOutcome{catchIndex(i),1} = 'failure';
        elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) < tooFastTime
            catchTrialOutcome{catchIndex(i),1} = 'failure';
        elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) > maxReactTime
            catchTrialOutcome{catchIndex(i),1} = 'CR';
        end
    end

end

clear dataTimecourse

try
    cd([fnout '\' runstr])
catch
    mkdir(fnout,runstr)
end

save([mouse '-' expt(iexp).date '-' runstr '-comboInputDataTC.mat'],'dataTC','input');
save([mouse '-' expt(iexp).date '-' runstr '-comboInputDataTCplusVar.mat']);
clear catchCycle catchDirectionDeg catchDirs isCatchTrial catchTrialOutcome catchIndex
end
