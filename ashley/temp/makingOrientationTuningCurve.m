%% 

% expt info
responseWindow = 13:18;
baselineWindow = 9:12;
downSampleRate = 10;

% load time-courses
clear all
FSAV_V1_100ms_naive_temp;
rc = behavConstsAV;
iexp = 1;

ms = expt(iexp).mouse;
dt = expt(iexp).date;

fn = fullfile(rc.ashleyAnalysis,ms,'two-photon imaging', dt, expt(iexp).dirtuning);

load(fullfile(fn,'timecourses_tun_cells.mat'))

% load vis stim info
visStimInfoName = ['data-i' ms '-' dt '-' expt(iexp).dirtuning_time];
load(fullfile(rc.behavData,visStimInfoName));
visStimInfo = input; clear input

% align neuronal data to the behavior/visual stim data
nOn = visStimInfo.nScansOn/downSampleRate;
nOff = visStimInfo.nScansOff/downSampleRate;
trialLengthFrames = nOn+nOff;

[nFrames,nCells] = size(data_tun_tc_subnp);

nTrials = nFrames/trialLengthFrames;

trialStarts = 1:trialLengthFrames:(nTrials*trialLengthFrames);

F_trials = nan(trialLengthFrames,nTrials,nCells);
for i = 1:nTrials
    ind = trialStarts(i):(trialStarts(i)+trialLengthFrames-1);
    F_trials(:,i,:) = data_tun_tc_subnp(ind,:);
end

% get each trial dF/F

F = mean(F_trials(baselineWindow,:,:),1);
dF = F_trials - F;
dFF = dF./F;

% find average responses to specific stim
trialDirections = cell2mat(visStimInfo.tGratingDirectionDeg);
directions = unique(trialDirections);



% generate tuning curves
% fit tuning curves
% quantify tuning - how well fit? OSI, DSI? 

