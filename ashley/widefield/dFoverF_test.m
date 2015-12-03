%% Reading files and loading image stack
cd('S:\Data\633\Widefield imaging\151110_633\633behavior_vistraining_CS_33exp_30Hz_1');
imageStack = '633behavior_vistraining_CS_33exp_30Hz_1_MMStack.ome';
data = double(readtiff([imageStack '.tif']));

mouse = '633';
date = '151110';
time = '1432';
cd('Z:\home\andrew\Behavior\Data');
load(['data-i' mouse '-' date '-' time]);

%% Initializing variables
frame_rate = 30;
pre_trial_time = 1000;
post_trial_time = 1000;
pre_trial_frames = ceil(pre_trial_time*(frame_rate./1000));
post_trial_frames = ceil(post_trial_time*(frame_rate./1000));
siz = size(data);
nTrials = length(input.cLeverDown);
trialStart = zeros(1,nTrials);
trialEnd = zeros(1,nTrials);
trial_dataF = zeros(siz(1),siz(2),nTrials);
trial_data = cell(1,nTrials);
trial_dataDF = cell(1,nTrials);
trial_dataDFoverF = cell(1,nTrials);

%% Calculating DFoverF for each frame
%Baseline F from 30 frames before each trial start
for i = 2:length(input.cLeverDown)
    trialStart(i) = input.cLeverDown{i} - 60;
    trialEnd(i) = input.cLeverUp{i} + 60;
    trial_data{i} = data(:,:,trialStart(i):trialEnd(i));
    trial_dataF(:,:,i) = mean(data(:,:,trialStart(i) + 10:trialStart(i) + 30),3);
end   
for i = 2:length(input.cLeverDown)
    trial_dataDF{i} = bsxfun(@minus, trial_data{i}, trial_dataF(:,:,i));
end
clear trial_data
for i = 2:length(input.cLeverDown)
    trial_dataDFoverF{i} = bsxfun(@rdivide, trial_dataDF{i}, trial_dataF(:,:,i));
end
clear trial_dataDF

%% Indexing success, false alarms and misses, hold times and react times
successIx = find(strcmp(input.trialOutcomeCell, 'success'));
earlyIx = find(strcmp(input.trialOutcomeCell, 'failure'));
missedIx = find(strcmp(input.trialOutcomeCell, 'ignore'));
holdTimes = cell2mat(input.holdTimesMs);
reactTimes = cell2mat(input.reactTimesMs);

%% Response around lever press
%60 frames before to 60 frames after
trialStart_win = 400;
ceil(trialStart_win*(frame_rate./1000));
trialStart_data = zeros(siz(1),siz(2),120,nTrials);
for i = 2:length(input.cLeverDown)
    trialStart_data(:,:,:,i) = trial_dataDFoverF{i}(:,:,1:120);
end
trialStart_data_avg = mean(trialStart_data,4);
%writetiff(trialStart_data_avg,'S:\Analysis\633\Widefield imaging\151108_633\leverPress2.tif');

%% Response around lever release in success trials
%60 frames before to 60 frames after
nSuccess = length(successIx);
sucTrialEnd_data = zeros(siz(1),siz(2),120,nSuccess);
sucTrial_dataDFoverF = cell(1,nSuccess);
for i = 2:nSuccess
    sucTrial_dataDFoverF{i} = trial_dataDFoverF{successIx(i)};
    trialL = size(sucTrial_dataDFoverF{i},3);
    sucTrialEnd_data(:,:,:,i) = sucTrial_dataDFoverF{i}(:,:,trialL-119:trialL);
end
sucTrialEnd_data_avg = mean(sucTrialEnd_data,4);
%writetiff(sucTrialEnd_data_avg,'S:\Analysis\633\Widefield imaging\151108_633\sucLeverRelease2.tif');

%% Response around lever release in early trials
%30 frames before to 30 frames after
nEarly = length(earlyIx);
earlyTrialEnd_data = zeros(siz(1),siz(2),120,nEarly);
earlyTrial_dataDFoverF = cell(1,nEarly);
for i = 2:nEarly
    earlyTrial_dataDFoverF{i} = trial_dataDFoverF{earlyIx(i)};
    trialL = size(earlyTrial_dataDFoverF{i},3);
    earlyTrialEnd_data(:,:,:,i) = earlyTrial_dataDFoverF{i}(:,:,trialL-119:trialL);
end
earlyTrialEnd_data_avg = mean(earlyTrialEnd_data,4);
%writetiff(earlyTrialEnd_data_avg,'S:\Analysis\633\Widefield imaging\151108_633\earlyLeverRelease2.tif');
    
%% Response around lever release in missed trials
%30 frames before to 30 frames after
nMissed = length(missedIx);
missedTrialEnd_data = zeros(siz(1),siz(2),120,nMissed);
missedTrial_dataDFoverF = cell(1,nMissed);
for i = 2:nMissed
    missedTrial_dataDFoverF{i} = trial_dataDFoverF{missedIx(i)};
    trialL = size(missedTrial_dataDFoverF{i},3);
    missedTrialEnd_data(:,:,:,i) = missedTrial_dataDFoverF{i}(:,:,trialL-119:trialL);
end
missedTrialEnd_data_avg = mean(missedTrialEnd_data,4);
%writetiff(missedTrialEnd_data_avg,'S:\Analysis\633\Widefield imaging\151108_633\missedLeverRelease2.tif');    
    
%% Response around orientation change in success trials using pre-target F
%around_target_frames = 30;

%for i = 1:nSuccess
%    sucTargetFrame(i) = input.cTargetOn(successIx(i));
%    sucTargetBef(i) = sucTargetFrame{i} - around_target_frames;
%    sucTargetAft(i) = sucTargetFrame{i} + around_target_frames;
%    sucTrialTarget_data{i} = data(:,:,sucTargetBef(i):sucTargetAft(i));
%    sucTrialTarget_dataF(:,:,i) = mean(data(:,:,sucTargetBef(i):sucTargetBef(i) + 20),3);
%end
%for i = 1:nSuccess
%    sucTrialTarget_dataDF{i} = bsxfun(@minus, sucTrialTarget_data{i}, sucTrialTarget_dataF(:,:,i));
%end
%clear sucTrialTarget_data
%for i = 1:nSuccess
%    sucTrialTarget_dataDFoverF{i} = bsxfun(@rdivide, sucTrialTarget_dataDF{i}, sucTrialTarget_dataF(:,:,i));
%end
%clear sucTrialTarget_dataDF
%for i = 1:nSuccess
%    sucTrialTarget_data(:,:,:,i) = sucTrialTarget_dataDFoverF{i}(:,:,around_target_frames - 29:around_target_frames + 30);
%end

%sucTrialTarget_data_avg = mean(sucTrialTarget_data,4);
%writetiff(sucTrialTarget_data_avg,'S:\Analysis\633\Widefield imaging\151105_633\sucTrialTarget_preTargetF.tif');

%% Response around orientation change in success trials using pre-trial F

around_target_frames = 75;
sucTargetFrame = cell(1,nSuccess);
sucTargetBef = zeros(1,nSuccess);
sucTargetAft = zeros(1,nSuccess);
sucTrialTarget_data = cell(1,nSuccess);
sucTrialTarget_dataDF = cell(1,nSuccess);
sucTrialTarget_dataDFoverF = cell(1,nSuccess);
sucTrialTarget_dataset = zeros(siz(1),siz(2),150,nSuccess);

for k = 2:nSuccess
    sucTargetFrame(k) = input.cTargetOn(successIx(k));
    sucTargetBef(k) = sucTargetFrame{k} - around_target_frames;
    sucTargetAft(k) = sucTargetFrame{k} + around_target_frames;
    sucTrialTarget_data{k} = data(:,:,sucTargetBef(k):sucTargetAft(k));
end
for k = 2:nSuccess
    sucTrialTarget_dataDF{k} = bsxfun(@minus, sucTrialTarget_data{k}, trial_dataF(:,:,k));
end
clear sucTrialTarget_data
for k = 2:nSuccess
    sucTrialTarget_dataDFoverF{k} = bsxfun(@rdivide, sucTrialTarget_dataDF{k}, trial_dataF(:,:,k));
end
clear sucTrialTarget_dataDF
for k = 2:nSuccess
    sucTrialTarget_dataset(:,:,:,k) = sucTrialTarget_dataDFoverF{k}(:,:,around_target_frames - 74:around_target_frames + 75);
end

sucTrialTarget_data_avg = mean(sucTrialTarget_dataset,4);
%writetiff(sucTrialTarget_data_avg,'S:\Analysis\633\Widefield imaging\151108_633\sucTrialTarget_preTrialF2.tif');

%% Response around orientation change in missed trials using pre-trial F

around_target_frames = 75;
missTargetFrame = cell(1,nMissed);
missTargetBef = zeros(1,nMissed);
missTargetAft = zeros(1,nMissed);
missTrialTarget_data = cell(1,nMissed);
missTrialTarget_dataDF = cell(1,nMissed);
missTrialTarget_dataDFoverF = cell(1,nMissed);
missTrialTarget_dataset = zeros(siz(1),siz(2),150,nMissed);

for k = 2:nMissed
    missTargetFrame(k) = input.cTargetOn(missedIx(k));
    missTargetBef(k) = missTargetFrame{k} - around_target_frames;
    missTargetAft(k) = missTargetFrame{k} + around_target_frames;
    missTrialTarget_data{k} = data(:,:,missTargetBef(k):missTargetAft(k));
end
for k = 2:nMissed
    missTrialTarget_dataDF{k} = bsxfun(@minus, missTrialTarget_data{k}, trial_dataF(:,:,k));
end
clear missTrialTarget_data
for k = 2:nMissed
    missTrialTarget_dataDFoverF{k} = bsxfun(@rdivide, missTrialTarget_dataDF{k}, trial_dataF(:,:,k));
end
clear missTrialTarget_dataDF
for k = 2:nMissed
    missTrialTarget_dataset(:,:,:,k) = missTrialTarget_dataDFoverF{k}(:,:,around_target_frames - 74:around_target_frames + 75);
end

missTrialTarget_data_avg = mean(missTrialTarget_dataset,4);
%writetiff(missTrialTarget_data_avg,'S:\Analysis\633\Widefield imaging\151108_633\missedTrialTarget_preTrialF2.tif');

%% Plot timecourses for events for single days based on ROIs

leverPressTC_AL = stackGetTimeCourses(trialStart_data_avg, mask_cell_AL);
leverPressTC_LM = stackGetTimeCourses(trialStart_data_avg, mask_cell_LM);
leverPressTC_RL = stackGetTimeCourses(trialStart_data_avg, mask_cell_RL);
leverPressTC_V = stackGetTimeCourses(trialStart_data_avg, mask_cell_V);
leverPressTC_CS = stackGetTimeCourses(trialStart_data_avg, mask_cell_CS);
leverPressTC_CV = stackGetTimeCourses(trialStart_data_avg, mask_cell_CV);

sucLeverReleaseTC_AL = stackGetTimeCourses(sucTrialEnd_data_avg, mask_cell_AL);
sucLeverReleaseTC_LM = stackGetTimeCourses(sucTrialEnd_data_avg, mask_cell_LM);
sucLeverReleaseTC_RL = stackGetTimeCourses(sucTrialEnd_data_avg, mask_cell_RL);
sucLeverReleaseTC_V = stackGetTimeCourses(sucTrialEnd_data_avg, mask_cell_V);
sucLeverReleaseTC_CS = stackGetTimeCourses(sucTrialEnd_data_avg, mask_cell_CS);
sucLeverReleaseTC_CV = stackGetTimeCourses(sucTrialEnd_data_avg, mask_cell_CV);

earlyLeverReleaseTC_AL = stackGetTimeCourses(earlyTrialEnd_data_avg, mask_cell_AL);
earlyLeverReleaseTC_LM = stackGetTimeCourses(earlyTrialEnd_data_avg, mask_cell_LM);
earlyLeverReleaseTC_RL = stackGetTimeCourses(earlyTrialEnd_data_avg, mask_cell_RL);
earlyLeverReleaseTC_V = stackGetTimeCourses(earlyTrialEnd_data_avg, mask_cell_V);
earlyLeverReleaseTC_CS = stackGetTimeCourses(earlyTrialEnd_data_avg, mask_cell_CS);
earlyLeverReleaseTC_CV = stackGetTimeCourses(earlyTrialEnd_data_avg, mask_cell_CV);

missedLeverReleaseTC_AL = stackGetTimeCourses(missedTrialEnd_data_avg, mask_cell_AL);
missedLeverReleaseTC_LM = stackGetTimeCourses(missedTrialEnd_data_avg, mask_cell_LM);
missedLeverReleaseTC_RL = stackGetTimeCourses(missedTrialEnd_data_avg, mask_cell_RL);
missedLeverReleaseTC_V = stackGetTimeCourses(missedTrialEnd_data_avg, mask_cell_V);
missedLeverReleaseTC_CS = stackGetTimeCourses(missedTrialEnd_data_avg, mask_cell_CS);
missedLeverReleaseTC_CV = stackGetTimeCourses(missedTrialEnd_data_avg, mask_cell_CV);

sucTrialTargetTC_AL = stackGetTimeCourses(sucTrialTarget_data_avg, mask_cell_AL);
sucTrialTargetTC_LM = stackGetTimeCourses(sucTrialTarget_data_avg, mask_cell_LM);
sucTrialTargetTC_RL = stackGetTimeCourses(sucTrialTarget_data_avg, mask_cell_RL);
sucTrialTargetTC_V = stackGetTimeCourses(sucTrialTarget_data_avg, mask_cell_V);
sucTrialTargetTC_CS = stackGetTimeCourses(sucTrialTarget_data_avg, mask_cell_CS);
sucTrialTargetTC_CV = stackGetTimeCourses(sucTrialTarget_data_avg, mask_cell_CV);

missTrialTargetTC_AL = stackGetTimeCourses(missTrialTarget_data_avg, mask_cell_AL);
missTrialTargetTC_LM = stackGetTimeCourses(missTrialTarget_data_avg, mask_cell_LM);
missTrialTargetTC_RL = stackGetTimeCourses(missTrialTarget_data_avg, mask_cell_RL);
missTrialTargetTC_V = stackGetTimeCourses(missTrialTarget_data_avg, mask_cell_V);
missTrialTargetTC_CS = stackGetTimeCourses(missTrialTarget_data_avg, mask_cell_CS);
missTrialTargetTC_CV = stackGetTimeCourses(missTrialTarget_data_avg, mask_cell_CV);