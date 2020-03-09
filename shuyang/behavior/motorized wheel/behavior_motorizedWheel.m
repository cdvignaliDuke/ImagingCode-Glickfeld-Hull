%%This script is to extract frame index of different speed modes.
% during motorized treadmill experiment in 2P imaging
% run_slow_mat, run_fast_mat, and stationary_mat are maxtrix containing
% frame index for slow, fast running trials and stationary trials.
% each row is a trial and the first second of each trial is being thrown
% away

% matrix for running triggered average and runing offset average:
%

%%
clear;
folder = 'Z:\Data\behavior\motorizedWheel\';
sessionID = '1049-200225';% if there're more than 1 sessions on a single subject on the same day, put the time in here after the date
filename = dir([folder 'data-i' '*' sessionID  '*' ]);
for i = 1: size(filename,1)
    behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\running\behavioral_analysis\' sessionID '_' num2str(i)];
    if ~exist(behav_dest)
        mkdir(behav_dest);
    end
end

%% matrix for each behavioral state
%load data
file = [folder, filename(i).name];
load(file);
countTime = cell2mat(input.counterTimesUs);
countValue = input.counterValues;
nstart = find(countValue{1} == 1, 1, 'last');
countValue{1} = countValue{1}(nstart:end);
% if all trials are the same length:
ntrials = length(countValue);
commandVolt = cell2mat(input.commandVoltage);
nframe_trial_run = length(countValue{find(commandVolt~=0,1)});%all running trials are the same length so just get index the first one to get a trial length
nframe_trial_stay = length(countValue{find(commandVolt==0,1)});
% when commnadVolt == 0, speed == 0 
% when commandVolt == 2, speed == 12cm/s
% when commandVolt == 3, speed == 18cm/s
% count how many trials there are for each condition 
nslowTrials = sum(commandVolt == 2);
nfastTrials = sum(commandVolt == 3);
nstayTrials = sum(commandVolt == 0);

% through out the first second(30 frames) of each trial
lenframe_include_run = nframe_trial_run - 30; % get rid of the first second
lenframe_include_stay = nframe_trial_stay - 30;
run_slow_mat = zeros(nslowTrials,lenframe_include_run); % each row is a trial
run_fast_mat = zeros(nfastTrials,lenframe_include_run);
staytionary_mat = zeros(nstayTrials,lenframe_include_stay); 

t1 = 0; % slow
t2 = 0; % fast
t3 = 0; % stationary
for t = 1:ntrials
    if     commandVolt(t) == 2
        t1 = t1 + 1;
        run_slow_mat(t1,:) = countValue{t}(31:end); % each row is a trial
    elseif commandVolt(t) == 3
        t2 = t2 + 1;
        run_fast_mat(t2,:) = countValue{t}(31:end);
    elseif commandVolt(t) == 0
        t3 = t3 + 1;
        staytionary_mat(t3,:) = countValue{t}(31:end);       
    end
end

% save matrices
save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
    'run_slow_mat','run_fast_mat','staytionary_mat');

%% matrix for running trig and running offset
% run trig slow
% run trig fast
% run off slow
% run off fast

t1 = 0; 
t2 = 0;
t3 = 0;
t4 = 0; 
t5 = 0; 
t6 = 0;

len = 75; 
runtrig_slow_mat = zeros(nslowTrials,len); % each row is a trial 
runoff_slow_mat = zeros(nslowTrials,len);
runtrig_fast_mat = zeros(nfastTrials,len);
runoff_fast_mat = zeros(nfastTrials,len);
spd_decrease_mat = zeros(nfastTrials,len);
spd_increase_mat = zeros(nfastTrials,len);
%run trig: stationary ---> running trial
%run off: running trial ---> stationary
%speed change: increase & decrease

for t = 1:ntrials-1
    if     commandVolt(t) == 0 && commandVolt(t+1) == 2 % runtrig slow
        t1 = t1 + 1;
        runtrig_slow_mat(t1,:) = countValue{t}(end)-29:countValue{t+1}(45); % each row is a trial
    elseif commandVolt(t) == 0 && commandVolt(t+1) == 3 % runtrig fast
        t2 = t2 + 1;
        runtrig_fast_mat(t2,:) = countValue{t}(end)-29:countValue{t+1}(45);
    elseif commandVolt(t) == 2 && commandVolt(t+1) == 0 % runoff slow
        t3 = t3 + 1;
        runoff_slow_mat(t3,:) = countValue{t}(end)-29:countValue{t+1}(45);    
    elseif commandVolt(t) == 3 && commandVolt(t+1) == 0 % runoff fast
        t4 = t4 + 1;
        runoff_fast_mat(t4,:) = countValue{t}(end)-29:countValue{t+1}(45);
    elseif commandVolt(t) == 3 && commandVolt(t+1) == 2 % speed decrease
        t5 = t5 + 1;
        spd_decrease_mat(t5,:) = countValue{t}(end)-29:countValue{t+1}(45);
    elseif commandVolt(t) == 2 && commandVolt(t+1) == 3 % speed increase
        t6 = t6 + 1;
        spd_increase_mat(t6,:) = countValue{t}(end)-29:countValue{t+1}(45);
     
    end
end
% delete rows with zeros because not all fast/slow speed trials will meet the criteria in the for loop
runtrig_slow_mat = runtrig_slow_mat(any(runtrig_slow_mat,2),:);
runtrig_fast_mat = runtrig_fast_mat(any(runtrig_fast_mat,2),:);
runoff_slow_mat = runoff_slow_mat(any(runoff_slow_mat,2),:);
runoff_fast_mat = runoff_fast_mat(any(runoff_fast_mat,2),:);
spd_decrease_mat = spd_decrease_mat(any(spd_decrease_mat,2),:);
spd_increase_mat = spd_increase_mat(any(spd_increase_mat,2),:);

runtrig_mat = cat(1,runtrig_slow_mat,runtrig_fast_mat);
runoff_mat = cat(1,runoff_slow_mat,runoff_fast_mat);

save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
    'runtrig_slow_mat','runtrig_fast_mat','runoff_slow_mat','runoff_fast_mat',...
    'spd_decrease_mat','spd_increase_mat','runtrig_mat','runoff_mat','-append');

%% write a vector that has commandvoltages and is the same length as frames 
% need a for loop
speed_inVolt  = [];
for v = 1: length(commandVolt)
    speed_inVolt(countValue{v}) = commandVolt(v);
end

figure; plot(speed_inVolt);
ylim([-0.5 3.5]);ylabel('Volts');

save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
    'speed_inVolt','-append');




