%%This script is to extract frame index of different speed modes.
% during motorized treadmill experiment in 2P imaging
% runoff_fast_mat, runoff_slow_mat,runtrig_fast_mat, runtrig_slow_mat,
% run_off_mat, runtrig_mat are maxtrix containing
% frame index for running triggered average and runing offset average.
% each row is a trial 
% matrix for tone trigger during different behavioral states:
% tone_fast_mat,tone_slow_mat,tone_stay_mat,tone_run_mat
% lastly, stay_woTone includes all stationary frames excluding the first second 
% of stationary trial and 1s after tone delivery.

%%
clear;
folder = 'Z:\Data\behavior\motorizedWheel\tone\';
sessionID = '1064-200319';% if there're more than 1 sessions on a single subject on the same day, put the time in here after the date
filename = dir([folder 'data-i' '*' sessionID  '*' ]);
for i = 1: size(filename,1)
    behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\tone\behavioral_analysis\' sessionID '_' num2str(i)];
    if ~exist(behav_dest)
        mkdir(behav_dest);
    end
end

%% load data
file = [folder, filename(i).name];
load(file);
countTime = cell2mat(input.counterTimesUs);
countValue = input.counterValues;
nstart = find(countValue{1} == 1, 1, 'last');
countValue{1} = countValue{1}(nstart:end);

tonesOn = input.cTonesOn;
tonesOn{1} = tonesOn{1}(2:end); %there's always a zero at the beginning of the experiment, so delete this
tonesOn_vec = cell2mat(tonesOn);
tonesOn_vec = double(tonesOn_vec);
commandVolt = cell2mat(input.commandVoltage);

% if all running trials and stationary trials are the same length:
ntrials = length(countValue);
% when commnadVolt == 0, speed == 0 
% when commandVolt == 2, speed == 9.6-12cm/s
% when commandVolt == 3, speed == 14.4-18cm/s
% count how many trials there are for each condition 
nslowTrials = sum(commandVolt == 2);
nfastTrials = sum(commandVolt == 3);
nstayTrials = sum(commandVolt == 0);

%% matrix for running trig and running offset, speed change, and tone trigger
%--------------------------------------running-----------------------------------------------------------
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

%--------------------------------------------tone-----------------------------------------------------
t1 = 0; 
t2 = 0;
t3 = 0;
tonesOn_stay = cell(1,nstayTrials);
tonesOn_slow = cell(1,nslowTrials);
tonesOn_fast = cell(1,nfastTrials);

for t = 1: ntrials
    if     commandVolt(t) == 0
        t1 = t1 + 1;
        tonesOn_stay{t1} = tonesOn{t};
    elseif commandVolt(t) == 2
        t2 = t2 + 1;
        tonesOn_slow{t2} = tonesOn{t};
    elseif commandVolt(t) == 3
        t3 = t3 + 1;
        tonesOn_fast{t3} = tonesOn{t};
    end
end
tonesOn_stay = cell2mat(tonesOn_stay);
tonesOn_slow = cell2mat(tonesOn_slow);
tonesOn_fast = cell2mat(tonesOn_fast);

len2 = 60;
tone_slow_mat = zeros(2*nslowTrials,len2); % each row is an tone trial (each running trial has 2 tones)
tone_fast_mat = zeros(2*nfastTrials,len2);
tone_stay_mat = zeros(2*nstayTrials,len2);
for a = 1:length(tonesOn_slow)
    tone_slow_mat(a,:) =  tonesOn_slow(a)-29:tonesOn_slow(a)+30; % the 30th column is the frame when tone happens
end
for a = 1:length(tonesOn_fast)
    tone_fast_mat(a,:) =  tonesOn_fast(a)-29:tonesOn_fast(a)+30; % the 30th column is the frame when tone happens
end
for a = 1:length(tonesOn_stay)
    tone_stay_mat(a,:) =  tonesOn_stay(a)-29:tonesOn_stay(a)+30; % the 30th column is the frame when tone happens
end
tone_run_mat = cat(1,tone_slow_mat,tone_fast_mat);

save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
    'tonesOn','tonesOn_vec', 'runtrig_slow_mat','runtrig_fast_mat','runoff_slow_mat',...
    'runoff_fast_mat','spd_decrease_mat','spd_increase_mat','runtrig_mat',...
    'runoff_mat','tonesOn_stay','tonesOn_slow','tonesOn_fast',...
    'tone_slow_mat','tone_fast_mat','tone_stay_mat','tone_run_mat');

%% write a vector that has commandvoltages and is the same length as frames 
% need a for loop
speed_inVolt  = [];
for v = 1: length(commandVolt)
    speed_inVolt(countValue{v}) = commandVolt(v);
end

figure; plot(speed_inVolt);
ylim([-0.5 3.5]);ylabel('Volts');xlabel('frames');
vline(tonesOn_vec,'red');

savefig([behav_dest '\' sessionID '_' num2str(i) '_speed.fig']);

save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
    'speed_inVolt','-append');

%% frames during stationary without tone ---- useful for deconvolution, sort out bad cells
stay_all_frames = [];
for c = 1: length(commandVolt)
    if commandVolt(c) == 0
        stay_all_frames = cat(2,stay_all_frames,countValue{c}(31:end)); % delete the first second because the mice might be transitioning from the last state
    end
end
tone_trig_stay = [];
for a = 1: length(tonesOn_stay)
    tone_trig_stay = cat(2,tone_trig_stay,tonesOn_stay(a)-1:tonesOn_stay(a)+30);
end
stay_woTone = setdiff(stay_all_frames,tone_trig_stay);
stay_woTone = double(stay_woTone);

save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
    'stay_woTone','-append');




