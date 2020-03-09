%

clear;
folder = 'Z:\Data\behavior\motorizedWheel\';
sessionID = '1041-200116';% if there're more than 1 sessions on a single subject on the same day, put the time in here after the date
filename = dir([folder 'data-i' '*' sessionID  '*' ]);
for i = 1: size(filename,1)
    behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\behavioral_analysis\' sessionID '_' num2str(i)];
    if ~exist(behav_dest)
        mkdir(behav_dest);
    end
end

%%
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

lenframe_include_run = nframe_trial_run - 33; % get rid of the first second
lenframe_include_stay = nframe_trial_stay - 33;
run_slow = zeros(nslowTrials,lenframe_include_run); % each row is a trial
run_fast = zeros(nfastTrials,lenframe_include_run);
staytionary = zeros(nstayTrials,lenframe_include_stay); 

% have to do a for loop if don't use cellfun
t1 = 0; % slow
t2 = 0; % fast
t3 = 0; % stationary
for t = 1:ntrials
    if     commandVolt(t) == 2
        t1 = t1 + 1;
        run_slow(t1,:) = countValue{t}(34:end); % each row is a trial
    elseif commandVolt(t) == 3
        t2 = t2 + 1;
        run_fast(t2,:) = countValue{t}(34:end);
    elseif commandVolt(t) == 0
        t3 = t3 + 1;
        staytionary(t3,:) = countValue{t}(34:end);       
    end
end

% save matrices
