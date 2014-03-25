%% load MWorks files

mworks = 'data-i023-140317-1616';
load (mworks);

mworks2 = 'data-i023-140317-1640';
load(mworks2);

%% Parameters that may have changed
trial_Dir_2 = cell2mat(input.tGratingDirectionDeg);
trial_Dir_All = cat(2,trial_Dir,trial_Dir_2);


%% LOCOMOTION SECTION
%% Sum-up locomotion data for each trial and calculate speed (Stim ON only)

%create a matrix including all counters
nCounter = floor(((input.nScansOn./frame_rate).*1000)./double(input.speedTimerIntervalMs));
trial_loc_mat2 = zeros(nCounter,input.trialSinceReset);

for icount = 1:nCounter
    trial_loc_mat2(icount, :) = cell2mat_padded(eval(['input.spCounter' num2str(icount)]));
end;

%sum locomotion for each trial (ON stimuli only)
trial_totalloc = sum(trial_loc_mat2);

% Calculate speed for each trial in pulses/second
trialon2_timeS = input.nScansOn./frame_rate;
trialon_avgspeed_mat2 = trial_totalloc./trialon2_timeS;

trialon_avgspeed_matAll = cat(2,trialon_avgspeed_mat,trialon_avgspeed_mat2);

%% Set threshold for locomotion and sort trials as running and not running
trialNumbers = 1:length(trialon_avgspeed_matAll);
trialon_run = (trialon_avgspeed_matAll>4);
runTrials = find(trialon_run);
trialon_LocMatrix = [trialNumbers', trialon_avgspeed_matAll',trialon_run'];
trialon_norun = (trialon_avgspeed_matAll<4);
norunTrials = find(trialon_norun);
trialon_whichrun = trialon_avgspeed_matAll(trialon_run); 

%Go back to other script

%% 2P SECTION
%% load second 2P data sets

data = readrawfile;
data2 = readrawfile;

%% reshape data
%average signals in time

data2_down = stackGroupProject(data2,down);
clear data2

%remove negative data (by addition)
data2_sub = data2_down-min(min(min(data2_down,[],1),[],2),[],3);
clear data2_down

%% register
data2_avg = mean(data2_sub(:,:,300:310),3);
figure; imagesq(data2_avg); colormap(gray)

[out data2_reg] = stackRegister(data2_sub, data2_avg);
clear data2_sub

%% concatinate data_reg and data2_reg
data_regAll = cat(3,data_reg,data2_reg);
clear data_reg
clear data2_reg

nRep = size(data_regAll,3)./((nON+nOFF)*nStim);

%% create dF/F stack

%find off and on frames
nOFF_ind = zeros(1,(nOFF*nStim*nRep));
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

nON_ind = setdiff(1:size(data_regAll,3),nOFF_ind);
nON_avg = mean(data_regAll(:,:,nON_ind),3);
nOFF_avg = mean(data_regAll(:,:,nOFF_ind),3);

% find on indices for the first frame of each stimulus start period and iti (Off) period
for itrial = 1:(nStim*nRep)
    nON_ind_firsts(itrial) = nON_ind(1+(nON*(itrial-1)));
end
for itrial = 1:(nStim*nRep)
    nOFF_ind_firsts(itrial) = nOFF_ind(1+(nOFF*(itrial-1)));
end

%average nON and nOFF for specific trials
siz = size(data_regAll);
trial_nOFF_avg = zeros(siz(1), siz(2) ,size(nOFF_ind_firsts,2));
trial_nON_avg = zeros(siz(1), siz(2) ,size(nON_ind_firsts,2));
for itrial = 1:(nStim*nRep)
    trial_nOFF_avg(:,:,itrial) = mean(data_regAll(:,:,(nOFF_ind_firsts(itrial): (nOFF_ind_firsts(itrial)+nOFF)-1)),3);
end
for itrial = 1:(nStim*nRep)
    trial_nON_avg(:,:,itrial) = mean(data_regAll(:,:,(nON_ind_firsts(itrial): (nON_ind_firsts(itrial)+nON)-1)),3);
end

%dF/F
dF_data = bsxfun(@minus,data_regAll, nOFF_avg);
dFoverF_data = bsxfun(@rdivide, dF_data, nOFF_avg);
max_dF = max(dFoverF_data,[],3);
figure; imagesq(max_dF); colormap(gray)

%index stim types
stim_ind = zeros(nStim,nRep);
for idir = 1:nStim
    stim_ind(idir,:) = find(trial_Dir_All==Dir(idir));
end

%%Go back to other script.
