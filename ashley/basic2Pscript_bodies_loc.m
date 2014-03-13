%% load MWorks file

mworks = 'data-i004-140306-1743';
load (mworks);

%% Parameters
frame_rate = input.frameImagingRateMs;
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = 150./down;
nOFF = 150./down;
nStim = 8;
Az = [-30];
El = [15];
Dirs = [0 45 90 135 180 225 270 315];
date = '140302';
mouse = 'G004';
%% Sum-up locomotion data for each trial and calculate speed (Stim ON only)

%create a matrix including all counters
nCounter = floor(((input.nScansOn./frame_rate).*1000)./double(input.speedTimerIntervalMs));
trial_loc_mat = zeros(nCounter,input.trialSinceReset);

for icount = 1:nCounter
    trial_loc_mat(icount, :) = cell2mat_padded(eval(['input.spCounter' num2str(icount)]));
end;

%sum locomotion for each trial (ON stimuli only)
trial_totalloc = sum(trial_loc_mat);

% Calculate speed for each trial in pulses/second
trialon_timeS = input.nScansOn./frame_rate;
trialon_avgspeed_mat = trial_totalloc./trialon_timeS;

%% Set threshold for locomotion and sort trials as running and not running
trialon_run = find(trialon_avgspeed_mat>0.5);
trialon_stat = find(trialon_avgspeed_mat<0.5);
trialon_whichrun = zeros(size(trialon_avgspeed_mat));
trialon_whichrun(trialon_run) = 1;

%% create indicies for stimulus type
nRep = 5;
stim_mat = zeros(nStim,nRep,(nON+nOFF));
start = 1;
for iRep = 1:nRep
    for iStim = 1:nStim   
        stim_mat(iStim,iRep,:) = 1+((iStim-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nStim)): nON + nOFF + ((iStim-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nStim));
    end
    start= start+nON+nOFF;
end

%% load 2P imaging data
data = readrawfile;

%% reshape data
%average signals in time

data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

%% register
data_avg = mean(data_sub(:,:,300:310),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%% create dF/F stack
nRep = size(data_reg,3)./((nON+nOFF)*nStim);

%find off and on frames
nOFF_ind = zeros(1,(nOFF*nStim*nRep));
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
nON_avg = mean(data_reg(:,:,nON_ind),3);
nOFF_avg = mean(data_reg(:,:,nOFF_ind),3);

%dF/F
dF_data = bsxfun(@minus,data_reg, nOFF_avg);
dFoverF_data = bsxfun(@rdivide, dF_data, nOFF_avg);
max_dF = max(dFoverF_data,[],3);
figure; imagesq(max_dF); colormap(gray)

%% use max dF/F to find ROIS

bwout = imCellEditInteractive(max_dF);
mask_cell = bwlabel(bwout);

%timecourses
data_TC = stackGetTimeCourses(dFoverF_data,mask_cell);
figure; tcOffsetPlot(data_TC)

%% analyze by stimulus type
%find indices for each stim type
stim_mat = zeros(nStim,nRep,(nON+nOFF));
start = 1;
for iRep = 1:nRep
    for iStim = 1:nStim   
        stim_mat(iStim,iRep,:) = 1+((iStim-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nStim)): nON + nOFF + ((iStim-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nStim));
    end
    start= start+nON+nOFF;
end

%sort by locomotion/stationary
