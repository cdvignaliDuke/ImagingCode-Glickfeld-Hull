%% load MWorks file

mworks = 'data-i023-140317-1616';
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
Dir = [0 45 90 135 180 225 270 315];
date = '140317';
mouse = 'G023';
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
trialNumbers = 1:length(trialon_avgspeed_mat);
trialon_run = (trialon_avgspeed_mat>4);
runTrials = find(trialon_run);
trialon_LocMatrix = [trialNumbers', trialon_avgspeed_mat',trialon_run'];
trialon_norun = (trialon_avgspeed_mat<4);
norunTrials = find(trialon_norun);
trialon_whichrun = trialon_avgspeed_mat(trialon_run);

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

%% add on indices to locomotion matrix
for itrial = 1:(nStim*nRep)
    nON_ind_firsts(itrial) = nON_ind(1+(nON*(itrial-1)));
end
trialon_LocMatrix = [trialNumbers', trialon_avgspeed_mat',trialon_run',nON_ind_firsts'];

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

%plot data
for iCell = 3;
    figure;
    for iStim = 1:nStim
        subplot(2,3,iStim)
        rep_mat = zeros(nON+nOFF,nRep);
        for iRep = 1:nRep
            plot(1-nOFF:nON,data_TC(squeeze(stim_mat(iStim,iRep,:))',iCell), 'k');
            hold on
            rep_mat(:,iRep) = data_TC(squeeze(stim_mat(iStim,iRep,:))',iCell);
        end
        plot(1-nOFF:nON, mean(rep_mat,2), 'r');
        hold on
        ylim([min(data_TC(:,iCell),[],1) max(data_TC(:,iCell),[],1)])
        stim_Az = rem(iStim,size(Az,2));
        if stim_Az == 0
            stim_Az = size(Az,2);
        end
        stim_El= ceil(iStim./size(Az,2));
        title(['Az = ' num2str(Az(1,stim_Az)) ' El = ' num2str(El(1,stim_El))]);
    end
end
%% find start off and on frames for locomotion pulses
data_TC
sitm_mat
trialon_avgspeed_mat
stim_mat_locperrep
stim_mat_locperon
%find locomotion for each stim type and each repeat
stim_mat_locperrep = zeros(nStim,1,nRep);
for iRep = 1:nRep
    for iStim = 1:nStim   
        stim_mat_locperrep(iStim,:,iRep) = trialon_avgspeed_mat(iStim + (nStim*(iRep-1)));
    end
end

stim_mat_locperon = zeros((nStim*nRep),nON);

%1+((iStim-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nStim)): nON + nOFF + ((iStim-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nStim));
