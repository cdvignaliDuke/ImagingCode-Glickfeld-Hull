% Code by Tiago Marques, edited by rest of Cortical Circuits lab
% This script needs to be run once for each set of repetitions per
% stimulus. We are currently running 4 sets of stimuli: 
% (Vr vertical right, VL - vertical left, HD - horizontal down; HU-
% horizontal Up)

%% Experiment parameters 
clear all
close all
clc
%
animalName = 'i1303';
date = '191118';
time = '1042';
cond = {'VR', 'HD', 'VL', 'HU'};
stim_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
database_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\Widefield_images';
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Widefield_imaging';

fprintf(['Create ' animalName ' ' date ' relative maps \n' ])
%% stimulus data

stimdata_fn = ['data-' animalName '-' date '-' time '.mat'];
load(fullfile(stim_fn,stimdata_fn))

fCamera = double(input.frameRateHz);

fprintf(['Imported stim data \n' ])

%% Script variables
maxSize = 10; % Max size in Gb to load the videos - RAM dependent
save_flag = 1;

% FOR DEBUGGING - DON'T CHANGE
tbin=1;

%% Import section
% Paths for stimulus information and camera videos
stimNum = celleqel2mat_padded(input.tStimulusNumber);
nCond = length(unique(celleqel2mat_padded(input.tStimulusNumber)));
cITIStart = celleqel2mat_padded(input.cITIStart);
cStimOn = celleqel2mat_padded(input.cStimOn);
expt_fn = fullfile(database_fn,[date '_' animalName],[date '_' animalName '_movingBar_1'],[date '_' animalName '_movingBar_1_MMStack_Pos0.ome.tif']);
data = readtiff(expt_fn,'single');
vessels = mean(data(:,:,1:100),3);

fprintf(['Imported image data \n' ])

%% measure phase and amplitude for each condition
for iCond = 1:nCond
    
    fprintf(['Analyzing ' cond{iCond}  '\n' ])
    ind = find(stimNum == iCond-1);
    sz = size(data);
    framesPerBin = cITIStart(iCond+1)-cStimOn(iCond)+ fCamera + fCamera;
    stimPeriod = framesPerBin./fCamera;
    rep_bin = zeros(sz(1),sz(2),stimPeriod.*fCamera,length(ind));
    for it = 1:length(ind)
        if ind(it)+1 <= length(cITIStart)
            rep_bin(:,:,:,it) = data(:,:,cStimOn(ind(it))-fCamera:cITIStart(ind(it)+1)-1+fCamera);
        else
            rep_bin(:,:,:,it) = data(:,:,cStimOn(ind(it))-fCamera:input.counterValues{end}(end)-1+fCamera);
        end
    end
    frame = reshape(rep_bin,[sz(1) sz(2) framesPerBin.*length(ind)]);
    rep_bin = mean(rep_bin,4);
    framesTotal = size(frame,3);

% Determines number of image splits
    totalSize = framesTotal*size(frame,1)*size(frame,2)*8/10^9;
    nSplit = ceil(totalSize/maxSize);
    xSplit = floor(size(frame,1)/nSplit);

    fFourier = fCamera/tbin;

   % FFT parameters
    tFourier = 1/fFourier;
    fStim = 1/stimPeriod;
    L = framesTotal;
    NFFT = 2^nextpow2(L);
    t = (0:L-1)*tFourier;
    f = fFourier/2*linspace(0,1,(L/2+1));
    [fClo, fIndex] = min(abs(f-fStim));

    phase = zeros(size(frame,1),size(frame,2));
    phase_bef=zeros(size(frame,1),size(frame,2));
    phase_aft=zeros(size(frame,1),size(frame,2));

    amplitude = zeros(size(frame,1),size(frame,2));

    % For each image split
    for i =1:nSplit
        if i<nSplit % All except last, split size is the standard
            xSize = xSplit;
        else    % Last split, split size is the remaining part of the image
            xSize = size(frame,1)-(nSplit-1)*xSplit;
        end
        % Matrix to store the video for the current split
        pxTraces = zeros(xSize,size(frame,2), framesTotal);

        %%%%%%%%%% Load file section

        for k = 1:length(ind)
            for j=1:framesPerBin
                % Adds to the current video
                pxTraces(:,:,(k-1)*framesPerBin+j)= rep_bin((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:,j);
            end
        end


        %%%%%%%%%% Fourier Transform Section

        % Performs fast fourier transform for current split
        imageFourier = fft(pxTraces,[],3);
        imageFourier = imageFourier(:,:,1:L/2+1); % Trims to positive frequencies
        % Gets phase and amplitude of the fft
        phase((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:) = angle(imageFourier(:,:,fIndex));
        phase_bef((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:) = angle(imageFourier(:,:,fIndex-1));
        phase_aft((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:) = angle(imageFourier(:,:,fIndex+11));
        amplitude((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:) = 2*abs(imageFourier(:,:,fIndex));

    end

    % Figures outputs
    figure
    imagesc(fliplr(phase)); colormap(hsv)
    title(['Phase- ' cond{iCond}])

%     figure
%     imagesc(fliplr(phase_bef)); colormap(hsv)
%     title('Phase -df')
% 
%     figure
%     imagesc(fliplr(phase_aft)); colormap(hsv)
%     title('Phase +df')
% 
%     figure
%     imagesc(fliplr(amplitude)); colormap(jet)
% 
%     figure
%     imagesc(rot90(vessels,-1))
%     colormap('gray')
% 
%     figure
%     imagesc(fliplr(amplitude)); colormap(jet)



    % Saving section
    input_sub = trialChopper(input,ind);
    if ~exist(fullfile(fn_out, animalName, [date '_' animalName]))
        mkdir(fullfile(fn_out, animalName, [date '_' animalName]))
    end
    if save_flag
        save(fullfile(fn_out, animalName, [date '_' animalName], [date '_' animalName '_relative_' cond{iCond} '.mat']), 'phase','amplitude','vessels','input_sub');
    end
end