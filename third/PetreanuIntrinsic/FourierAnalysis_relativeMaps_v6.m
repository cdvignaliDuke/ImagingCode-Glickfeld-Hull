% Code by Tiago Marques, edited by rest of Cortical Circuits lab
% This script needs to be run once for each set of repetitions per
% stimulus. We are currently running 4 sets of stimuli: 
% (Vr vertical right, VL - vertical left, HD - horizontal down; HU-
% horizontal Up)

set(0,'defaultfigureposition',[100,200,480,360])

%% Experiment parameters 
saveName = 'RD10161_vr';
fCamera = 5;

%% Script variables
maxSize = 10; % Max size in Gb to load the videos - RAM dependent
save_flag = 1;

%% FOR DEBUGGING - DON'T CHANGE
tbin=1;

%% Import section
% Paths for stimulus information and camera videos
[expNumStim pathStimulus]=uigetfile('.mat','Select Experiment Stimulus File');
cd(pathStimulus);
cd ..

[fileVessels pathVessels]=uigetfile('.fig','Select Vessels Image'); % picture of skull / craniotomy with green LED
cd(pathStimulus);

cd ..
pathVideos=[cd filesep 'CAMERA'];

% Imports stimulus parameters
stimData = load(fullfile(pathStimulus,expNumStim));
vessels = getimage(openfig([pathVessels filesep fileVessels],'reuse','invisible'));

if isfield(stimData,'stimPeriod')==0
    stimData.stimPeriod = stimData.stimPeriod;
end

% Gets videos to analyze
d = dir(pathVideos);
str = {d.name};
str = sortrows({d.name}');
[s,v] = listdlg('PromptString','Select files to analyze:', 'OKString', 'OK',...
    'SelectionMode','multiple',...
    'ListString', str, 'Name', 'Select a File');
names = str(s);
numFiles = size(names, 1);

frame = read_qcamraw([pathVideos filesep names{1}],1);

%% Determines frames per file
framesPerFile = stimData.stimPeriod*fCamera;
framesPerBin = framesPerFile/tbin;
framesTotal = framesPerBin*numFiles;

%% Determines number of image splits
totalSize = framesTotal*size(frame,1)*size(frame,2)*8/1024/1024/1024;
nSplit = ceil(totalSize/maxSize);
xSplit = floor(size(frame,1)/nSplit);

fFourier = fCamera/tbin;

%% FFT parameters
tFourier = 1/fFourier;
fStim = 1/stimData.stimPeriod;
L = framesPerBin * numFiles;
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
    pxTraces = zeros(xSize,size(frame,2), framesPerBin * numFiles);
    %%
    %%%%%%%%%% Load file section
    %%
    for k = 1:(numFiles)
        names{k}
        % Loads the file
        rep = read_qcamraw([pathVideos filesep names{k}],1:framesPerFile);
        % Matrix to store in different compartments the different temporal phases
        repBin = zeros(size(frame,1),size(frame,2),framesPerBin,tbin);
        for p=1:tbin % Each temporal phase
            % Gets the 
            repBin(:,:,:,p) = rep(:,:,p:tbin:framesPerFile);
        end
        % Averages the temporal phases across the compartments to temporally bin the video
        repBin = mean(repBin,4);
        for j=1:framesPerBin
            % Adds to the current video
            pxTraces(:,:,(k-1)*framesPerBin+j)= repBin((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:,j);
        end
    end
    
    %%
    %%%%%%%%%% Fourier Transform Section
    %%
    %% Performs fast fourier transform for current split
    imageFourier = fft(pxTraces,[],3);
    imageFourier = imageFourier(:,:,1:L/2+1); % Trims to positive frequencies
    % Gets phase and amplitude of the fft
    phase((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:) = angle(imageFourier(:,:,fIndex));
    phase_bef((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:) = angle(imageFourier(:,:,fIndex-1));
    phase_aft((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:) = angle(imageFourier(:,:,fIndex+11));
    amplitude((xSplit*(i-1)+1):(xSplit*(i-1)+xSize),:) = 2*abs(imageFourier(:,:,fIndex));
    
end

%% Figures outputs
figure
imagesc(fliplr(phase)); colormap(hsv)
title('Phase')

figure
imagesc(fliplr(phase_bef)); colormap(hsv)
title('Phase -df')

figure
imagesc(fliplr(phase_aft)); colormap(hsv)
title('Phase +df')

figure
imagesc(fliplr(amplitude)); colormap(jet)

figure
imagesc(rot90(vessels,-1))
colormap('gray')

figure
imagesc(fliplr(amplitude)); colormap(jet)



%% Saving section
if save_flag
    save([cd filesep 'ANALYSIS' filesep saveName], 'phase','amplitude','vessels','stimData');
end
