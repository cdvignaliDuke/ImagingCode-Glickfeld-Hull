tic
%% load 1 behavior dataset
% SubNum = '626';
% date = '160203';
% taskTime = '1053';
% taskFolder = '001';
% mouse = 'AW26';
% taskFName = '001_000_000';
% dirFolder = '005';

% load MWorks file
CD = 'Y:\home\andrew\Behavior\Data';
mworks = ['data-' 'i' SubNum '-' date '-' taskTime]; 
load (fullfile(CD,mworks));

% Set current directory to imaging data location
CD = ['Z:\data\' mouse '\two-photon imaging\' date '\' taskFolder];
cd(CD);
imgMatFile = [taskFName '.mat'];
load(imgMatFile);

nframes = info.config.frames;
tic
data = sbxread(taskFName,0,nframes);
toc
% pmt = 1; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
data = squeeze(data);

%% which frames from this dataset are important? - create a frame index
%   ****** for FSAV, the anticipation frames from successful trials that are at least 6
%   cycles long - in V1, potentially other areas, the average response
%   should be phasic, synced to the timing of each flashing grating
CYC = 6;
prePushFrames = 30;

frameRate = input.frameRateHz;
cycles = find(cell2mat(input.tCyclesOn) == CYC);
trialOutcome_S = find(strcmp(input.trialOutcomeCell,'success'));
cLeverDown = cell2mat(input.cLeverDown);
cycTimeFrames = input.nFramesOn+input.nFramesOff;
trCycLengthFrames = cycTimeFrames*CYC;
trCycLengthMs = (trCycLengthFrames/frameRate)*1000;
nTrials_SCYC = length(cLeverDown(intersect(trialOutcome_S,cycles)));
frameInd_trMat = linspaceNDim(double(cLeverDown(intersect(trialOutcome_S,cycles))-prePushFrames)', double(cLeverDown(intersect(trialOutcome_S,cycles))+trCycLengthFrames)',prePushFrames+trCycLengthFrames+1);
frameInd = reshape(frameInd_trMat',1,size(frameInd_trMat,1)*size(frameInd_trMat,2));
%% register those frames of interest, find timecourse, dF/F
% load registration img and mask from direction tuning
fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
cd(fileDirMasks);
load('regImg.mat');
load('mask&TCDir.mat');

% register img
[out data_reg] = stackRegister(data(:,:,frameInd), data_avg);
clear data

% dF = bsxfun(@minus,double(data_reg),double(mean(data_reg,3)));
% dFoverF = bsxfun(@rdivide,dF,double(mean(data_reg,3)));

maxF = max(data_reg,[],3);
% maxdFoverF = max(dFoverF,[],3);
figure;
imagesc(maxF); colormap gray
% figure;
% imagesc(maxdFoverF); colormap gray
writetiff(maxF,fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAV'));


% get timecourses
dataTC = stackGetTimeCourses(data_reg, mask_cell);

dataTC_trMat = reshape(dataTC,prePushFrames+trCycLengthFrames+1,nTrials_SCYC,size(dataTC,2));

F = mean(dataTC_trMat(1:prePushFrames,:,:),1);
dF = bsxfun(@minus, dataTC_trMat,F);
dFoverF = bsxfun(@rdivide,dF,F);
dFoverF_trMean = squeeze(mean(dFoverF,2));

%%
baseRespCells = find(any(dFoverF_trMean(1:prePushFrames,:) > taskRespCutoff,1));
stimRespCells = find(any(dFoverF_trMean(prePushFrames:size(dFoverF_trMean,1),:) > taskRespCutoff,1));
drivenCellsInd = setdiff(stimRespCells,baseRespCells);

%% plot average timecourse
figure;
plot(mean(dFoverF_trMean(:,drivenCellsInd),2),'LineWidth',3,'Color','k')
hold on
vline(prePushFrames,'k')
for i = 1:CYC-1
   vline((cycTimeFrames*i)+prePushFrames,'k:')
   hold on
end
xlabel('frames')
ylabel('dF/F')
title({'task driven resp - driven cells'; [mouse '-' date]; [num2str(trCycLengthMs) 'ms/tr; respCutoff=' num2str(taskRespCutoff)]; ['nCells = ' num2str(length(drivenCellsInd)) '/' num2str(size(dataTC,2))]; [num2str(nTrials_SCYC) '/' num2str(length(trialOutcome_S)) ' success trials']})
%%
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAV'), '-dpdf');
%%
figure;
plot(mean(dFoverF_trMean,2),'LineWidth',3)
hold on
vline(prePushFrames,'k')
for i = 1:CYC-1
   vline((cycTimeFrames*i)+prePushFrames,'k:')
   hold on
end
xlabel('frames')
ylabel('dF/F')
title({'task driven resp - all cells'; [mouse '-' date]; [num2str(trCycLengthMs) 'ms/tr; respCutoff=' num2str(taskRespCutoff)]; ['nCells = ' num2str(length(drivenCellsInd)) '/' num2str(size(dataTC,2))]; [num2str(nTrials_SCYC) '/' num2str(length(trialOutcome_S)) ' success trials']})
%%
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAV_all'), '-dpdf');
%%
toc