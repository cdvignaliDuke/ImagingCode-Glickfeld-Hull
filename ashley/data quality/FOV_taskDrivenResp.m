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

% %% which frames from this dataset are important? - create a frame index
% %   ****** for FSAV, the anticipation frames from successful trials that are at least 6
% %   cycles long - in V1, potentially other areas, the average response
% %   should be phasic, synced to the timing of each flashing grating
% CYC = 4;
% prePushFrames = 30;
% preTargetFrames = 5;
% postTargetFrames = input.nFramesOn+30;
% 
% frameRate = input.frameRateHz;
% cycles = find(cell2mat(input.tCyclesOn) == CYC);
% trialOutcome_S = find(strcmp(input.trialOutcomeCell,'success'));
% cLeverDown = cell2mat(input.cLeverDown);
% cTargetOn = cell2mat_padded(input.cTargetOn);
% cycTimeFrames = input.nFramesOn+input.nFramesOff;
% trCycLengthFrames = cycTimeFrames*CYC;
% trCycLengthMs = (trCycLengthFrames/frameRate)*1000;
% nTrials_SCYC = length(cLeverDown(intersect(trialOutcome_S,cycles)));
% frameInd_trMat_ant = linspaceNDim(double(cLeverDown(intersect(trialOutcome_S,cycles))-prePushFrames)', double(cLeverDown(intersect(trialOutcome_S,cycles))+trCycLengthFrames)',prePushFrames+trCycLengthFrames+1);
% frameInd_ant = reshape(frameInd_trMat_ant',1,size(frameInd_trMat_ant,1)*size(frameInd_trMat_ant,2));
% frameInd_trMat_tar = linspaceNDim(double(cTargetOn(intersect(trialOutcome_S,cycles))-preTargetFrames)', double(cTargetOn(intersect(trialOutcome_S,cycles))+postTargetFrames)',preTargetFrames+postTargetFrames+1);
% frameInd_tar = reshape(frameInd_trMat_tar',1,size(frameInd_trMat_tar,1)*size(frameInd_trMat_tar,2));

%% create a frame index for all trials from first 3 cycles, to see if cells are driven by the flashing stimulus
% first 4 cycles of all trials that are at least that long
% for target response, find all success trials, regardless of cycle length.
taskAz = num2str(input.gratingAzimuthDeg);
taskEl = num2str(input.gratingElevationDeg);
posStr = [taskAz '/' taskEl ' - Az/El'];
sizeStr = [num2str(input.gratingHeightDeg) 'x' num2str(input.gratingWidthDeg)];
    
CYC = 4:input.maxCyclesOn;
prePushFrames = 30;
preTargetFrames = 5;
postTargetFrames = input.nFramesOn+20;

frameRate = input.frameRateHz;
cycles = find(ismember(cell2mat(input.tCyclesOn),CYC));
trialOutcome_S = find(strcmp(input.trialOutcomeCell,'success'));
cLeverDown = cell2mat(input.cLeverDown);
cTargetOn = cell2mat_padded(input.cTargetOn);
cycTimeFrames = input.nFramesOn+input.nFramesOff;
trCycLengthFrames = cycTimeFrames*CYC(1);
trCycLengthMs = (trCycLengthFrames/frameRate)*1000;

nTrials_ant = length(cycles);
nTrials_tar = length(trialOutcome_S);

frameInd_trMat_ant = linspaceNDim(double(cLeverDown(cycles)-prePushFrames)', double(cLeverDown(cycles)+trCycLengthFrames)',prePushFrames+trCycLengthFrames+1);
frameInd_ant = reshape(frameInd_trMat_ant',1,size(frameInd_trMat_ant,1)*size(frameInd_trMat_ant,2));
frameInd_trMat_tar = linspaceNDim(double(cTargetOn(trialOutcome_S)-preTargetFrames)', double(cTargetOn(trialOutcome_S)+postTargetFrames)',preTargetFrames+postTargetFrames+1);
frameInd_tar = reshape(frameInd_trMat_tar',1,size(frameInd_trMat_tar,1)*size(frameInd_trMat_tar,2));
if any(frameInd_tar > nframes)
    frameInd_tar = frameInd_tar(1,1:size(frameInd_tar,2)-(preTargetFrames+postTargetFrames));
    nTrials_tar = nTrials_tar-1;
end

%% register those frames of interest, find timecourse, dF/F
% load registration img and mask from direction tuning
% fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
% cd(fileDirMasks);
% load('regImg.mat');
% load('mask&TCDir.mat');

% register img
data_avg = mean(data(:,:,dataAvgFrames),3);
figure; imagesc(data_avg);colormap gray

% beep 
% pause(60)

save(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'regImg'),'data_avg');

[out data_reg_ant] = stackRegister(data(:,:,frameInd_ant), data_avg);
[out data_reg_tar] = stackRegister(data(:,:,frameInd_tar), data_avg);
clear data

%max dF/F image for anticipation period frames
dF_reg = bsxfun(@minus,double(data_reg_ant),double(mean(data_reg_ant,3)));
dFoverF_reg = bsxfun(@rdivide,dF_reg,double(mean(data_reg_ant,3)));

maxdFoverF = max(dFoverF_reg,[],3);

IMG = maxdFoverF;
vals = unique(IMG(:));
maxCutoff_ind = ceil(0.9*length(vals));
maxCutoff = vals(maxCutoff_ind);
IMG(IMG < maxCutoff) = 0;
filterCutoff = 0.0015;
K = filter2(fspecial('average',5),IMG)/255;
Kmask_ant = K;
Kmask_ant(Kmask_ant < filterCutoff) = 0;
Kmask_ant(Kmask_ant > 0) = 1;
figure; 
subplot(1,4,1)
imagesc(maxdFoverF); colormap gray
title('anticipation maxdFoverF')
subplot(1,4,2)
imagesc(IMG); colormap gray
title('cutoff maxdFoverF')
subplot(1,4,3)
imagesc(K);colormap gray
title('smoothed cutoff')
subplot(1,4,4)
imagesc(Kmask_ant); colormap gray
title('mask')

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAVant_mask'), '-dpdf');

writetiff(maxdFoverF,fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAVant'));

%max dF/F image for anticipation period frames
dF_reg = bsxfun(@minus,double(data_reg_tar),double(mean(data_reg_tar,3)));
dFoverF_reg = bsxfun(@rdivide,dF_reg,double(mean(data_reg_tar,3)));

maxdFoverF = max(dFoverF_reg,[],3);

IMG = maxdFoverF;
vals = unique(IMG(:));
maxCutoff_ind = ceil(0.9*length(vals));
maxCutoff = vals(maxCutoff_ind);
IMG(IMG < maxCutoff) = 0;
filterCutoff = 0.0015;
K = filter2(fspecial('average',5),IMG)/255;
Kmask_tar = K;
Kmask_tar(Kmask_tar < filterCutoff) = 0;
Kmask_tar(Kmask_tar > 0) = 1;
figure; 
subplot(1,4,1)
imagesc(maxdFoverF); colormap gray
title('target maxdFoverF')
subplot(1,4,2)
imagesc(IMG); colormap gray
title('cutoff maxdFoverF')
subplot(1,4,3)
imagesc(K);colormap gray
title('smoothed cutoff')
subplot(1,4,4)
imagesc(Kmask_tar); colormap gray
title('mask')

print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAVtar_mask'), '-dpdf');

writetiff(maxdFoverF,fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAVtar'));

save(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'taskDrivenPixelsMask_targetandanticipation'), 'Kmask_tar','Kmask_ant')

%% anticipation period response
% get timecourses
% dataTC_ant = stackGetTimeCourses(data_reg_ant, mask_cell);
dataTC_ant = stackGetTimeCourses(data_reg_ant, Kmask_ant);
dataTC_trMat_ant = reshape(dataTC_ant,prePushFrames+trCycLengthFrames+1,nTrials_ant,size(dataTC_ant,2));

F_ant = mean(dataTC_trMat_ant(1:prePushFrames,:,:),1);
dF_ant = bsxfun(@minus, dataTC_trMat_ant,F_ant);
dFoverF_ant = bsxfun(@rdivide,dF_ant,F_ant);
dFoverF_trMean_ant = squeeze(mean(dFoverF_ant,2));

% %%
% baseRespCells = find(any(dFoverF_trMean_ant(1:prePushFrames,:) > taskRespCutoff,1));
% stimRespCells = find(any(dFoverF_trMean_ant(prePushFrames:size(dFoverF_trMean_ant,1),:) > taskRespCutoff,1));
% drivenCellsInd_ant = setdiff(stimRespCells,baseRespCells);
% 
% %% plot average timecourse
% figure;
% plot(mean(dFoverF_trMean_ant(:,drivenCellsInd_ant),2),'LineWidth',3,'Color','k')
% hold on
% vline(prePushFrames,'k')
% for i = 1:CYC(1)-1
%    vline((cycTimeFrames*i)+prePushFrames,'k:')
%    hold on
% end
% xlabel('frames')
% ylabel('dF/F')
% title({'task driven resp - driven cells'; [mouse '-' date]; [num2str(trCycLengthMs) 'ms/tr; respCutoff=' num2str(taskRespCutoff)]; ['nCells = ' num2str(length(drivenCellsInd_ant)) '/' num2str(size(dataTC_ant,2))]; [num2str(nTrials_ant) ' trials']})
% %%
% set(0,'defaultfigurepaperorientation','portrait');
% set(0,'defaultfigurepapersize',[8.5 11]);
% set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
% 
% print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAVant'), '-dpdf');
%%
figure;
plot(mean(dFoverF_trMean_ant,2),'LineWidth',3)
hold on
vline(prePushFrames,'k')
for i = 1:CYC(1)-1
   vline((cycTimeFrames*i)+prePushFrames,'k:')
   hold on
end
xlabel('frames')
ylabel('dF/F')
title({'task driven resp'; [mouse '-' date]; [num2str(trCycLengthMs) 'ms/tr;']; [num2str(nTrials_ant) ' trials'];[posStr '; ' sizeStr]})
%%
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAVant_all'), '-dpdf');






%% ********* target response
% get timecourses
% dataTC_tar = stackGetTimeCourses(data_reg_tar, mask_cell);
dataTC_tar = stackGetTimeCourses(data_reg_tar, Kmask_tar);

dataTC_trMat_tar= reshape(dataTC_tar,preTargetFrames+postTargetFrames+1,nTrials_tar,size(dataTC_tar,2));

F_tar = mean(dataTC_trMat_tar(1:preTargetFrames,:,:),1);
dF_tar = bsxfun(@minus, dataTC_trMat_tar,F_tar);
dFoverF_tar = bsxfun(@rdivide,dF_tar,F_tar);
dFoverF_trMean_tar = squeeze(mean(dFoverF_tar,2));

%%
% baseRespCells = find(any(dFoverF_trMean_tar(1:preTargetFrames,:) > taskRespCutoff,1));
% % stimRespCells = find(any(dFoverF_trMean_tar(preTargetFrames:size(dFoverF_trMean_tar,1),:) > taskRespCutoff,1));
% stimRespCells = find(any(dFoverF_trMean_tar(preTargetFrames+4:preTargetFrames+8,:) > taskRespCutoff,1));
% drivenCellsInd_tar = setdiff(stimRespCells,baseRespCells);

% %% plot average timecourse
% figure;
% plot(mean(dFoverF_trMean_tar(:,drivenCellsInd_tar),2),'LineWidth',3,'Color','k')
% hold on
% vline(preTargetFrames,'k')
% hold on
% yminmax = [min(mean(dFoverF_trMean_tar(:,drivenCellsInd_tar),2)) max(mean(dFoverF_trMean_tar(:,drivenCellsInd_tar),2))];
% resppatch = patch([preTargetFrames+4 preTargetFrames+8 preTargetFrames+8 preTargetFrames+4],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
% set(resppatch,'FaceAlpha',0.15);
% set(resppatch,'EdgeColor','none');
% xlabel('frames')
% ylabel('dF/F')
% title({'target driven resp - driven cells'; [mouse '-' date]; ['respCutoff=' num2str(taskRespCutoff)]; ['nCells = ' num2str(length(drivenCellsInd_tar)) '/' num2str(size(dataTC_tar,2))]; [num2str(nTrials_tar) ' success trials']})
% %%
% set(0,'defaultfigurepaperorientation','portrait');
% set(0,'defaultfigurepapersize',[8.5 11]);
% set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
% 
% print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAVtar'), '-dpdf');
%%
figure;
plot(mean(dFoverF_trMean_tar,2),'LineWidth',3)
hold on
vline(preTargetFrames,'k')
xlabel('frames')
ylabel('dF/F')
title({'target driven resp'; [mouse '-' date]; [num2str(trCycLengthMs) 'ms/tr;']; [num2str(nTrials_tar) ' success trials'];[posStr '; ' sizeStr]})
%%
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'cellsDrivenByFSAVtar_all'), '-dpdf');
%%

toc