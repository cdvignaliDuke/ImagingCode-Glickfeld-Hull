%experiment info
mouse = 'AW58';
SubNum = '657';
expdate = '160902';
mw_time = '1450';
filename = 'pairedPulseAud_1.tif';
foldername = 'pairedPulseAud_1';
interval_ms = 250;

frameRate = 1000/interval_ms;

fnin = fullfile('Z:\data',mouse,'widefield imaging',[mouse '_' expdate]);
fnout = fullfile('Z:\Analysis',mouse,'widefield imaging',expdate);
mkdir(fnout);

%load data
data = readtiff(fullfile(fnin,foldername,filename));
data = double(data);

load(['Y:\home\andrew\Behavior\Data\data-i' SubNum '-' expdate '-' mw_time])

%important variables
AVtrials = find(cell2mat(input.tBlock2TrialNumber) == 0);
Atrials = find(cell2mat(input.tBlock2TrialNumber) == 1);
ntrials = input.trialSinceReset;
preTrialFrames = 15;
stimON = input.nFramesOn;
stimOFF = input.nFramesOff;
trialFrames = (stimON+stimOFF)*2;
trialStart_ind = cell2mat(input.cLeverDown);
pulse1_win = preTrialFrames+1:preTrialFrames+stimON;
pulse2_win = preTrialFrames+stimON+stimOFF:preTrialFrames+(2*stimON)+stimOFF;

%find trial-by-trial dF/F
dFoverF = zeros(size(data,1),size(data,2),preTrialFrames+trialFrames,ntrials);
for i = 1:ntrials
    F = nanmean(data(:,:,trialStart_ind(i)-preTrialFrames:trialStart_ind(i)),3);
    dFoverF(:,:,:,i) = bsxfun(@rdivide, bsxfun(@minus, data(:,:,trialStart_ind(i)-preTrialFrames:trialStart_ind(i)+trialFrames-1),F),F);
end
writetiff(squeeze(mean(dFoverF,3)),[fnout '\pairedpulse_dFoverF']);
  
%first pulse response
allTrials_avg_1p = nanmean(nanmean(dFoverF(:,:,pulse1_win,:),4),3);
AV_avg_1p = nanmean(nanmean(dFoverF(:,:,pulse1_win,AVtrials),4),3);
A_avg_1p = nanmean(nanmean(dFoverF(:,:,pulse1_win,Atrials),4),3);
writetiff(AV_avg_1p,[fnout '\pairedpulse_AV_avg_1p'])
writetiff(A_avg_1p,[fnout '\pairedpulse_A_avg_1p'])

%second pulse response
allTrials_avg_2p = nanmean(nanmean(dFoverF(:,:,pulse2_win,:),4),3);
AV_avg_2p = nanmean(nanmean(dFoverF(:,:,pulse2_win,AVtrials),4),3);
A_avg_2p = nanmean(nanmean(dFoverF(:,:,pulse2_win,Atrials),4),3);
writetiff(AV_avg_2p,[fnout '\pairedpulse_AV_avg_2p'])
writetiff(A_avg_2p,[fnout '\pairedpulse_A_avg_2p'])

%%
% plot images
figure;
suptitle('mean of response during 500 ms stimulus window')
subplot(2,3,1)
imagesc(allTrials_avg_1p)
caxis([-0.01 0.005])

    axis square
title('all trials')
ylabel('1st pulse')
subplot(2,3,2)
imagesc(AV_avg_1p)
caxis([-0.01 0.005])

    axis square
title('vis + aud')
subplot(2,3,3)
imagesc(A_avg_1p)
caxis([-0.01 0.005])

    axis square
title('aud only')
subplot(2,3,4)
imagesc(allTrials_avg_2p)
caxis([-0.01 0.005])

    axis square
ylabel('2nd pulse')
subplot(2,3,5)
imagesc(AV_avg_2p)
caxis([-0.01 0.005])

    axis square
subplot(2,3,6)
imagesc(A_avg_2p)
caxis([-0.01 0.005])
axis square
%%
%plot subtractions
figure;
suptitle('response minus no-stim control')
subplot(1,3,1)
imagesq(A_avg_1p-A_avg_2p);
caxis([-0.01 0.005])
colorbar
title('auditory only')
axis square
subplot(1,3,2)
imagesq(AV_avg_2p--A_avg_2p)
caxis([-0.0001 0.0001])
colorbar
title('visual only')
axis square
subplot(1,3,3)
imagesq(AV_avg_1p--A_avg_2p)
caxis([-0.0001 0.0001])
colorbar
title('visual+auditory only')
axis square

%
figure
subplot(1,2,1)
imagesc(AV_avg_1p-A_avg_1p)
caxis([-0.0001 0.0001])
colorbar
axis square
title('vis+aud minus aud only')
subplot(1,2,2)
imagesc(AV_avg_1p-AV_avg_2p)
caxis([-0.0001 0.0001])
colorbar
axis square
title('vis+aud minus vis only')


