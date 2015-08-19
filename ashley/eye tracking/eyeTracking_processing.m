SubNum = '613';
date = '150810';
ImgFolder = '003';
mouse = 'AW13';
fName = '003_000_000';
frame_rate = 15;
calib = 1/26.6; %mm per pixel

fnout = ['Z:\home\lindsey\Analysis\Behavior\EyeTracking\' mouse '-' date '-' ImgFolder];
% % load MWorks file
% CD = ['Z:\data\' mouse '\mworks\' date];
% cd(CD);
% mworks = ['data-' 'i' SubNum '-' date '-' time]; 
% load (mworks);

% Set current directory to crash folder
CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\eye tracking\' date '\' ImgFolder];
cd(CD);

%% 
fn = [fName '_eye.mat'];
load(fn);          % should be a '*_eye.mat' file
    
data = squeeze(data);      % the raw images...
xc = size(data,2)/2;       % image center
yc = size(data,1)/2;
W=40;

rad_range = [15 25];
data = data(yc-W:yc+W,xc-W:xc+W,:);
warning off;

A = cell(size(data,3),1);
B = cell(size(data,3),1);
for n = 1:size(data,3)
    A{n} = [0,0];
    B{n} = [0];
end
eye = struct('Centroid',A,'Area',B);

for n = 1:size(data,3)
    [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.9);
    if(isempty(center))
        eye(n).Centroid = [NaN NaN];    % could not find anything...
        eye(n).Area = NaN;
    else
        [~,idx] = max(metric);          % pick the circle with best score
        eye(n).Centroid = center(idx,:);
        eye(n).Area = pi*radii(idx)^2;
    end
    if mod(n,100)==0
        fprintf('Frame %d/%d\n',n,size(data,3));
    end
end
Centroid = cell2mat({eye.Centroid}');
Area = cell2mat({eye.Area}');

%% load mworks file
time = '1516';
fn_mworks = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' SubNum '-' date '-' time '.mat'];
load(fn_mworks);
min_hold = 2000;
prepush_frames = 15;
postpush_frames = ceil(min_hold*(frame_rate/1000));
leverDownFrame = cell2mat(input.cLeverDown);
ntrials = length(input.trialOutcomeCell);

%% extract trial by trial data 
%locked to press
area_mat_down = zeros(prepush_frames+postpush_frames, ntrials);
centroid_mat_down = zeros(prepush_frames+postpush_frames,2, ntrials);
rad = sqrt(Area./(4*pi));
for itrial = 1:ntrials
    area_mat_down(:,itrial) = rad(1+leverDownFrame(itrial)-prepush_frames:leverDownFrame(itrial)+postpush_frames,:);
    centroid_mat_down(:,:,itrial) = Centroid(1+leverDownFrame(itrial)-prepush_frames:leverDownFrame(itrial)+postpush_frames,:);
    if (sum(isnan(area_mat_down(:,itrial)),1) > 1) & (sum(isnan(area_mat_down(:,itrial)),1) < size(area_mat_down,1))
        nan_ind = find(isnan(area_mat_down(:,itrial)));
        data_ind = find(~isnan(area_mat_down(:,itrial)));
        if length(find(tsmovavg(isnan(area_mat_down(:,itrial)), 's', 5, 1) == 1))>0
            area_mat_down(:,itrial) = NaN(prepush_frames+postpush_frames, 1);
            centroid_mat_down(:,:,itrial) = NaN(prepush_frames+postpush_frames,2,1);
        else
            for inan = 1:length(nan_ind)
                gap = min(abs(nan_ind(inan)-data_ind),[],1);
                good_ind = find(abs(nan_ind(inan)-data_ind) == gap);
                area_mat_down(nan_ind(inan),itrial) = mean(area_mat_down(data_ind(good_ind),itrial),1);
                centroid_mat_down(nan_ind(inan),:,itrial) = mean(centroid_mat_down(data_ind(good_ind),itrial),1);
            end
        end
    end
end

rad_mat_down = bsxfun(@times, area_mat_down, calib);
centroid_mat_down = bsxfun(@times,centroid_mat_down,calib);
%locked to target
area_mat_target = zeros(pretarget_frames+posttarget_frames, ntrials);
centroid_mat_target = zeros(pretarget_frames+posttarget_frames,2, ntrials);
for itrial = 1:ntrials
    if targetOnFrame(itrial)>0
        area_mat_target(:,itrial) = rad(1+targetOnFrame(itrial)-pretarget_frames:targetOnFrame(itrial)+posttarget_frames,:);
        centroid_mat_target(:,:,itrial) = Centroid(1+targetOnFrame(itrial)-pretarget_frames:targetOnFrame(itrial)+posttarget_frames,:);
        if (sum(isnan(area_mat_target(:,itrial)),1) > 1) & (sum(isnan(area_mat_target(:,itrial)),1) < size(area_mat_target,1))
            nan_ind = find(isnan(area_mat_target(:,itrial)));
            data_ind = find(~isnan(area_mat_target(:,itrial)));
            if length(find(tsmovavg(isnan(area_mat_target(:,itrial)), 's', 5, 1) == 1))>0
                area_mat_target(:,itrial) = NaN(pretarget_frames+posttarget_frames, 1);
                centroid_mat_target(:,:,itrial) = NaN(pretarget_frames+posttarget_frames,2,1);
            else
                for inan = 1:length(nan_ind)
                    gap = min(abs(nan_ind(inan)-data_ind),[],1);
                    good_ind = find(abs(nan_ind(inan)-data_ind) == gap);
                    area_mat_target(nan_ind(inan),itrial) = mean(area_mat_target(data_ind(good_ind),itrial),1);
                    centroid_mat_target(nan_ind(inan),:,itrial) = mean(centroid_mat_target(data_ind(good_ind),itrial),1);
                end
            end
        end
    else
        area_mat_target(:,itrial) = NaN(pretarget_frames+posttarget_frames, 1);
        centroid_mat_target(:,:,itrial) = NaN(pretarget_frames+posttarget_frames,2, 1);
    end
end
rad_mat_target = bsxfun(@times, area_mat_target, calib);
centroid_mat_target = bsxfun(@times,centroid_mat_target,calib);

%% plot eye traces locked to events
close all
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

holdIx = find(cell2mat(input.holdTimesMs)>min_hold);
rad_mat_down_base = bsxfun(@rdivide, rad_mat_down, mean(rad_mat_down(1:15,:),1));
centroid_mat_down_base = bsxfun(@rdivide, centroid_mat_down, mean(centroid_mat_down(1:15,:,:),1));
%plot change in eye area
figure;
tt = (1-prepush_frames:postpush_frames)*(1000/frame_rate);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_down(:,holdIx),2), nanstd(rad_mat_down(:,holdIx),[],2)./sqrt(length(holdIx)));
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_down(1,holdIx),2)).*[0.8 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,holdIx),2), nanstd(rad_mat_down_base(:,holdIx),[],2)./sqrt(length(holdIx)));
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,holdIx),3)), squeeze(nanstd(centroid_mat_down(:,1,holdIx),[],3))./sqrt(length(holdIx)));
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_down(1,1,holdIx),3))).*[0.9 1.1]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,holdIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,holdIx),[],3))./sqrt(length(holdIx)));
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.1])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,holdIx),3)), squeeze(nanstd(centroid_mat_down(:,2,holdIx),[],3))./sqrt(length(holdIx)));
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_down(1,2,holdIx),3))).*[0.9 1.1]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,holdIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,holdIx),[],3))./sqrt(length(holdIx)));
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.1])
suptitle('Align to lever down')
print([fnout '_avg_pressalign.pdf'], '-dpdf');

%plot change in Pupil radius by block type- success only
b1Ix = intersect(holdIx, intersect(find(strcmp(input.trialOutcomeCell,'success')), find(cell2mat(input.tBlock2TrialNumber)==0)));
b2Ix = intersect(holdIx, intersect(find(strcmp(input.trialOutcomeCell,'success')), find(cell2mat(input.tBlock2TrialNumber))));

figure;
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_down(:,b1Ix),2), nanstd(rad_mat_down(:,b1Ix),[],2)./sqrt(length(b1Ix)), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_down(:,b2Ix),2), nanstd(rad_mat_down(:,b2Ix),[],2)./sqrt(length(b2Ix)), '-g');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_down(1,[b1Ix b2Ix]),2)).*[0.8 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,b1Ix),2), nanstd(rad_mat_down_base(:,b1Ix),[],2)./sqrt(length(b1Ix)), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,b2Ix),2), nanstd(rad_mat_down_base(:,b2Ix),[],2)./sqrt(length(b2Ix)), '-g');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_down(:,1,b1Ix),[],3))./sqrt((length(b1Ix))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_down(:,1,b2Ix),[],3))./sqrt((length(b2Ix))), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_down(1,1,[b1Ix b2Ix]),3))).*[0.9 1.1]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_down_base(:,1,b1Ix),[],3))./sqrt((length(b1Ix))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_down_base(:,1,b2Ix),[],3))./sqrt((length(b2Ix))), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.1])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_down(:,2,b1Ix),[],3))./sqrt((length(b1Ix))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_down(:,2,b2Ix),[],3))./sqrt((length(b2Ix))), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_down(1,2,[b1Ix b2Ix]),3))).*[0.9 1.1]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_down_base(:,2,b1Ix),[],3))./sqrt((length(b1Ix))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_down_base(:,2,b2Ix),[],3))./sqrt((length(b2Ix))), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.1])
suptitle('Align to lever down- Black: visual; Green: auditory')
print([fnout '_avg_pressalign_AV.pdf'], '-dpdf');

%plot change in Pupil radius by outcome type
successIx = find(strcmp(input.trialOutcomeCell,'success'));
missedIx = find(strcmp(input.trialOutcomeCell,'ignore'));
figure;
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_down(:,successIx),2), nanstd(rad_mat_down(:,successIx),[],2)./sqrt(length(successIx)), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_down(:,missedIx),2), nanstd(rad_mat_down(:,missedIx),[],2)./sqrt(length(missedIx)), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_down(1,:),2)).*[0.8 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,successIx),2), nanstd(rad_mat_down_base(:,successIx),[],2)./sqrt(length(successIx)), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,missedIx),2), nanstd(rad_mat_down_base(:,missedIx),[],2)./sqrt(length(missedIx)), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down(:,1,successIx),[],3))./sqrt((length(successIx))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,1,missedIx),[],3))./sqrt((length(missedIx))), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_down(1,1,:),3))).*[0.9 1.1]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,successIx),[],3))./sqrt((length(successIx))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,missedIx),[],3))./sqrt((length(missedIx))), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.1])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down(:,2,successIx),[],3))./sqrt((length(successIx))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,2,missedIx),[],3))./sqrt((length(missedIx))), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_down(1,2,:),3))).*[0.9 1.1]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,successIx),[],3))./sqrt((length(successIx))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,missedIx),[],3))./sqrt((length(missedIx))), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.1])
suptitle('Align to lever down- Black: success; Red: missed')
print([fnout '_avg_pressalign_SM.pdf'], '-dpdf');

%plot change in Pupil radius locked to target
rad_mat_target_base = bsxfun(@rdivide, rad_mat_target, mean(rad_mat_down(1:15,:),1));
centroid_mat_target_base = bsxfun(@rdivide, centroid_mat_target, mean(centroid_mat_down(1:15,:,:),1));
figure;
tt = (1-pretarget_frames:posttarget_frames)*(1000/frame_rate);
nonan_trials = sum(~isnan(rad_mat_target(1,:)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_target,2), nanstd(rad_mat_target,[],2)./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.9 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base,2), nanstd(rad_mat_target_base,[],2)./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.9 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,:),3)), squeeze(nanstd(centroid_mat_target(:,1,:),[],3))./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.9 1.3]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,:),3)), squeeze(nanstd(centroid_mat_target_base(:,1,:),[],3))./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.3])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,:),3)), squeeze(nanstd(centroid_mat_target(:,2,:),[],3))./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
suptitle('Align to target')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.9 1.3]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,:),3)), squeeze(nanstd(centroid_mat_target_base(:,2,:),[],3))./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.3])
print([fnout '_avg_targetalign.pdf'], '-dpdf');

%plot change in Pupil radius by block type
b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
b2Ix = find(cell2mat(input.tBlock2TrialNumber));

figure;
targetTrB1 = sum(~isnan(rad_mat_target(1,b1Ix)),2);
targetTrB2 = sum(~isnan(rad_mat_target(1,b2Ix)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_target(:,b1Ix),2), nanstd(rad_mat_target(:,b1Ix),[],2)./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target(:,b2Ix),2), nanstd(rad_mat_target(:,b2Ix),[],2)./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.9 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b1Ix),2), nanstd(rad_mat_target_base(:,b1Ix),[],2)./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b2Ix),2), nanstd(rad_mat_target_base(:,b2Ix),[],2)./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.9 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.9 1.3]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.3])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.9 1.3]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.3])
suptitle('Align to target- Black: visual; Green: auditory')
print([fnout '_avg_targetalign_AV.pdf'], '-dpdf');

%plot change in Pupil radius by block type- success only
b1Ix = intersect(find(strcmp(input.trialOutcomeCell,'success')), find(cell2mat(input.tBlock2TrialNumber)==0));
b2Ix = intersect(find(strcmp(input.trialOutcomeCell,'success')), find(cell2mat(input.tBlock2TrialNumber)));

figure;

targetTrB1 = sum(~isnan(rad_mat_target(1,b1Ix)),2);
targetTrB2 = sum(~isnan(rad_mat_target(1,b2Ix)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_target(:,b1Ix),2), nanstd(rad_mat_target(:,b1Ix),[],2)./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target(:,b2Ix),2), nanstd(rad_mat_target(:,b2Ix),[],2)./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.9 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b1Ix),2), nanstd(rad_mat_target_base(:,b1Ix),[],2)./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b2Ix),2), nanstd(rad_mat_target_base(:,b2Ix),[],2)./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.9 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.9 1.3]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.3])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.9 1.3]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.3])
suptitle('Align to target- Success only Black: visual; Green: auditory')
print([fnout '_avg_targetalign_AV_Sonly.pdf'], '-dpdf');


%plot change in Pupil radius by block type
successIx = find(strcmp(input.trialOutcomeCell,'success'));
missedIx = find(strcmp(input.trialOutcomeCell,'ignore'));

figure;

targetTrS = sum(~isnan(rad_mat_target(1,successIx)),2);
targetTrM = sum(~isnan(rad_mat_target(1,missedIx)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_target(:,successIx),2), nanstd(rad_mat_target(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target(:,missedIx),2), nanstd(rad_mat_target(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.9 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.9 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.9 1.3]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.3])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.9 1.3]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.3])
suptitle('Align to target- Black: success; Red: missed')
print([fnout '_avg_targetalign_SM.pdf'], '-dpdf');

%hit and miss for V trials only
b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
b2Ix = find(cell2mat(input.tBlock2TrialNumber));
successIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'success')));
missedIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'ignore')));

figure;

targetTrS = sum(~isnan(rad_mat_target(1,successIx)),2);
targetTrM = sum(~isnan(rad_mat_target(1,missedIx)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_target(:,successIx),2), nanstd(rad_mat_target(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target(:,missedIx),2), nanstd(rad_mat_target(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.9 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.9 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.9 1.3]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.3])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.9 1.3]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.3])
suptitle('Align to target- Visual only- Black: success; Red: missed')
print([fnout '_avg_targetalign_SM_Vonly.pdf'], '-dpdf');

%hit and miss for A trials only
successIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
missedIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'ignore')));

figure;

targetTrS = sum(~isnan(rad_mat_target(1,successIx)),2);
targetTrM = sum(~isnan(rad_mat_target(1,missedIx)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_target(:,successIx),2), nanstd(rad_mat_target(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target(:,missedIx),2), nanstd(rad_mat_target(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.9 1.1]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.9 1.1])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.9 1.3]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.9 1.3])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.9 1.3]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.9 1.3])
suptitle('Align to target- Auditory only- Black: success; Red: missed')
print([fnout '_avg_targetalign_SM_Aonly.pdf'], '-dpdf');

%% plot pupil dynamics
%normalize to max area (for whole experiment)
rad_max = sqrt(max(Area,[],1)/(4*pi))*calib;
rad_mat_target_norm = rad_mat_target./rad_max;
rad_mat_down_norm = rad_mat_down./rad_max;

b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
b2Ix = find(cell2mat(input.tBlock2TrialNumber));
reactTimes = cell2mat(input.reactTimesMs);
figure;
min_x = min(nanmean(rad_mat_target_norm(1:15,:),1),[],2)*0.9;
max_x = max(nanmean(rad_mat_target_norm(1:15,:),1),[],2)*1.1;
subplot(2,2,1)
scatter(nanmean(rad_mat_target_norm(1:15,b1Ix),1), reactTimes(:,b1Ix),'ok');
hold on
successIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
scatter(nanmean(rad_mat_target_norm(1:15,b2Ix),1), reactTimes(:,b2Ix),'og');
ylim([0 600])
xlim([min_x max_x])
title('Pupil size before target')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
subplot(2,2,2)
scatter(nanmean(rad_mat_down_norm(1:15,b1Ix),1), reactTimes(:,b1Ix),'ok');
hold on
scatter(nanmean(rad_mat_down_norm(1:15,b2Ix),1), reactTimes(:,b2Ix),'og');
ylim([0 600])
xlim([min_x max_x])
title('Pupil size before press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
[n, edges] = histcounts([nanmean(rad_mat_target_norm(1:15,success1Ix),1) nanmean(rad_mat_target_norm(1:15,success2Ix),1)],8);
[n_b1, bin_b1] = histc(nanmean(rad_mat_target_norm(1:15,success1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_target_norm(1:15,success2Ix),1),edges);
subplot(2,2,3)
for i = 1:length(edges)
    ind = find(bin_b1 == i-1);
    ploterr(nanmean(nanmean(rad_mat_target_norm(1:15,success1Ix(ind)),1),2), mean(reactTimes(:,success1Ix(ind)),2), std(nanmean(rad_mat_target_norm(1:15,success1Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success1Ix(ind)),[],2)./sqrt(length(ind)), 'ok')
    hold on;
    ind = find(bin_b2 == i-1);
    ploterr(nanmean(nanmean(rad_mat_target_norm(1:15,success2Ix(ind)),1),2), mean(reactTimes(:,success2Ix(ind)),2), std(nanmean(rad_mat_target_norm(1:15,success2Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success2Ix(ind)),[],2)./sqrt(length(ind)), 'og')
end
ylim([0 600])
xlim([min_x max_x])
ylabel('React time (ms)')
xlabel('Normalized pupil size')
title('Average before target')
success1Ix = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'success')));
success2Ix = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
[n, edges] = histcounts([nanmean(rad_mat_down_norm(1:15,success1Ix),1) nanmean(rad_mat_down_norm(1:15,success2Ix),1)],8);
[n_b1, bin_b1] = histc(nanmean(rad_mat_down_norm(1:15,success1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_down_norm(1:15,success2Ix),1),edges);
subplot(2,2,4)
for i = 1:length(edges)
    ind = find(bin_b1 == i-1);
    ploterr(nanmean(nanmean(rad_mat_down_norm(1:15,success1Ix(ind)),1),2), mean(reactTimes(:,success1Ix(ind)),2), std(nanmean(rad_mat_down_norm(1:15,success1Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success1Ix(ind)),[],2)./sqrt(length(ind)), 'ok')
    hold on;
    ind = find(bin_b2 == i-1);
    ploterr(nanmean(nanmean(rad_mat_down_norm(1:15,success2Ix(ind)),1),2), mean(reactTimes(:,success2Ix(ind)),2), std(nanmean(rad_mat_down_norm(1:15,success2Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success2Ix(ind)),[],2)./sqrt(length(ind)), 'og')
end
ylim([0 600])
xlim([min_x max_x])
title('Average before press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
suptitle([mouse ' ' date ' Pupil size and react time'])
print([fnout '_PupilvsReact.pdf'], '-dpdf');    


%distribution of pupil size
rad_all = sqrt(Area/(4*pi))*calib;
cleverup = cell2mat(input.cLeverUp);
cleverdown = cell2mat(input.cLeverDown);
tr_frames = [];
for itrial = 1:ntrials
    tr_frames = [tr_frames double([cleverdown(itrial):cleverup(itrial)])];
end
min_x = min(rad_all,[],1)*0.9;
max_x = max(rad_all,[],1)*1.1;
figure;
subplot(2,2,1)
hist(rad_all,100)
title('All Frames')
xlim([min_x max_x])
subplot(2,2,2)
hist(rad_all(tr_frames,:),50)
title('Trial Frames')
iti_frames = find(ismember(1:size(rad_all,1), tr_frames)==0);
xlim([min_x max_x])
subplot(2,2,3)
hist(rad_all(iti_frames,:),100)
title('ITI Frames')
xlim([min_x max_x])
subplot(2,2,4)
hist(rad_all(tr_frames,:)./max(rad_all,[],1),50)
xlim([0 1])
title('Trial Frames- norm to max pupil size')
print([fnout '_pupil_size_hist.pdf'], '-dpdf');
save([fnout '_pupil.mat'], 'Area', 'Centroid', 'frame_rate')