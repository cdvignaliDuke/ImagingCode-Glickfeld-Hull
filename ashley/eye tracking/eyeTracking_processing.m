SubNum = '614';
date = '150623';
ImgFolder = '004';
mouse = 'AW14';
fName = '004_000_000';

fnout = '\\CRASH.dhe.duke.edu\data\home\lindsey
% % load MWorks file
% CD = ['Z:\data\' mouse '\mworks\' date];
% cd(CD);
% mworks = ['data-' 'i' SubNum '-' date '-' time]; 
% load (mworks);

% Set current directory to crash folder
CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

%% 
fn = [fName '_eye.mat'];
load(fn);          % should be a '*_eye.mat' file
rad_range = [8 15];
    
data = squeeze(data);      % the raw images...
xc = size(data,2)/2;       % image center
yc = size(data,1)/2;
W=40;

data = data(yc-W:yc+W,xc-W:xc+W,:);
warning off;

%Edit Patrick: changed this so that struct doesn't change size every time
%while maintaining backwards compatibility

A = cell(size(data,3),1);
B = cell(size(data,3),1);
for n = 1:size(data,3)
    A{n} = [0,0];
    B{n} = [0];
end

eye = struct('Centroid',A,'Area',B);
for(n=1:size(data,3))
    [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.9);
    if(isempty(center))
        eye(n).Centroid = [NaN NaN];    % could not find anything...
        eye(n).Area = NaN;
    else
        [~,idx] = max(metric);          % pick the circle with best score
        eye(n).Centroid = center(idx,:);
        eye(n).Area = 4*pi*radii(idx)^2;
    end
    
    if mod(n,100)==0
        fprintf('Frame %d/%d\n',n,size(data,3));
    end
end

Centroid = cell2mat({eye.Centroid}');
Area = cell2mat({eye.Area}');

%% load mworks file
time = '1217';
fn_mworks = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\mworks\' date '\data-i' SubNum '-' date '-' time '.mat'];
load(fn_mworks);
prepush_frames = 15;
postpush_frames = 45;
leverDownFrame = cell2mat(input.cLeverDown);
ntrials = length(input.trialOutcomeCell);
area_mat_down = zeros(prepush_frames+postpush_frames, ntrials);
centroid_mat_down = zeros(prepush_frames+postpush_frames,2, ntrials);
for itrial = 1:ntrials
    area_mat_down(:,itrial) = Area(1+leverDownFrame(itrial)-prepush_frames:leverDownFrame(itrial)+postpush_frames,:);
    centroid_mat_down(:,:,itrial) = Centroid(1+leverDownFrame(itrial)-prepush_frames:leverDownFrame(itrial)+postpush_frames,:);
end

%plot change in eye area
figure;
tt = 1-prepush_frames:postpush_frames;
subplot(3,1,1)
shadedErrorBar(tt,nanmean(area_mat_down,2), nanstd(area_mat_down,[],2)./sqrt(ntrials));
xlabel('Frames')
ylabel('Eye area')
subplot(3,1,2)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,:),3)), squeeze(nanstd(centroid_mat_down(:,1,:),[],3))./sqrt(ntrials));
xlabel('Frames')
ylabel('Eye position- Horizontal')
subplot(3,1,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,:),3)), squeeze(nanstd(centroid_mat_down(:,2,:),[],3))./sqrt(ntrials));
xlabel('Frames')
ylabel('Eye position- Vertical')
suptitle('Align to lever down')

%plot change in eye area by block type
b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
b2Ix = find(cell2mat(input.tBlock2TrialNumber));

figure;
tt = 1-prepush_frames:postpush_frames;
subplot(3,1,1)
shadedErrorBar(tt,nanmean(area_mat_down(:,b1Ix),2), nanstd(area_mat_down(:,b1Ix),[],2)./sqrt(length(b1Ix)), '-k');
hold on
shadedErrorBar(tt,nanmean(area_mat_down(:,b2Ix),2), nanstd(area_mat_down(:,b2Ix),[],2)./sqrt(length(b2Ix)), '-g');
xlabel('Frames')
ylabel('Eye area')
subplot(3,1,2)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_down(:,1,b1Ix),[],3))./sqrt((length(b1Ix))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_down(:,1,b2Ix),[],3))./sqrt((length(b2Ix))), '-g');
xlabel('Frames')
ylabel('Eye position- Horizontal')
subplot(3,1,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_down(:,2,b1Ix),[],3))./sqrt((length(b1Ix))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_down(:,2,b2Ix),[],3))./sqrt((length(b2Ix))), '-g');
xlabel('Frames')
ylabel('Eye position- Vertical')
suptitle('Align to lever down')

%change with target
targetOnFrame = cell2mat_padded(input.cTargetOn);
targetInd = find(targetOnFrame>0);
targetOnFrame(find(targetOnFrame==0)) = [];
targetTrials = size(targetOnFrame,1);
area_mat_target = zeros(prepush_frames+postpush_frames, targetTrials);
centroid_mat_target = zeros(prepush_frames+postpush_frames,2, targetTrials);
for itrial = 1:targetTrials
    area_mat_target(:,itrial) = Area(1+targetOnFrame(itrial)-prepush_frames:targetOnFrame(itrial)+postpush_frames,:);
    centroid_mat_target(:,:,itrial) = Centroid(1+targetOnFrame(itrial)-prepush_frames:targetOnFrame(itrial)+postpush_frames,:);
end

%plot change in eye area
figure;
tt = 1-prepush_frames:postpush_frames;
subplot(3,1,1)
shadedErrorBar(tt,nanmean(area_mat_target,2), nanstd(area_mat_target,[],2)./sqrt(targetTrials));
xlabel('Frames')
ylabel('Eye area')
subplot(3,1,2)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,:),3)), squeeze(nanstd(centroid_mat_target(:,1,:),[],3))./sqrt(targetTrials));
xlabel('Frames')
ylabel('Eye position- Horizontal')
subplot(3,1,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,:),3)), squeeze(nanstd(centroid_mat_target(:,2,:),[],3))./sqrt(targetTrials));
xlabel('Frames')
ylabel('Eye position- Vertical')
suptitle('Align to target')

%plot change in eye area by block type
b1Ix = intersect(targetInd, find(cell2mat(input.tBlock2TrialNumber)==0));
b2Ix = intersect(targetInd, find(cell2mat(input.tBlock2TrialNumber)));

figure;
tt = 1-prepush_frames:postpush_frames;
subplot(3,1,1)
shadedErrorBar(tt,nanmean(area_mat_target(:,b1Ix),2), nanstd(area_mat_target(:,b1Ix),[],2)./sqrt(length(b1Ix)), '-k');
hold on
shadedErrorBar(tt,nanmean(area_mat_target(:,b2Ix),2), nanstd(area_mat_target(:,b2Ix),[],2)./sqrt(length(b2Ix)), '-g');
xlabel('Frames')
ylabel('Eye area')
subplot(3,1,2)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b1Ix),[],3))./sqrt((length(b1Ix))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b2Ix),[],3))./sqrt((length(b2Ix))), '-g');
xlabel('Frames')
ylabel('Eye position- Horizontal')
subplot(3,1,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b1Ix),[],3))./sqrt((length(b1Ix))), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b2Ix),[],3))./sqrt((length(b2Ix))), '-g');
xlabel('Frames')
ylabel('Eye position- Vertical')
suptitle('Align to target')