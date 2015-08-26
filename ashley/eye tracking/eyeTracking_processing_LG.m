clear all
close all
SubNum = '614';
date = '150818';
runs = ['001'; '002'];
time_mat = ['1141'; '1216'];
mouse = 'AW14';
frame_rate = 15;
calib = 1/26.6; %mm per pixel
nrun = size(runs,1);
%% load and combine mworks files
for irun = 1:nrun
    time = time_mat(irun,:);
    fn_mworks = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' SubNum '-' date '-' time '.mat'];
    if irun == 1
        input = mwLoadData(fn_mworks, [], []);
    else
        input = [input mwLoadData(fn_mworks, [], [])];
    end
end
input = concatenateDataBlocks(input);
    
runstr = runs(1,:);
if nrun>1
    for irun = 2:nrun
        runstr = [runstr '-' runs(irun,:)];
    end
end
fnout = ['Z:\home\lindsey\Analysis\Behavior\EyeTracking\' mouse '-' date '\' mouse '-' date '-' runstr];

%% 
min_hold = 2000;
prepush_frames = 15;
postpush_frames = ceil(min_hold*(frame_rate/1000));
prerelease_frames = 15;
postrelease_frames = ceil(1500*(frame_rate/1000));
pretarget_frames = 15;
posttarget_frames = ceil(4000*(frame_rate/1000));

%% Load and combine eye tracking data
% Set current directory to crash folder
Area = {};
Centroid = {};
Eye_data = {};
for irun =  1:nrun
    CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\eye tracking\' date '\' runs(irun,:)];
    cd(CD);
    fn = [runs(irun,:) '_000_000_eye.mat'];
    load(fn);          % should be a '*_eye.mat' file
    
    data = squeeze(data);      % the raw images...
    xc = size(data,2)/2;       % image center
    yc = size(data,1)/2;
    W=40;

    rad_range = [15 30];
    data = data(yc-W:yc+W,xc-W:xc+W,:);
    warning off;
    
    A = cell(size(data,3),1);
    B = cell(size(data,3),1);
    for n = 1:size(data,3)
        A{n} = [0,0];
        B{n} = [0];
    end
    eye = struct('Centroid',A,'Area',B);
    radii = [];
    for n = 1:size(data,3)
        if ~isempty(radii)
            rad_range = [radii-3 radii+3];
            rad_range(find(rad_range<15)) = 15;
            rad_range(find(rad_range>30)) = 30;
        else
            rad_range = [15 30];
        end
        [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),round(rad_range),'Sensitivity',0.9,'Method', 'TwoStage');
        if(isempty(center))
            eye(n).Centroid = [NaN NaN];    % could not find anything...
            eye(n).Area = NaN;
        else
            [~,idx] = max(metric);          % pick the circle with best score
            eye(n).Centroid = center(idx,:);
            eye(n).Area = pi*radii(idx)^2;
            radii = radii(idx);
        end
        if mod(n,100)==0
            fprintf('Frame %d/%d\n',n,size(data,3));
        end
    end
    Centroid{irun} = cell2mat({eye.Centroid}');
    Area{irun} = cell2mat({eye.Area}');
    Eye_data{irun} = data;
end

%% reset frame counter
run_trials = input.trialsSinceReset;
cLeverDown = cell2mat(input.cLeverDown);
cLeverUp = cell2mat(input.cLeverUp);
cTargetOn = celleqel2mat_padded(input.cTargetOn);
cItiStart = cell2mat(input.cItiStart);
Area_temp = [];
Centroid_temp = [];
Eye_data_temp = [];
for irun = 1:nrun
    if irun < nrun
        offset = size(Area{irun},1);
        startTrial = run_trials(irun)+1;
        endTrial = run_trials(irun)+run_trials(irun+1);
        cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
        cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
        cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
        cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
    end
    Area_temp = [Area_temp; Area{irun}];
    Centroid_temp = [Centroid_temp; Centroid{irun}];
    %Eye_data_temp = cat(3, Eye_data_temp, Eye_data{irun});
end
clear Eye_data;
ntrials = length(input.trialOutcomeCell);

%% no measurement frames
figure; 
x = find((Area_temp>1200));
start = 1;
frames = sort(randsample(length(x),100));
for i = 1:100
    subplot(10,10,start);
    imagesq(Eye_data_temp(:,:,x(frames(i)))); 
    title(x(frames(i)))
    hold on
    plot(Centroid_temp(x(frames(i)),1), Centroid_temp(x(frames(i)),2), 'ok')
    hold on
    plot(Centroid_temp(x(frames(i)),1)+sqrt(Area_temp(x(frames(i)),1)/pi), Centroid_temp(x(frames(i)),2), 'ok')
    start = start+1;
end
print([fnout '_nanframes.pdf'], '-dpdf');

%% Remove NaNs if sparse and align to start and target
nanrun = ceil(500*(frame_rate/1000));
Rad_temp = sqrt(Area_temp./pi);
rad_mat_down = zeros(prepush_frames+postpush_frames, ntrials);
centroid_mat_down = zeros(prepush_frames+postpush_frames,2, ntrials);
rad_mat_up = zeros(prerelease_frames+postrelease_frames, ntrials);
centroid_mat_up = zeros(prerelease_frames+postrelease_frames,2, ntrials);
nframes = size(Rad_temp,1);
for itrial = 1:ntrials
    if itrial == ntrials
        crange = [double(input.cItiStart{itrial}):nframes];
    else
        if double(input.cItiStart{itrial})< 1
            crange = [1:double(input.cItiStart{itrial+1}-1)];
        else
            crange = [double(input.cItiStart{itrial}): double(input.cItiStart{itrial+1}-1)];
        end
        if sum(isnan(Rad_temp(crange,1)),2)>0
            if length(find(tsmovavg(isnan(Rad_temp(crange,1)), 's', nanrun, 1) == 1))> 0
                Rad_temp(crange,1) = NaN(length(crange),1);
            else
                nanind = find(isnan(Rad_temp(crange,1)));
                dataind = find(~isnan(Rad_temp(crange,1)));
                for inan = 1:length(nan_ind)
                    gap = min(abs(nan_ind(inan)-data_ind),[],1);
                    good_ind = find(abs(nan_ind(inan)-data_ind) == gap);
                    Rad_temp(nan_ind(inan),1) = mean(Rad_temp(data_ind(good_ind),1),1);
                    Centroid_temp(nan_ind(inan),:) = mean(Centroid_temp(data_ind(good_ind),:),1);
                end
            end
        end
    end
    rad_mat_down(:,itrial) = Rad_temp(1+cLeverDown(itrial)-prepush_frames:cLeverDown(itrial)+postpush_frames,:);
    centroid_mat_down(:,:,itrial) = Centroid_temp(1+cLeverDown(itrial)-prepush_frames:cLeverDown(itrial)+postpush_frames,:);
    rad_mat_up(:,itrial) = Rad_temp(1+cLeverUp(itrial)-prerelease_frames:cLeverUp(itrial)+postrelease_frames,:);
    centroid_mat_up(:,:,itrial) = Centroid_temp(1+cLeverUp(itrial)-prerelease_frames:cLeverUp(itrial)+postrelease_frames,:);
    if cTargetOn(itrial)>0 & cTargetOn(itrial)+posttarget_frames < nframes
        rad_mat_target(:,itrial) = Rad_temp(1+cTargetOn(itrial)-pretarget_frames:cTargetOn(itrial)+posttarget_frames,:);
        centroid_mat_target(:,:,itrial) = Centroid_temp(1+cTargetOn(itrial)-pretarget_frames:cTargetOn(itrial)+posttarget_frames,:);
    else
        rad_mat_target(:,itrial) = NaN(pretarget_frames+posttarget_frames, 1);
        centroid_mat_target(:,:,itrial) = NaN(pretarget_frames+posttarget_frames,2, 1);
    end
end
rad_mat_down = bsxfun(@times, rad_mat_down, calib);
centroid_mat_down = bsxfun(@times,centroid_mat_down,calib);
rad_mat_up = bsxfun(@times, rad_mat_up, calib);
centroid_mat_up = bsxfun(@times,centroid_mat_up,calib);
rad_mat_target = bsxfun(@times, rad_mat_target, calib);
centroid_mat_target = bsxfun(@times,centroid_mat_target,calib);        

%% plot eye traces align to press
close all
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

holdIx = find(cell2mat(input.holdTimesMs)>min_hold);
rad_mat_down_base = bsxfun(@rdivide, rad_mat_down, mean(rad_mat_down(1:15,:),1));
centroid_mat_down_base = bsxfun(@rdivide, centroid_mat_down, mean(centroid_mat_down(1:15,:,:),1));

%plot change in eye area align to press
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

%hit and miss for V trials only
b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
b2Ix = find(cell2mat(input.tBlock2TrialNumber));
successIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'success')));
missedIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'ignore')));

figure;

downTrS = sum(~isnan(rad_mat_down(1,successIx)),2);
downTrM = sum(~isnan(rad_mat_down(1,missedIx)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_down(:,successIx),2), nanstd(rad_mat_down(:,successIx),[],2)./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_down(:,missedIx),2), nanstd(rad_mat_down(:,missedIx),[],2)./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_down(1,:),2)).*[0.8 1.4]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,successIx),2), nanstd(rad_mat_down_base(:,successIx),[],2)./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,missedIx),2), nanstd(rad_mat_down_base(:,missedIx),[],2)./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.4])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down(:,1,successIx),[],3))./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,1,missedIx),[],3))./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_down(1,1,:),3))).*[0.8 1.4]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,successIx),[],3))./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,missedIx),[],3))./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.8 1.4])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down(:,2,successIx),[],3))./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,2,missedIx),[],3))./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_down(1,2,:),3))).*[0.8 1.4]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,successIx),[],3))./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,missedIx),[],3))./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.8 1.4])
suptitle('Align to press- Visual only- Black: success; Red: missed')
print([fnout '_avg_pressalign_SM_Vonly.pdf'], '-dpdf');

%hit and miss for A trials only
successIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
missedIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'ignore')));

figure;

downTrS = sum(~isnan(rad_mat_down(1,successIx)),2);
downTrM = sum(~isnan(rad_mat_down(1,missedIx)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_down(:,successIx),2), nanstd(rad_mat_down(:,successIx),[],2)./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_down(:,missedIx),2), nanstd(rad_mat_down(:,missedIx),[],2)./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_down(1,:),2)).*[0.8 1.4]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,successIx),2), nanstd(rad_mat_down_base(:,successIx),[],2)./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_down_base(:,missedIx),2), nanstd(rad_mat_down_base(:,missedIx),[],2)./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.4])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down(:,1,successIx),[],3))./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,1,missedIx),[],3))./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_down(1,1,:),3))).*[0.8 1.4]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,successIx),[],3))./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,missedIx),[],3))./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.8 1.4])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down(:,2,successIx),[],3))./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,2,missedIx),[],3))./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_down(1,2,:),3))).*[0.8 1.4]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,successIx),[],3))./sqrt(downTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,missedIx),[],3))./sqrt(downTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.8 1.4])
suptitle('Align to press- Auditory only- Black: success; Red: missed')
print([fnout '_avg_pressalign_SM_Aonly.pdf'], '-dpdf');

%% plot change in pupil radius locked to lever up
rad_mat_up_base = bsxfun(@rdivide, rad_mat_up, mean(rad_mat_down(1:15,:),1));
centroid_mat_up_base = bsxfun(@rdivide, centroid_mat_up, mean(centroid_mat_down(1:15,:,:),1));

tt = (1-prerelease_frames:postrelease_frames).*(1000/frame_rate);

FIx = find(strcmp(input.trialOutcomeCell, 'failure'));
SIx = find(strcmp(input.trialOutcomeCell, 'success'));
MIx = find(strcmp(input.trialOutcomeCell, 'ignore'));
FIxlong = intersect(find(cell2mat(input.tCyclesOn)>3), FIx);
SIxlong = intersect(find(cell2mat(input.tCyclesOn)>3), SIx);
b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
b2Ix = find(cell2mat(input.tBlock2TrialNumber)==1);
Fb1Ix = intersect(b1Ix, FIxlong);
Fb2Ix = intersect(b2Ix, FIxlong);
Sb1Ix = intersect(b1Ix, SIxlong);
Sb2Ix = intersect(b2Ix, SIxlong);

figure;
subplot(2,2,1)
shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Fb1Ix),2), nanstd(rad_mat_up_base(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Sb1Ix),2), nanstd(rad_mat_up_base(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k'); 
title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
xlim([-500 1500])
subplot(2,2,2)
shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Fb2Ix),2), nanstd(rad_mat_up_base(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Sb2Ix),2), nanstd(rad_mat_up_base(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
xlim([-500 1500])
subplot(2,2,3)
shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Fb1Ix),2), nanstd(rad_mat_up_base(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Fb2Ix),2), nanstd(rad_mat_up_base(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'g'); 
title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
xlim([-500 1500])
subplot(2,2,4)
shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Sb1Ix),2), nanstd(rad_mat_up_base(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Sb2Ix),2), nanstd(rad_mat_up_base(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'g'); 
alignYaxes
title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
xlim([-500 1500])
print([fnout '_avg_releasealign_SM_AV.pdf'], '-dpdf');

%% plot change in Pupil radius locked to target
rad_mat_target_base = bsxfun(@rdivide, rad_mat_target, mean(rad_mat_down(1:15,:),1));
centroid_mat_target_base = bsxfun(@rdivide, centroid_mat_target, mean(centroid_mat_down(1:15,:,:),1));
figure;
tt = (1-pretarget_frames:posttarget_frames)*(1000/frame_rate);
nonan_trials = sum(~isnan(rad_mat_target(1,:)),2);
subplot(3,2,1)
shadedErrorBar(tt,nanmean(rad_mat_target,2), nanstd(rad_mat_target,[],2)./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base,2), nanstd(rad_mat_target_base,[],2)./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.4])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,:),3)), squeeze(nanstd(centroid_mat_target(:,1,:),[],3))./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,:),3)), squeeze(nanstd(centroid_mat_target_base(:,1,:),[],3))./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.8 1.4])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,:),3)), squeeze(nanstd(centroid_mat_target(:,2,:),[],3))./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
suptitle('Align to target')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,:),3)), squeeze(nanstd(centroid_mat_target_base(:,2,:),[],3))./sqrt(nonan_trials));
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.8 1.4])
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
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b1Ix),2), nanstd(rad_mat_target_base(:,b1Ix),[],2)./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b2Ix),2), nanstd(rad_mat_target_base(:,b2Ix),[],2)./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.4])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.8 1.4])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.8 1.4])
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
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b1Ix),2), nanstd(rad_mat_target_base(:,b1Ix),[],2)./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b2Ix),2), nanstd(rad_mat_target_base(:,b2Ix),[],2)./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.4])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.8 1.4])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b1Ix),[],3))./sqrt(targetTrB1), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-g');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.8 1.4])
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
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.4])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.8 1.4])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.8 1.4])
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
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.4])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.8 1.4])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.8 1.4])
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
ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
subplot(3,2,2)
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Pupil radius')
ylim([0.8 1.4])
subplot(3,2,3)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
subplot(3,2,4)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Horizontal')
ylim([0.8 1.4])
subplot(3,2,5)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
subplot(3,2,6)
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
hold on
shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
xlabel('Time (ms)')
ylabel('Eye position- Vertical')
ylim([0.8 1.4])
suptitle('Align to target- Auditory only- Black: success; Red: missed')
print([fnout '_avg_targetalign_SM_Aonly.pdf'], '-dpdf');

%% plot pupil dynamics
%normalize to max area (for whole experiment)
rad_max = sqrt(max(Area_temp,[],1)/(pi))*calib;
rad_mat_target_norm = rad_mat_target./rad_max;
rad_mat_down_norm = rad_mat_down./rad_max;

b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
b2Ix = find(cell2mat(input.tBlock2TrialNumber));
reactTimes = cell2mat(input.reactTimesMs);
figure;
min_x = min(nanmean(rad_mat_target_norm(1:15,:),1),[],2)*0.9;
max_x = max(nanmean(rad_mat_target_norm(1:15,:),1),[],2)*1.1;
subplot(2,3,3)
scatter(nanmean(rad_mat_target_norm(1:15,b1Ix),1), reactTimes(:,b1Ix),'ok');
hold on
successIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
scatter(nanmean(rad_mat_target_norm(1:15,b2Ix),1), reactTimes(:,b2Ix),'og');
ylim([0 600])
xlim([min_x max_x])
title('Pupil size before target')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
subplot(2,3,1)
scatter(nanmean(rad_mat_down_norm(1:15,b1Ix),1), reactTimes(:,b1Ix),'ok');
hold on
scatter(nanmean(rad_mat_down_norm(1:15,b2Ix),1), reactTimes(:,b2Ix),'og');
ylim([0 600])
xlim([min_x max_x])
title('Pupil size before press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
subplot(2,3,2)
scatter(nanmean(rad_mat_down_norm(16:30,b1Ix),1), reactTimes(:,b1Ix),'ok');
hold on
scatter(nanmean(rad_mat_down_norm(16:30,b2Ix),1), reactTimes(:,b2Ix),'og');
ylim([0 600])
xlim([min_x max_x])
title('Pupil size after press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')

success1Ix = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'success')));
success2Ix = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
[n, edges] = histcounts([nanmean(rad_mat_target_norm(1:15,success1Ix),1) nanmean(rad_mat_target_norm(1:15,success2Ix),1)],8);
[n_b1, bin_b1] = histc(nanmean(rad_mat_target_norm(1:15,success1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_target_norm(1:15,success2Ix),1),edges);
subplot(2,3,6)
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
[n, edges] = histcounts([nanmean(rad_mat_down_norm(1:15,success1Ix),1) nanmean(rad_mat_down_norm(1:15,success2Ix),1)],8);
[n_b1, bin_b1] = histc(nanmean(rad_mat_down_norm(1:15,success1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_down_norm(1:15,success2Ix),1),edges);
subplot(2,3,4)
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
[n, edges] = histcounts([nanmean(rad_mat_down_norm(16:30,success1Ix),1) nanmean(rad_mat_down_norm(16:30,success2Ix),1)],8);
[n_b1, bin_b1] = histc(nanmean(rad_mat_down_norm(16:30,success1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_down_norm(16:30,success2Ix),1),edges);
subplot(2,3,5)
for i = 1:length(edges)
    ind = find(bin_b1 == i-1);
    ploterr(nanmean(nanmean(rad_mat_down_norm(16:30,success1Ix(ind)),1),2), mean(reactTimes(:,success1Ix(ind)),2), std(nanmean(rad_mat_down_norm(16:30,success1Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success1Ix(ind)),[],2)./sqrt(length(ind)), 'ok')
    hold on;
    ind = find(bin_b2 == i-1);
    ploterr(nanmean(nanmean(rad_mat_down_norm(16:30,success2Ix(ind)),1),2), mean(reactTimes(:,success2Ix(ind)),2), std(nanmean(rad_mat_down_norm(16:30,success2Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success2Ix(ind)),[],2)./sqrt(length(ind)), 'og')
end
ylim([0 600])
xlim([min_x max_x])
title('Average after press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
suptitle([mouse ' ' date ' Pupil size and react time'])
print([fnout '_PupilvsReact.pdf'], '-dpdf');    

%hit rate vs pupil size- before press
[n, edges] = histcounts([nanmean(rad_mat_down_norm(1:15,:),1)],5);
smIx = strcmp(input.trialOutcomeCell,'success') + strcmp(input.trialOutcomeCell,'ignore');
successIx = strcmp(input.trialOutcomeCell,'success');
missedIx = strcmp(input.trialOutcomeCell,'ignore');
sm1Ix = intersect(find(~isnan(rad_mat_down_norm(1,:))), intersect(b1Ix, find(smIx)));
sm2Ix = intersect(find(~isnan(rad_mat_down_norm(1,:))), intersect(b2Ix, find(smIx)));
[n_b1, bin_b1] = histc(nanmean(rad_mat_down_norm(1:15,sm1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_down_norm(1:15,sm2Ix),1),edges);
figure;
subplot(3,2,1)
for i = 1:length(edges)
    ind1 = find(bin_b1 == i);
    if (sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2))>5
        [pct1, err1] = binofit(sum(successIx(sm1Ix(ind1)),2), sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        ploterr(nanmean(nanmean(rad_mat_down_norm(1:15,sm1Ix(ind1)),1),2), pct1, std(nanmean(rad_mat_down_norm(1:15,sm1Ix(ind1)),1),[],2)./sqrt(length(ind1)), pct1-err1(1), 'ok')
        hold on;
    end
    ind2 = find(bin_b2 == i);
    if sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2)>5
        [pct2, err2] = binofit(sum(successIx(sm2Ix(ind2)),2), sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        ploterr(nanmean(nanmean(rad_mat_down_norm(1:15,sm2Ix(ind2)),1),2), pct2, std(nanmean(rad_mat_down_norm(1:15,sm2Ix(ind2)),1),[],2)./sqrt(length(ind2)), pct2-err2(1), 'og')
        hold on
    end
end
xlim([0 1])
ylim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')

tGratingDirection = cell2mat(input.tGratingDirectionDeg);
Oris = unique(tGratingDirection);
nOri = length(Oris);
colmat = strvcat('k', 'b', 'g', 'y', 'r', 'm');
ori_rad_mat = zeros(nOri,length(edges));
ori_hit_mat = zeros(nOri,length(edges));
ori_n_mat = zeros(nOri,length(edges));
for iOri = 2:nOri
    ori_ind = find(tGratingDirection(sm1Ix) == Oris(iOri));
    for  i = 1:length(edges)
        ind1 = intersect(ori_ind, find(bin_b1 == i));
        ori_rad_mat(iOri,i) = nanmean(nanmean(rad_mat_down_norm(1:15,sm1Ix(ind1)),1),2);
        ori_hit_mat(iOri,i) = sum(successIx(sm1Ix(ind1)),2)/(sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        x = (sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        if isempty(x)
            x = 0;
        end
        ori_n_mat(iOri,i) = x;
    end
    subplot(3,2,3)
    plot(ori_rad_mat(iOri,:), ori_hit_mat(iOri,:),['-o' colmat(iOri-1,:)])
    hold on
    subplot(3,2,4)
    plot(ori_rad_mat(iOri,:), ori_n_mat(iOri,:),['-o' colmat(iOri-1,:)])
    hold on;
end
subplot(3,2,3)
ylim([0 1])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')
title('Visual')
legend(num2str(chop(Oris(2:end)',2)))
subplot(3,2,4)
ylim([0 max(max(ori_n_mat,[],1),[],2)])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Number of trials')

tVolume = chop(double(celleqel2mat_padded(input.tSoundTargetAmplitude)),2);
Vols = unique(tVolume);
nVol = length(Vols);

vol_rad_mat = zeros(nVol,length(edges));
vol_hit_mat = zeros(nVol,length(edges));
vol_n_mat = zeros(nVol,length(edges));
for iVol = 2:nVol
    vol_ind = find(tVolume(sm2Ix) == Vols(iVol));
    for  i = 1:length(edges)
        ind2 = intersect(vol_ind, find(bin_b2 == i));
        vol_rad_mat(iVol,i) = nanmean(nanmean(rad_mat_down_norm(1:15,sm2Ix(ind2)),1),2);
        vol_hit_mat(iVol,i) = sum(successIx(sm2Ix(ind2)),2)/(sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        x = (sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        if isempty(x)
            x = 0;
        end
        vol_n_mat(iVol,i) = x;
    end
    subplot(3,2,5)
    plot(vol_rad_mat(iVol,:), vol_hit_mat(iVol,:),['-o' colmat(iVol-1,:)])
    hold on
    subplot(3,2,6)
    plot(vol_rad_mat(iVol,:), vol_n_mat(iVol,:),['-o' colmat(iVol-1,:)])
    hold on;
end
subplot(3,2,5)
ylim([0 1])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')
title('Auditory')
legend(num2str(chop(Vols(2:end)',2)))
subplot(3,2,6)
ylim([0 max(max(vol_n_mat,[],1),[],2)])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Number of trials')

subplot(3,2,2)
for iVol = 2:nVol
    vol_ind = find(tVolume == Vols(iVol));
    plot(Vols(iVol),sum(successIx(vol_ind),2)/(sum(successIx(vol_ind),2) + sum(missedIx(vol_ind),2)), 'og');
    hold on
end
for iOri = 2:nOri
    ori_ind = find(tGratingDirection == Oris(iOri));
    plot(Oris(iOri)./Oris(end),sum(successIx(ori_ind),2)/(sum(successIx(ori_ind),2) + sum(missedIx(ori_ind),2)), 'ok');
    hold on
end
ylim([0 1])
xlim([0 1])
xlabel('Change (%)')
ylabel('Hit rate')
suptitle([mouse ' ' date ' Pupil size (before press) and Hit rate'])
print([fnout '_prepress_pupilvHR.pdf'], '-dpdf');

%hit rate vs pupil size- before target
[n, edges] = histcounts([nanmean(rad_mat_target_norm(1:15,:),1)],5);
smIx = strcmp(input.trialOutcomeCell,'success') + strcmp(input.trialOutcomeCell,'ignore');
successIx = strcmp(input.trialOutcomeCell,'success');
missedIx = strcmp(input.trialOutcomeCell,'ignore');
sm1Ix = intersect(find(~isnan(rad_mat_target_norm(1,:))), intersect(b1Ix, find(smIx)));
sm2Ix = intersect(find(~isnan(rad_mat_target_norm(1,:))), intersect(b2Ix, find(smIx)));
[n_b1, bin_b1] = histc(nanmean(rad_mat_target_norm(1:15,sm1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_target_norm(1:15,sm2Ix),1),edges);
figure;
subplot(3,2,1)
for i = 1:length(edges)
    ind1 = find(bin_b1 == i);
    if (sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2))>5
        [pct1, err1] = binofit(sum(successIx(sm1Ix(ind1)),2), sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        ploterr(nanmean(nanmean(rad_mat_target_norm(1:15,sm1Ix(ind1)),1),2), pct1, std(nanmean(rad_mat_target_norm(1:15,sm1Ix(ind1)),1),[],2)./sqrt(length(ind1)), pct1-err1(1), 'ok')
        hold on;
    end
    ind2 = find(bin_b2 == i);
    if sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2)>5
        [pct2, err2] = binofit(sum(successIx(sm2Ix(ind2)),2), sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        ploterr(nanmean(nanmean(rad_mat_target_norm(1:15,sm2Ix(ind2)),1),2), pct2, std(nanmean(rad_mat_target_norm(1:15,sm2Ix(ind2)),1),[],2)./sqrt(length(ind2)), pct2-err2(1), 'og')
        hold on
    end
end
xlim([0 1])
ylim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')

tGratingDirection = cell2mat(input.tGratingDirectionDeg);
Oris = unique(tGratingDirection);
nOri = length(Oris);
colmat = strvcat('k', 'b', 'g', 'y', 'r', 'm');
ori_rad_mat = zeros(nOri,length(edges));
ori_hit_mat = zeros(nOri,length(edges));
ori_n_mat = zeros(nOri,length(edges));
for iOri = 2:nOri
    ori_ind = find(tGratingDirection(sm1Ix) == Oris(iOri));
    for  i = 1:length(edges)
        ind1 = intersect(ori_ind, find(bin_b1 == i));
        ori_rad_mat(iOri,i) = nanmean(nanmean(rad_mat_target_norm(1:15,sm1Ix(ind1)),1),2);
        ori_hit_mat(iOri,i) = sum(successIx(sm1Ix(ind1)),2)/(sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        x = (sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        if isempty(x)
            x = 0;
        end
        ori_n_mat(iOri,i) = x;
    end
    subplot(3,2,3)
    plot(ori_rad_mat(iOri,:), ori_hit_mat(iOri,:),['-o' colmat(iOri-1,:)])
    hold on
    subplot(3,2,4)
    plot(ori_rad_mat(iOri,:), ori_n_mat(iOri,:),['-o' colmat(iOri-1,:)])
    hold on;
end
subplot(3,2,3)
ylim([0 1])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')
title('Visual')
legend(num2str(chop(Oris(2:end)',2)))
subplot(3,2,4)
ylim([0 max(max(ori_n_mat,[],1),[],2)])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Number of trials')

vol_rad_mat = zeros(nVol,length(edges));
vol_hit_mat = zeros(nVol,length(edges));
vol_n_mat = zeros(nVol,length(edges));
for iVol = 2:nVol
    vol_ind = find(tVolume(sm2Ix) == Vols(iVol));
    for  i = 1:length(edges)
        ind2 = intersect(vol_ind, find(bin_b2 == i));
        vol_rad_mat(iVol,i) = nanmean(nanmean(rad_mat_target_norm(1:15,sm2Ix(ind2)),1),2);
        vol_hit_mat(iVol,i) = sum(successIx(sm2Ix(ind2)),2)/(sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        x = (sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        if isempty(x)
            x = 0;
        end
        vol_n_mat(iVol,i) = x;
    end
    subplot(3,2,5)
    plot(vol_rad_mat(iVol,:), vol_hit_mat(iVol,:),['-o' colmat(iVol-1,:)])
    hold on
    subplot(3,2,6)
    plot(vol_rad_mat(iVol,:), vol_n_mat(iVol,:),['-o' colmat(iVol-1,:)])
    hold on;
end
subplot(3,2,5)
ylim([0 1])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')
title('Auditory')
legend(num2str(chop(Vols(2:end)',2)))
subplot(3,2,6)
ylim([0 max(max(vol_n_mat,[],1),[],2)])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Number of trials')

subplot(3,2,2)
for iVol = 2:nVol
    vol_ind = find(tVolume == Vols(iVol));
    plot(Vols(iVol),sum(successIx(vol_ind),2)/(sum(successIx(vol_ind),2) + sum(missedIx(vol_ind),2)), 'og');
    hold on
end
for iOri = 2:nOri
    ori_ind = find(tGratingDirection == Oris(iOri));
    plot(Oris(iOri)./Oris(end),sum(successIx(ori_ind),2)/(sum(successIx(ori_ind),2) + sum(missedIx(ori_ind),2)), 'ok');
    hold on
end
ylim([0 1])
xlim([0 1])
xlabel('Change (%)')
ylabel('Hit rate')
suptitle([mouse ' ' date ' Pupil size (before target) and Hit rate'])
print([fnout '_pretarget_pupilvHR.pdf'], '-dpdf');


%distribution of pupil size
rad_all = (sqrt(Area_temp./pi))*calib;
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

%% saving
save([fnout '_pupil.mat'], 'Area', 'Centroid', 'frame_rate', 'rad_mat_down','centroid_mat_down','rad_mat_target', 'centroid_mat_target' )