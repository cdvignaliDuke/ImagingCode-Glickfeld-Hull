%% extract eye data
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
calib = 1/26.6; %mm per pixel
pre_frames = ceil(1000./(1000/frame_rate));
post_frames = ceil(1000./(1000/frame_rate));

% Load and combine eye tracking data
    % Set current directory to crash folder
Area = {};
Centroid = {};
Eye_data = {};
for irun =  1:nrun
    %CD = [LG_base '\Data\2P_images\' mouse '\' date '_' mouse '\' ImgFolder(irun,:)];
    CD = ['\\crash.dhe.duke.edu\data\home\ashley\data\AW71\two-photon imaging\' date '\' ImgFolder(irun,:)];
    cd(CD);
    if irun == 2 & strcmp(date,'170923') & strcmp(mouse,'i698')
    fn = [ImgFolder(irun,:) '_000_001_eye.mat'];
    else
    fn = [ImgFolder(irun,:) '_000_000_eye.mat'];
    end
    data = load(fn);          % should be a '*_eye.mat' file

    data = squeeze(data.data);      % the raw images...
    xc = size(data,2)/2;       % image center
    yc = size(data,1)/2;
    W=40;

    rad_range = [12 22];
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
    Centroid{irun} = cell2mat({eye.Centroid}');
    Area{irun} = cell2mat({eye.Area}');
    Eye_data{irun} = data;
end
    
%% align eye data
Centroid_temp = [];
Area_temp = [];
Eye_data_temp = [];
for irun = 1:nrun
    %CD = [LG_base '\Data\2P_images\' mouse '\' date '_' mouse '\' ImgFolder(irun,:)];
    CD = ['\\crash.dhe.duke.edu\data\home\ashley\data\AW71\two-photon imaging\' date '\' ImgFolder(irun,:)];
    cd(CD);
    if irun == 2 & strcmp(date,'170923') & strcmp(mouse,'i698')
    imgMatFile = [ImgFolder(irun,:) '_000_001.mat'];
    else
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    end
    load(imgMatFile);
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];
    
    c = Centroid{irun};
    Centroid_temp = [Centroid_temp; c(1:min(nframes),:)];
    a = Area{irun};
    Area_temp = [Area_temp; a(1:min(nframes),:)];
    e = Eye_data{irun};
    Eye_data_temp = cat(3,Eye_data_temp, e(:,:,1:min(nframes)));
end
clear c a e Eye_data Centroid Area

% no measurement frames
figure; 
hist(sqrt(Area_temp./pi));
figure;
x = find(isnan(Area_temp));
if length(x)>25
    minx = 25;
else
    minx = length(x);
end
start = 1;
frames = sort(randsample(length(x),minx));
for i = 1:minx
    subplot(5,5,start);
    imagesq(Eye_data_temp(:,:,x(frames(i)))); 
    title(x(frames(i)))
    start = start+1;
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eyeTCs.mat']),'Area_temp','Centroid_temp')
%% sort eye data by cycle
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
nOn = input.nFramesOn;
if iscell(nOn)
    nOn = unique(celleqel2mat_padded(nOn));
end

area_trial = nan(pre_frames+post_frames,maxCyc+1,nTrials);
area_target = nan(pre_frames+post_frames,nTrials);
targetFramesOff = nan(1,nTrials);
for itrial = 1:nTrials
    if ~isnan(cStart(itrial))
        for icyc = 1:nCyc(itrial)
            if icyc > 1
                cyc_add = ((icyc-1)*nOn)+sum(tFramesOff(itrial,1:icyc-1));
            else
                cyc_add = 0;
            end
            if cStart(itrial)+post_frames+cyc_add <= size(Area_temp,1)
                area_trial(:,icyc,itrial) = Area_temp(cStart(itrial)-pre_frames+cyc_add:cStart(itrial)+post_frames+cyc_add-1,:);
            else
                area_trial(:,icyc,itrial) = NaN(pre_frames+post_frames,1);
            end
            if SIx(itrial) || MIx(itrial)
                if icyc == nCyc(itrial)
                    targetFramesOff(1,itrial) = tFramesOff(itrial,icyc-1);
                    if cStart(itrial)+post_frames+cyc_add <= size(Area_temp,1)
                        area_target(:,itrial) = Area_temp(cStart(itrial)-pre_frames+cyc_add:cStart(itrial)+post_frames+cyc_add-1,:);
                    else
                        area_target(:,itrial) = NaN(pre_frames+post_frames,1);
                    end
                end
            end
        end
    else
        area_trial(:,icyc,itrial) = NaN(pre_frames+post_frames,1);
    end
end    

radius_temp = sqrt(Area_temp./pi)*calib;
radius_trial = sqrt(area_trial./pi).*calib;
radius_target = sqrt(area_target./pi).*calib;
radius_trial_maxnorm = radius_trial./max(radius_temp,[],1);
radius_target_maxnorm = radius_target./max(radius_temp,[],1)';

%% figures
%plot radius by outcome
figure;
tt = [1-pre_frames:post_frames].*(1000./frame_rate);
pre_win = ceil(250./(1000/frame_rate));
edges = 0:0.1:1;
subplot(2,3,1)
shadedErrorBar(tt,nanmean(radius_trial_maxnorm(:,1,MIx),3),nanstd(radius_trial_maxnorm(:,1,MIx),[],3)./sqrt(length(MIx)),'r');
hold on
shadedErrorBar(tt,nanmean(radius_trial_maxnorm(:,1,SIx),3),nanstd(radius_trial_maxnorm(:,1,SIx),[],3)./sqrt(length(SIx)),'k');
ylim([0 1])
ylabel('% of Max Radius')
xlabel('Time from trial start (ms)')
title('Pre-trial')
subplot(2,3,2)
histogram(squeeze(nanmean(radius_trial_maxnorm(pre_frames-pre_win:pre_frames-1,1,SIx),1)),edges,'EdgeColor','k','FaceColor','none')
hold on
histogram(squeeze(nanmean(radius_trial_maxnorm(pre_frames-pre_win:pre_frames-1,1,MIx),1)),edges,'EdgeColor','r','FaceColor','none')
xlabel('% of Max Radius')
ylabel('Trials')
xlim([0 1])
subplot(2,3,3)
errorbar(1, nanmean(nanmean(radius_trial_maxnorm(pre_frames-pre_win:pre_frames-1,1,MIx),1),3),nanstd(nanmean(radius_trial_maxnorm(pre_frames-pre_win:pre_frames-1,1,MIx),1),[],3)./sqrt(length(MIx)),'or');
hold on
errorbar(1, nanmean(nanmean(radius_trial_maxnorm(pre_frames-pre_win:pre_frames-1,1,SIx),1),3),nanstd(nanmean(radius_trial_maxnorm(pre_frames-pre_win:pre_frames-1,1,SIx),1),[],3)./sqrt(length(SIx)),'ok');
ylim([0 1])
ylabel('% of Max Radius')
[h p] = ttest2(nanmean(radius_trial_maxnorm(pre_frames-pre_win:pre_frames-1,1,SIx),1),nanmean(radius_trial_maxnorm(pre_frames-pre_win:pre_frames-1,1,MIx),1));
title(['p = ' num2str(chop(p,2))])
subplot(2,3,4)
shadedErrorBar(tt,nanmean(radius_target_maxnorm(:,MIx),2),nanstd(radius_target_maxnorm(:,MIx),[],2)./sqrt(length(MIx)),'r');
hold on
shadedErrorBar(tt,nanmean(radius_target_maxnorm(:,SIx),2),nanstd(radius_target_maxnorm(:,SIx),[],2)./sqrt(length(SIx)),'k');
ylim([0 1])
ylabel('% of Max Radius')
xlabel('Time from target on (s)')
title('Pre-target')
subplot(2,3,5)
histogram(squeeze(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,SIx),1)),edges,'EdgeColor','k','FaceColor','none')
hold on
histogram(squeeze(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,MIx),1)),edges,'EdgeColor','r','FaceColor','none')
xlabel('% of Max Radius')
ylabel('Trials')
xlim([0 1])
subplot(2,3,6)
errorbar(1, nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,MIx),1),2),nanstd(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,MIx),1),[],2)./sqrt(length(MIx)),'or');
hold on
errorbar(1, nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,SIx),1),2),nanstd(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,SIx),1),[],2)./sqrt(length(SIx)),'ok');
ylim([0 1])
ylabel('% of Max Radius')
[h p] = ttest2(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,SIx),1),nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,MIx),1));
title(['p = ' num2str(chop(p,2))])
suptitle([mouse ' ' date ': Normalized radius by outcome (Black- hit, Red- miss)'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_radByOutcome.pdf']),'-dpdf','-bestfit')


%plot radius by outcome and ISI for target
MIxOSIx = find(MIx+SIx);
figure
col_mat = strvcat('b','r','y');
radius_target_all = [];
radius_target_hit = [];
radius_target_miss = [];
avg_radius_target_all = zeros(noff,2);
avg_radius_target_hit = zeros(noff,2);
avg_radius_target_miss = zeros(noff,2);
subplot(3,3,1)
for ioff = 1:noff
    ind = intersect(MIxOSIx,find(targetFramesOff == offs(ioff)));
    n(ioff) = sum(~isnan(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1)));
    shadedErrorBar(tt,nanmean(radius_target_maxnorm(:,ind),2),nanstd(radius_target_maxnorm(:,ind),[],2)./sqrt(n(ioff)),col_mat(ioff));
    avg_radius_target_all(ioff,:) = [nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1),2) nanstd(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1),[],2)./sqrt(n(ioff))];
    radius_target_all = [radius_target_all [nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1); ioff.*ones(1,length(ind))]];
    hold on
end
ylim([0 1])
ylabel('% of Max Radius')
xlabel('Time from target on (s)')
title(['All outcome- ' num2str(n)])
subplot(3,3,2)
for ioff = 1:noff
    ind = intersect(MIxOSIx,find(targetFramesOff == offs(ioff)));
    histogram(squeeze(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1)),edges,'EdgeColor',col_mat(ioff),'FaceColor','none')
    hold on
end
xlim([0 1])
ylabel('Trials')
xlabel('% of Max Radius')
subplot(3,3,3)
errorbar(offs.*(1000/frame_rate), avg_radius_target_all(:,1),avg_radius_target_all(:,2),'-o');
p = anova1(radius_target_all(1,:), radius_target_all(2,:),'off');
ylim([0 1])
ylabel('% of Max Radius')
xlabel('ISI (ms)')
title(['p = ' num2str(chop(p,2))])
subplot(3,3,4)
for ioff = 1:noff
    ind = intersect(find(SIx),find(targetFramesOff == offs(ioff)));
    n(ioff) = sum(~isnan(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1)));
    shadedErrorBar(tt,nanmean(radius_target_maxnorm(:,ind),2),nanstd(radius_target_maxnorm(:,ind),[],2)./sqrt(n(ioff)),col_mat(ioff));
    avg_radius_target_hit(ioff,:) = [nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1),2) nanstd(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1),[],2)./sqrt(n(ioff))];
    radius_target_hit = [radius_target_hit [nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1); ioff.*ones(1,length(ind))]];
    hold on
end
ylim([0 1])
ylabel('% of Max Radius')
xlabel('Time from target on (s)')
title(['Success- ' num2str(n)])
subplot(3,3,5)
for ioff = 1:noff
    ind = intersect(find(SIx),find(targetFramesOff == offs(ioff)));
    histogram(squeeze(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1)),edges,'EdgeColor',col_mat(ioff),'FaceColor','none')
    hold on
end
xlim([0 1])
ylabel('Trials')
xlabel('% of Max Radius')
subplot(3,3,6)
errorbar(offs.*(1000/frame_rate), avg_radius_target_hit(:,1),avg_radius_target_hit(:,2),'-o');
p = anova1(radius_target_hit(1,:), radius_target_hit(2,:),'off');
ylim([0 1])
ylabel('% of Max Radius')
xlabel('ISI (ms)')
title(['p = ' num2str(chop(p,2))])
subplot(3,3,7)
for ioff = 1:noff
    ind = intersect(find(MIx),find(targetFramesOff == offs(ioff)));
    n(ioff) = sum(~isnan(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1)));
    shadedErrorBar(tt,nanmean(radius_target_maxnorm(:,ind),2),nanstd(radius_target_maxnorm(:,ind),[],2)./sqrt(n(ioff)),col_mat(ioff));
    avg_radius_target_miss(ioff,:) = [nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1),2) nanstd(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1),[],2)./sqrt(n(ioff))];
    radius_target_miss = [radius_target_miss [nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1); ioff.*ones(1,length(ind))]];
    hold on
end
ylim([0 1])
ylabel('% of Max Radius')
xlabel('Time from target on (s)')
title(['Miss- ' num2str(n)])
subplot(3,3,8)
for ioff = 1:noff
    ind = intersect(find(MIx),find(targetFramesOff == offs(ioff)));
    histogram(squeeze(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1)),edges,'EdgeColor',col_mat(ioff),'FaceColor','none')
    hold on
end
xlim([0 1])
ylabel('Trials')
xlabel('% of Max Radius')
subplot(3,3,9)
errorbar(offs.*(1000/frame_rate), avg_radius_target_miss(:,1),avg_radius_target_miss(:,2),'-o');
p = anova1(radius_target_miss(1,:), radius_target_miss(2,:),'off');
ylim([0 1])
ylabel('% of Max Radius')
xlabel('ISI (ms)')
title(['p = ' num2str(chop(p,2))])
suptitle([mouse ' ' date ': Normalized radius by ISI and outcome pre-target (Blue- 250, Red- 500, Yellow- 750)'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_radByISI_pretarget.pdf']),'-dpdf','-bestfit')

%plot radius by ISI for all cyc
radius = cell(1,noff);
for ioff = 1:noff
    for itrial = 1:nTrials
        for icyc = 2:nCyc(itrial)-1
            if tFramesOff(itrial,icyc) == offs(ioff)
                radius{1,ioff} = [radius{1,ioff} radius_trial_maxnorm(:,icyc+1,itrial)];
            end
        end
    end
end
figure;
subplot(1,2,1)
avg_radius = [];
avg_radius_all = [];
for ioff = 1:noff
    r = radius{ioff};
    n(ioff) = sum(~isnan(nanmean(r(pre_frames-pre_win:pre_frames-1,:),1)));
    shadedErrorBar(tt,nanmean(r,2), nanstd(r,[],2)./sqrt(n(ioff)),col_mat(ioff));
    avg_radius = [avg_radius; [nanmean(nanmean(r(pre_frames-pre_win:pre_frames-1,:),1),2) nanstd(nanmean(r(pre_frames-pre_win:pre_frames-1,:),1),[],2)./sqrt(n(ioff))]];
    avg_radius_all = [avg_radius_all [nanmean(r(pre_frames-pre_win:pre_frames-1,:),1); ioff.*ones(1,size(r,2))]];
    hold on
end
ylim([0 1])
ylabel('% of Max Radius')
xlabel('Time from baseline on (s)')
title(['n = ' num2str(n)])
subplot(1,2,2)
errorbar(offs.*(1000/frame_rate), avg_radius(:,1), avg_radius(:,2),'-o')
p = anova1(avg_radius_all(1,:), avg_radius_all(2,:),'off');
ylim([0 1])
ylabel('% of Max Radius')
xlabel('ISI (ms)')
title(['p = ' num2str(chop(p,2))])
suptitle([mouse ' ' date ': Normalized radius by ISI pre-base (Blue- 250, Red- 500, Yellow- 750)'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_radByISI_prebase.pdf']),'-dpdf','-bestfit')


%plot pupil radius by outcome and target
avg_radius_miss = zeros(nDelta,2);
avg_radius_hit = zeros(nDelta,2);
radius_bytarg = [];
for itarg = 1:nDelta
    indm = intersect(find(MIx), find(tGratingDir == deltas(itarg)));
    inds = intersect(find(SIx), find(tGratingDir == deltas(itarg)));
    avg_radius_miss(itarg,:) = [nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,indm),1),2) nanstd(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,indm),1),[],2)./sqrt(length(indm))]; 
    avg_radius_hit(itarg,:) = [nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,inds),1),2) nanstd(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,inds),1),[],2)./sqrt(length(inds))]; 
    radius_bytarg = [radius_bytarg [[nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,indm),1); zeros(1,length(indm)); itarg.*ones(1,length(indm))] [nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,inds),1); ones(1,length(inds)); itarg.*ones(1,length(inds))]]];
end
figure;
errorbar(deltas, avg_radius_miss(:,1), avg_radius_miss(:,2),'-or')
hold on
errorbar(deltas, avg_radius_hit(:,1), avg_radius_hit(:,2),'-ok')
ylim([0 1])
ylabel('% of Max Radius')
xlabel('Target orientation (deg)')
p = anovan(radius_bytarg(1,:), {radius_bytarg(2,:),radius_bytarg(3,:)},'display','off');
title(['p= ' num2str(chop(p',2))])
suptitle([mouse ' ' date ': Normalized radius by outcome and Target (Black- Hit, Red- Miss)'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_radByTargByOutcome.pdf']),'-dpdf','-bestfit')

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eyeData.mat']), 'radius_target_maxnorm', 'radius_trial_maxnorm', 'radius', 'pre_frames', 'post_frames', 'pre_win', 'avg_radius_miss', 'avg_radius_hit','targetOffFrames')


