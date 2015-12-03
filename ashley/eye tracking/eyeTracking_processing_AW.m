clear all
close all
AWEyeDatasets_AW
calib = 1/26.6; %mm per pixel
min_hold = 2000;
pre_event_time = 1000;
post_release_time = 1500;
post_target_time = 4000;

for iexp = 1:size(expt,2)
    SubNum = expt(iexp).SubNum;
    date = expt(iexp).date;
    runs = expt(iexp).runs;
    time_mat = expt(iexp).time_mat;
    mouse = expt(iexp).mouse;
    eyeFolder = expt(iexp).folder;
    frame_rate = expt(iexp).frame_rate;
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
%     fnout = ['Z:\home\lindsey\Analysis\Behavior\EyeTracking\' mouse '-' date '\' mouse '-' date '-' runstr];
    fnout = ['Z:\Analysis\' mouse '\behavior\eye tracking\' date '\' mouse '-' date '-' runstr];

    %% 
    prepush_frames = ceil(pre_event_time*(frame_rate/1000));
    postpush_frames = ceil(min_hold*(frame_rate/1000));
    prerelease_frames = ceil(pre_event_time*(frame_rate/1000));
    postrelease_frames = ceil(post_release_time*(frame_rate/1000));
    pretarget_frames = ceil(pre_event_time*(frame_rate/1000));
    posttarget_frames = ceil(post_target_time*(frame_rate/1000));

    %% Load and combine eye tracking data
    % Set current directory to crash folder
    Area = {};
    Centroid = {};
    Eye_data = {};
    for irun =  1:nrun
        CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\' eyeFolder '\' date '\' runs(irun,:)];
        cd(CD);
        fn = [runs(irun,:) '_000_000_eye.mat'];
        load(fn);          % should be a '*_eye.mat' file

        data = squeeze(data);      % the raw images...
        xc = size(data,2)/2;       % image center
        yc = size(data,1)/2;
        W=40;

        rad_range = expt(iexp).eyeradrange;
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
            [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Method','TwoStage','Sensitivity',0.9);
            if(isempty(center))
                eye(n).Centroid = [NaN NaN];    % could not find anything...
                eye(n).Area = NaN;
            else
                [~,idx] = max(metric);          % pick the circle with best score
                eye(n).Centroid = center(idx,:);
                eye(n).Area = pi*radii(idx)^2;
                eye(n).Radius = radii;
            end
            if mod(n,100)==0
                fprintf('Frame %d/%d\n',n,size(data,3));
            end
        end
        Centroid{irun} = cell2mat({eye.Centroid}');
        Area{irun} = cell2mat({eye.Area}');
        Radius{irun} = cell2mat({eye.Radius}');
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
    Radius_temp = [];
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
        Radius_temp = [Radius_temp; Radius{irun}];
        Centroid_temp = [Centroid_temp; Centroid{irun}];
        Eye_data_temp = cat(3, Eye_data_temp, Eye_data{irun});
    end
    clear Eye_data;
    ntrials = length(input.trialOutcomeCell);

%% no measurement frames

    subplotsN = 25;
    subplotSize = ceil(sqrt(subplotsN));
    figure; 
    x = find(isnan(Area_temp));
    if length(x)>subplotsN
        minx = subplotsN;
    else
        minx = length(x);
    end
    start = 1;
    frames = sort(randsample(length(x),minx));
    for i = 1:minx
        subplot(subplotSize,subplotSize,start);
        imagesq(Eye_data_temp(:,:,x(frames(i)))); 
        title(x(frames(i)))
        hold on
%         plot(Centroid_temp(x(frames(i)),1), Centroid_temp(x(frames(i)),2), 'ok', 'MarkerSize', 2*sqrt(Area_temp(x(frames(i)),1)/pi))
        start = start+1;
    end
    print([fnout '_nanframes.pdf'], '-dpdf');


%% plot 25-100 random measured frames to check pupil looks good
    subplotsN = 25;
    subplotSize = ceil(sqrt(subplotsN));
    figure; 
    y = find(~isnan(Area_temp));
    if length(y)>subplotsN
        miny = subplotsN;
    else
        miny = length(y);
    end
    start = 1;
    frames = sort(randsample(length(y),miny));
    for i = 1:miny
        subplot(subplotSize,subplotSize,start);
        imagesq(Eye_data_temp(:,:,y(frames(i)))); 
        title(y(frames(i)))
        hold on
        viscircles(Centroid_temp(y(frames(i),:)),Radius_temp(y(frames(i))),'EdgeColor','w');
        viscircles(Centroid_temp(y(frames(i)),:),sqrt(Area_temp(y(frames(i)))/pi),'EdgeColor','w');
        plot(Centroid_temp(y(frames(i)),1), Centroid_temp(y(frames(i)),2), 'ok', 'MarkerSize', 2*sqrt(Area_temp(y(frames(i)),1)/pi))
        start = start+1;
    end
        print([fnout '_plotradcir.pdf'], '-dpdf');

     %% Remove NaNs if sparse and align to push, release and target
    nanrun = ceil(500*(frame_rate/1000));
    Rad_temp = sqrt(Area_temp./pi);
    rad_mat_down = zeros(prepush_frames+postpush_frames, ntrials);
    centroid_mat_down = zeros(prepush_frames+postpush_frames,2, ntrials);
    rad_mat_up = zeros(prerelease_frames+postrelease_frames, ntrials);
    centroid_mat_up = zeros(prerelease_frames+postrelease_frames,2, ntrials);
    rad_mat_target = zeros(pretarget_frames+posttarget_frames,ntrials);
    centroid_mat_target = zeros(pretarget_frames+posttarget_frames,2,ntrials);
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
    
    %% Saving
        save([fnout '_pupil.mat'], 'Area', 'Centroid', 'frame_rate', 'rad_mat_down','centroid_mat_down','rad_mat_target', 'centroid_mat_target','rad_mat_up','centroid_mat_up' );

end

