%clear all
close all
AWEyeDatasets
calib = 1/26.6; %mm per pixel
pre_event_time = 1000;
post_event_time = 4000;

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
        fn_mworks = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
        if irun == 1
            input = mwLoadData(fn_mworks, [], []);
        else
            input = [input mwLoadData(fn_mworks, [], [])];
        end
    end
    input = concatenateDataBlocks(input);

    runstr = catRunName(ImgFolder, nrun);
    fnout = ['Z:\home\lindsey\Analysis\Behavior\EyeTracking\' mouse '-' date '\' mouse '-' date '-' runstr];
    mkdir(fnout)
    %%
    preevent_frames = ceil(pre_event_time*(frame_rate/1000));
    postevent_frames = ceil(post_event_time*(frame_rate/1000));

    %% Load and combine eye tracking data
    % Set current directory to crash folder
    Area = {};
    Centroid = {};
    Eye_data = {};
    for irun =  1:nrun
        %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\' eyeFolder '\' date '\' runs(irun,:)];
        CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
        cd(CD);
        fn = [ImgFolder(irun,:) '_000_000_eye.mat'];
        data = load(fn);          % should be a '*_eye.mat' file

        data = squeeze(data.data);      % the raw images...
        xc = size(data,2)/2;       % image center
        yc = size(data,1)/2;
        W=40;

        rad_range = [6 15];
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

    %% reset frame counter
    run_trials = input.trialsSinceReset;
    
    cTrialStart = cell2mat(input.cTrialStart);
    cDecision = cell2mat(input.cDecision);
    
    Area_temp = [];
    Centroid_temp = [];
    Eye_data_temp = [];
    if nrun > 1
        for irun = 1:nrun
            if irun < nrun
                offset = size(Area{irun},1);
                startTrial = run_trials(irun)+1;
                endTrial = run_trials(irun)+run_trials(irun+1);
                cTrialStart(1,startTrial:endTrial) = cTrialStart(1,startTrial:endTrial)+offset;
                cDecision(1,startTrial:endTrial) = cDecision(1,startTrial:endTrial)+offset;
            end
            Area_temp = [Area_temp; Area{irun}];
            Centroid_temp = [Centroid_temp; Centroid{irun}];
            Eye_data_temp = cat(3, Eye_data_temp, Eye_data{irun});
        end
    else
        Area_temp = [Area_temp; Area{irun}];
        Centroid_temp = [Centroid_temp; Centroid{irun}];
        Eye_data_temp = cat(3, Eye_data_temp, Eye_data{irun});
    end
    clear Eye_data;
    ntrials = length(input.trialOutcomeCell);

%% no measurement frames
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
%         hold on
%         plot(Centroid_temp(x(frames(i)),1), Centroid_temp(x(frames(i)),2), 'ok', 'MarkerSize', 2*sqrt(Area_temp(x(frames(i)),1)/pi))
        start = start+1;
    end
    %print([fnout '_notnanframes.pdf'], '-dpdf');      
    
    %%
     nanrun = ceil(500*(frame_rate/1000));
    Rad_temp = sqrt(Area_temp./pi);
    sz = size(Eye_data_temp);
    rad_mat_start = zeros(preevent_frames+postevent_frames, ntrials);
    centroid_mat_start = zeros(preevent_frames+postevent_frames,2, ntrials);
    eye_mat_start = zeros(sz(1), sz(2), preevent_frames+postevent_frames, ntrials);
    rad_mat_decide = zeros(preevent_frames+postevent_frames, ntrials);
    centroid_mat_decide = zeros(preevent_frames+postevent_frames,2, ntrials);
    eye_mat_decide = zeros(sz(1), sz(2), preevent_frames+postevent_frames, ntrials);
    nframes = size(Rad_temp,1);
    for itrial = 1:ntrials
        if itrial == ntrials
            crange = [double(cTrialStart(itrial))-preevent_frames:nframes];
        else
            crange = [double(cTrialStart(itrial))-preevent_frames: double(cTrialStart(itrial+1)-preevent_frames-1)];
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
        if itrial < ntrials
            rad_mat_start(:,itrial) = Rad_temp(1+cTrialStart(itrial)-preevent_frames:cTrialStart(itrial)+postevent_frames,:);
            centroid_mat_start(:,:,itrial) = Centroid_temp(1+cTrialStart(itrial)-preevent_frames:cTrialStart(itrial)+postevent_frames,:);
            rad_mat_decide(:,itrial) = Rad_temp(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
            centroid_mat_decide(:,:,itrial) = Centroid_temp(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
            eye_mat_start(:,:,:,itrial) = Eye_data_temp(:,:,1+cTrialStart(itrial)-preevent_frames:cTrialStart(itrial)+postevent_frames);
            eye_mat_decide(:,:,:,itrial) = Eye_data_temp(:,:,1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames);
        else
            if (cTrialStart(itrial)+postevent_frames)<nframes
                rad_mat_start(:,itrial) = Rad_temp(1+cTrialStart(itrial)-preevent_frames:cTrialStart(itrial)+postevent_frames,:);
                centroid_mat_start(:,:,itrial) = Centroid_temp(1+cTrialStart(itrial)-preevent_frames:cTrialStart(itrial)+postevent_frames,:);
                eye_mat_start(:,:,:,itrial) = Eye_data_temp(:,:,1+cTrialStart(itrial)-preevent_frames:cTrialStart(itrial)+postevent_frames);
            else
                rad_mat_start(:,itrial) = nan(preevent_frames+postevent_frames,1);
                centroid_mat_start(:,:,itrial) = nan(preevent_frames+postevent_frames,2,1);
                eye_mat_start(:,:,:,itrial) = nan(sz(1),sz(2),preevent_frames+postevent_frames,1);
            end
            if (cDecision(itrial)+postevent_frames)<nframes
                rad_mat_decide(:,itrial) = Rad_temp(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
                centroid_mat_decide(:,:,itrial) = Centroid_temp(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
                eye_mat_decide(:,:,:,itrial) = Eye_data_temp(:,:,1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames);
            else
                rad_mat_decide(:,itrial) = nan(preevent_frames+postevent_frames,1);
                centroid_mat_decide(:,:,itrial) = nan(preevent_frames+postevent_frames,2,1);
                eye_mat_decide(:,:,:,itrial) = nan(sz(1),sz(2),preevent_frames+postevent_frames,1);
            end
        end
            
    end
    rad_mat_start = bsxfun(@times, rad_mat_start, calib);
    centroid_mat_start = bsxfun(@times,centroid_mat_start,calib);
    rad_mat_decide = bsxfun(@times, rad_mat_decide, calib);
    centroid_mat_decide = bsxfun(@times,centroid_mat_decide,calib);       
    %% Saving
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str  '_pupil.mat']), 'Area', 'Centroid', 'frame_rate' , 'rad_mat_start','centroid_mat_start','rad_mat_decide','centroid_mat_decide', 'input', 'cDecision', 'cTrialStart' );
end

