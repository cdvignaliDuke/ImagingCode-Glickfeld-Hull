clear all
close all

ds = 'szTuning_PM';
doRet = true;

eval(ds)

doManualSearchRegion = true;
doFindPupil = true;

slct_exp = 4:5;
%%
calib = 1/26.6; %mm per pixel

rc = behavConstsAV;
% imgParams_FSAV

for iexp = slct_exp
    if contains(ds,'szTuning')
        if doRet
            runs = [expt(iexp).retinotopyFolder];
            times = [expt(iexp).retinotopyTime];
        else
            runs = [expt(iexp).sizeTuningFolder];
            times = [expt(iexp).sizeTuningTime];
        end
    end

    expDate = expt(iexp).date;
    mouse = expt(iexp).mouse;
    subnum = mouse;
    frame_rate = params.frameRate;
    nrun = length(runs);
%% load and combine mworks files
if ~any(isnan(expt(iexp).eyeradrange))
    for irun = 1:nrun
        runFolder = runs{irun};
        runTime = times{irun};
        fName = [runFolder '_000_000'];
        fn = fullfile(rc.ashleyData,mouse,'two-photon imaging', expDate,runFolder);
        
        cd(fn)
        try
            eyeName = [runFolder '_000_000_eye.avi'];
            eyeObj = VideoReader(eyeName);
            data = read(eyeObj,[1 Inf]);
            data = squeeze(data(:,:,1,:));
        catch
            fprintf('loading...')
%             try
                eyeName = [runFolder '_000_000_eye.mat'];
                load(eyeName)
                data = squeeze(data(:,:,1,:));
%             catch
%                 continue
%             end
        end
        
        fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,runFolder);
        if ~exist(fnout,'dir') ~= 7
            mkdir(fnout)
        end
    %% Load and combine eye tracking data
    % Set current directory to crash folder
    Area = {};
    Centroid = {};
    Eye_data = {};
%     for irun =  1:nrun
% %         CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\' eyeFolder '\' date '\' runs(irun,:)];
%         CD = fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\data\',mouse,'two-photon imaging',expDate,runs(irun,:));
%         cd(CD);
%         fn = [runs(irun,:) '_000_000_eye.mat'];
%         load(fn);          % should be a '*_eye.mat' file

%         data = squeeze(data);      % the raw images...
        xc = size(data,2)/2;       % image center
        yc = size(data,1)/2;
        W=40;

        rad_range = expt(iexp).eyeradrange;
        if strcmp(subnum,'668')
            data = double(data(yc-W:yc+W,20:(W*2+20),:));
        else
            data = double(data(yc-W:yc+W,xc-W:xc+W,:));
        end
        warning off;
        if doManualSearchRegion
            region = imCellPolyEditInteractive(data(:,:,1)); % draw ROI around eyeball to avoid false circles on fur
        end
        
        fprintf('saving movie...\n')
        eye_mov = data(:,:,[1:500:size(data,3)]);
        writetiff(eye_mov,fullfile(fnout,'eye_movie'))
        
        if doFindPupil
            A = cell(size(data,3),1);
            B = cell(size(data,3),1);
            for n = 1:size(data,3)
                A{n} = [0,0];
                B{n} = [0];
            end
            eye = struct('Centroid',A,'Area',B);
            radii = [];
            for n = 1:size(data,3)
                [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range);%,'Method','TwoStage','Sensitivity',0.9
                if(isempty(center))
                    eye(n).Centroid = [NaN NaN];    % could not find anything...
                    eye(n).Area = NaN;
                    eye(n).Radius = NaN;
                else
                    if doManualSearchRegion
                        [~,metricOrder] = sort(metric,'descend'); % pick the circle with best score
                        isInRegion = false(length(metric),1); % check that it is within the eyeball
                        c = round(center);
                        for ir = 1:length(metric) 
                            isInRegion(ir) = region(c(ir,2),c(ir,1));
                        end
                        if any(isInRegion)
                            idx = metricOrder(find(ismember(metricOrder,find(isInRegion==1)),1));
                        else
                            idx = nan;
                        end
                    end
                    if isnan(idx)
                        eye(n).Centroid = [NaN NaN];    % could not find anything...
                        eye(n).Area = NaN;
                        eye(n).Radius = NaN;
                    else
                        eye(n).Centroid = center(idx,:);
                        eye(n).Area = pi*radii(idx)^2;
                        eye(n).Radius = radii;
                    end
                end
                if mod(n,100)==0
                    fprintf('Frame %d/%d\n',n,size(data,3));
                end
            end

            Centroid = cell2mat({eye.Centroid}');
            Area = cell2mat({eye.Area}');
            Radius = cell2mat({eye.Radius}');

            % plot missed frames
            subplotsN = 25;
            subplotSize = ceil(sqrt(subplotsN));
            figure; 
            x = find(isnan(Area));
            suptitle(sprintf('%s skipped frames',num2str(length(x))))
            if length(x)>subplotsN
                minx = subplotsN;
            else
                minx = length(x);
            end
            start = 1;
            frames = sort(randsample(length(x),minx));
            for i = 1:minx
                subplot(subplotSize,subplotSize,start);
                imagesq(data(:,:,x(frames(i)))); 
                title(x(frames(i)))
                hold on
                start = start+1;
            end
            print(fullfile(fnout,'skippedEyeFrames'),'-dpdf','-fillpage')

            % plot tracked frames
            subplotsN = 25;
            subplotSize = ceil(sqrt(subplotsN));
            figure; 
            y = find(~isnan(Area));
            if length(y)>subplotsN
                miny = subplotsN;
            else
                miny = length(y);
            end
            start = 1;
            frames = sort(randsample(length(y),miny));
            for i = 1:miny
                subplot(subplotSize,subplotSize,start);
                imagesq(data(:,:,y(frames(i)))); 
                title(y(frames(i)))
                hold on
                viscircles(Centroid(y(frames(i)),:),sqrt(Area(y(frames(i)))/pi),'EdgeColor','w');
                start = start+1;
            end
            print(fullfile(fnout,'trackedEyeFrames'),'-dpdf','-fillpage')

            % save
            save(fullfile(fnout,'eyeTC.mat'),'Centroid','Area','Radius')
        end
    end
end
end


