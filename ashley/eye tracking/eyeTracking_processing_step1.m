awFSAV_eye_naive100ms;
rc = behavConstsAV;
iexp = 6;


SubNum = expt(iexp).SubNum;
date = expt(iexp).date;
run = expt(iexp).runs(1,:);
time_mat = expt(iexp).time_mat;
mouse = expt(iexp).mouse;
eyeFolder = expt(iexp).folder;
frame_rate = expt(iexp).frame_rate;
time = time_mat(1,:);

fnout = fullfile('Z:\Analysis\', mouse,eyeFolder,date,[mouse '-' date '-' run '-eyedata-']);
        
load(['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' SubNum '-' date '-' time '.mat'])

load(fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\data',mouse,'two-photon imaging',date,run,[run '_000_000_eye.mat']))
load(fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\data',mouse,'two-photon imaging',date,run,[run '_000_000.mat']))

nframes = info.config.frames;
ntrials = input.trialSinceReset;

data = squeeze(data);      % the raw images...
xc = size(data,2)/2;       % image center
yc = size(data,1)/2;
W=40;

rad_range = expt(iexp).eyeradrange;
data = double(data(yc-W:yc+W,xc-W:xc+W,:));
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
    [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range);%,'Method','TwoStage','Sensitivity',0.9
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
Centroid = cell2mat({eye.Centroid}');
Area = cell2mat({eye.Area}');
Radius = cell2mat({eye.Radius}');

%% view image of selected pupil and non-selected pupil

subplotsN = 25;
subplotSize = ceil(sqrt(subplotsN));
x = find(isnan(Area));
if length(x)>subplotsN
    minx = subplotsN;
else
    minx = length(x);
end
frames = sort(randsample(length(x),minx));
figure; 
suptitle(['example from ' num2str(length(x)) ' nan frames'])
start = 1;
for i = 1:minx
    subplot(subplotSize,subplotSize,start);
    imagesq(data(:,:,x(frames(i)))); 
    title(x(frames(i)))
    hold on
%         plot(Centroid_temp(x(frames(i)),1), Centroid_temp(x(frames(i)),2), 'ok', 'MarkerSize', 2*sqrt(Area_temp(x(frames(i)),1)/pi))
    start = start+1;
end
print([fnout '_nanframes.pdf'], '-dpdf');
y = setdiff(1:length(Area),x);
if length(y)>subplotsN
    miny = subplotsN;
else
    miny = length(y);
end
frames = sort(randsample(length(y),miny));
figure; 
suptitle('analyzed frames')
start = 1;
for i = 1:miny
    subplot(subplotSize,subplotSize,start);
    imagesq(data(:,:,y(frames(i)))); 
    title(y(frames(i)))
    hold on
    viscircles(Centroid(y(frames(i)),:),sqrt(Area(y(frames(i)))/pi),'EdgeColor','w');
    start = start+1;
end
print([fnout '_plotradcir_measuredframes.pdf'], '-dpdf');

% save first 100 frames as tif
writetiff(data(:,:,1:100),[fnout '_eye100frames'])