% SubNum = '626';
% date = '160203';
% retTime = '1142';
% retFolder = '004';
% mouse = 'AW26';
% retFName = '004_000_000';
% dirFolder = '005';

% load MWorks file
CD = 'Y:\home\andrew\Behavior\Data';
mworks = ['data-' 'i' SubNum '-' date '-' retTime]; 
load (fullfile(CD,mworks));

down = 10;
nON = input.nScansOn/down;
nOFF = input.nScansOff/down;
nStim = input.gratingElevationStepN.*input.gratingAzimuthStepN;

try %loading saved, registered tif
   try
       data_reg = readtiff(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,retFolder,'Retinotopy_reg2dirtuning.tif'));
   catch
       data_reg = readtiff(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,retFolder,'Retinotopy_V1.tif'));
   end
catch %otherwise, load data, downsample, register img

% Set current directory to imaging data location
CD = ['Z:\data\' mouse '\two-photon imaging\' date '\' retFolder];
cd(CD);
imgMatFile = [retFName '.mat'];
load(imgMatFile);

nframes = info.config.frames;
tic
data = sbxread(retFName,0,nframes);
toc
% pmt = 1; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
data = squeeze(data);
   
 
down = 10;
nON = input.nScansOn/down;
nOFF = input.nScansOff/down;
nStim = input.gratingElevationStepN.*input.gratingAzimuthStepN;


data_down = stackGroupProject(data,down);
clear data

data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
data = data_sub;
clear data_sub

%get direction tuning registration image and get cells

% fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
% cd(fileDirMasks);
% load('regImg.mat');
% load('mask&TCDir.mat');
% clear data_TC

load(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'regImg.mat'))

[out data_reg] = stackRegister(data, data_avg);
clear data
end

% fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
% cd(fileDirMasks);
% load('mask&TCDir.mat');
% clear data_TC

%load mask from task dataset
load(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'taskDrivenPixelsMask_targetandanticipation.mat'));

%figure out new mask based on this dataset?
dataTC = stackGetTimeCourses(data_reg,Kmask_ant);

% dataTC = stackGetTimeCourses(data_reg,mask_cell);

%% get trials sorted by stim type
tAz = double(cell2mat(input.tGratingAzimuthDeg));
tEl = double(cell2mat(input.tGratingElevationDeg));
Az = unique(tAz);
El = unique(tEl);
nTrials = length(tAz);
%% get dF/F
if size(dataTC,1) == nTrials*(nON+nOFF)
    data = reshape(dataTC,nON+nOFF,nTrials,size(dataTC,2));
else
    dataTC = dataTC(1:nTrials*(nON+nOFF),:);
    data = reshape(dataTC,nON+nOFF,nTrials,size(dataTC,2));
end

F = mean(data(nOFF-4:nOFF,:,:),1);
dF = bsxfun(@minus,data,F);
dFoverF = bsxfun(@rdivide,dF,F);

%% find responsive cells
dFoverFmean = mean(dFoverF(:,:,:),3);
%% responsive cells average for all directions
gratingCoord = {};
meanResp2Coord = zeros(size(dFoverFmean,1),nStim);
start = 1;
for iA = 1:length(Az)
    for iE = 1:length(El)
        ind = intersect(find(tAz == Az(iA)), find(tEl == El(iE )));
        gratingCoord{start} = ['(' num2str(Az(iA)) ',' num2str(El(iE)) ')'];
        meanResp2Coord(:,start) = mean(dFoverFmean(:,ind),2);
        start = start +1;
    end
end

%% plot
colors = brewermap(double(nStim),'Spectral');
colors = num2cell(colors,2);
retFig = figure;
retFig = plot(meanResp2Coord,'LineWidth',3);
set(retFig, {'color'},colors);
legend(gratingCoord, 'Location', 'southeastoutside')
hold on
vline(nOFF,'k:')
xlabel('frames')
ylabel('dF/F')
title({'retinotopy driven resp'; [mouse '-' date]; 'anticipation mask used'})

%%
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'avgRespRet'), '-dpdf');