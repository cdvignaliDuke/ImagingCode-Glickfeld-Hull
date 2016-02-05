SubNum = '626';
date = '160203';
time = '1151';
ImgFolder = '005';
mouse = 'AW26';
fName = '005_000_000';
dirFolder = '005';

% load MWorks file
CD = 'Y:\home\andrew\Behavior\Data';
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (fullfile(CD,mworks));

down = 10;
nON = input.nScansOn/down;
nOFF = input.nScansOff/down;
nStim = input.gratingDirectionStepN;

data_reg = readtiff(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,ImgFolder,'DirectionTuning_V1.tif'));

fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
cd(fileDirMasks);
load('mask&TCDir.mat');
clear data_TC

%figure out new mask based on this dataset?

dataTC = stackGetTimeCourses(data_reg,mask_cell);

%% get trials sorted by stim type
tDir = double(cell2mat(input.tGratingDirectionDeg));
nTrials = length(tDir);
dir = unique(tDir);
%% get dF/F
data = reshape(dataTC,nON+nOFF,nTrials,size(dataTC,2));
F = mean(data(nOFF-4:nOFF,:,:),1);
dF = bsxfun(@minus,data,F);
dFoverF = bsxfun(@rdivide,dF,F);

%% responsive cells average for each direction
gratingDir = {};
resp2Dir = zeros(size(dFoverF,1),size(dataTC,2),nStim);
for iD = 1:nStim
    ind = find(tDir == dir(iD));
    gratingDir{iD} = [num2str(dir(iD)) sprintf('%c', char(176))];
    resp2Dir(:,:,iD) = squeeze(mean(dFoverF(:,ind,:),2));
end
%% find responsive cells
% respCutoff = 0.2;
% gratingRespCells = find(any(squeeze(mean(resp2Dir(nOFF+1:nON+nOFF,:,:),1)) > respCutoff,2));
[maxRespPerCell maxRespInd] = max(squeeze(mean(resp2Dir(nOFF+1:nON+nOFF,:,:),1)),[],2);
dFoverFmean = zeros(size(dFoverF,1),nStim);
for iD = 1:nStim 
   ind = find(maxRespInd == iD);
   dFoverFmean(:,iD) = squeeze(mean(resp2Dir(:,ind,iD),2));
end

%% plot
colors = brewermap(double(nStim),'Spectral');
colors = num2cell(colors,2);
dirFig = figure;
dirFig = plot(dFoverFmean,'LineWidth',3);
set(dirFig, {'color'},colors);
legend(gratingDir, 'Location', 'southeastoutside')
hold on
vline(nOFF,'k:')
xlabel('frames')
ylabel('dF/F')
title({[mouse '-' date]; ['nCells = ' num2str(size(dataTC,2))]})

%%
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'avgRespDir'), '-dpdf');