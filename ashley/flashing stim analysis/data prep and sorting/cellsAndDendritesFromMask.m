clear all
close all

%% load data & set path names
awFSAVdatasets_V1

iexp = 3;

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
date = expt(iexp).date;
time = expt(iexp).dirtuning_time;
ImgFolder = expt(iexp).dirtuning;
fName = [ImgFolder '_000_000'];
exptName = [mouse '-' date];

rc = behavConstsAV;
fnpath = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',date,ImgFolder);

mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (fullfile(rc.behavData,mworks));

data_reg = readtiff(fullfile(fnpath,'DirectionTuning_V1.tif'));

if exist(fullfile(fnpath,'maxDFoverF.tif'));
    maxDFoverF = readtiff(fullfile(fnpath,'maxDFoverF.tif'));
else
    if iscell(input.nScansOn)
        nON = unique(cell2mat(input.nScansOn))/10;
        nOFF = unique(cell2mat(input.nScansOff))/10;
    else
        nON = input.nScansOn/10;
        nOFF = input.nScansOff/10;
    end
    nTrials = input.trialSinceReset;
    tStartInd = double(1:nON+nOFF:size(data_reg,3));
    off_ind = linspaceNDim(tStartInd,tStartInd+(nOFF-1),nOFF);
    F = mean(data_reg(:,:,off_ind(:)),3);
    dF = bsxfun(@minus,data_reg,F);
    dFoverF = bsxfun(@rdivide,dF,F);
    maxDFoverF = max(dFoverF,[],3);
end

load(fullfile(fnpath,'mask&TCDir.mat'))
load(fullfile(fnpath,'neuropil.mat'))
tc = npSubTC;
%% select cells from mask_cell
nROI = length(unique(mask_cell))-1;
figure; imagesc(maxDFoverF); colormap gray; caxis([0 max(maxDFoverF(:))/4])

mask_ones = mask_cell;
mask_ones(mask_ones > 0) = 1;

bwout = imCellEditInteractive(mask_ones);
cellsOnly = bwlabel(bwout);
cellsOnlyInd = find(cellsOnly > 1);

cellsMatch = unique(mask_cell(cellsOnlyInd));
if any(cellsMatch == 0)
    cellsMatch = cellsMatch(2:end);
end
dendritesMatch = setdiff(1:nROI,cellsMatch);

%% separate time-courses

cell_tc = tc(:,cellsMatch);
den_tc = tc(:,dendritesMatch);

%% save
save(fullfile(fnpath,'cell&dendriteTC.mat'),'cell_tc','den_tc');
save(fullfile(fnpath,'cell&dendriteIndices.mat'),'cellsMatch','dendritesMatch');