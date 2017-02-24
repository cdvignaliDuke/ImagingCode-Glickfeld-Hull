clear all
close all

%% load data
awFSAVdatasets_V1

for iexp = [16]%1:size(expt,2)
    SubNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    time = expt(iexp).dirtuning_time;
    ImgFolder = expt(iexp).dirtuning;
    fName = [ImgFolder '_000_000'];
    exptName = [mouse '-' date];
    
    [input, data] = Load_SBXdataPlusMWorksData(SubNum,date,time,mouse,ImgFolder,fName);
    
%% analysis folder
rc = behavConstsAV;
fnpath = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',date,ImgFolder);

try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fnpath;
    cd(filedir);
end

%% figure handles and settings
setFigParams4Print('portrait')
ROIcorrFig = figure;
imgFig = figure;
maskFig = figure;
corrPropsFig = figure;

%% Parameters
frameRateHz = expt(iexp).frame_rate;
nON = double(input.nScansOn);
nOFF = double(input.nScansOff);
nStim = double(input.gratingDirectionStepN);
nRep = (size(data,3))./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
frames2S = frameRateHz*2;
tOffInd = double(1:nON+nOFF:(size(data,3)));
tOnInd = double(1:nON+nOFF:(size(data,3)))+nOFF;
itiFramesInd = linspaceNDim(tOnInd-frames2S,tOnInd-1,frames2S)';

%% select out only iti frames
data_iti = data(:,:,itiFramesInd(:));
%% register all iti data

%remove negative data
data_iti_sub = subNegData(data_iti);
clear data_iti

% register images to avgerage frame used across all analyses for this
% experiment
load([fnpath '\regImg.mat']);   

% first check if you've already done this step before...
if exist([fnpath '\itiFrames.tif'])
    data_reg = double(readtiff([fnpath '\itiFrames.tif']));
    data_stim_reg = double(readtiff([fnpath '\stimFrames.tif']));
    if size(data_reg,3) ~= length(itiFramesInd(:))
        [out data_reg] = stackRegister(data_iti_sub, data_avg);
        clear data_iti_sub

        writetiff(data_reg,'itiFrames');
    end       
else
    [out data_reg] = stackRegister(data_iti_sub, data_avg);
    clear data_iti_sub

    writetiff(data_reg,'itiFrames');
end

data_reg = double(data_reg);

% does the registered data look blurry?
figure;imagesc(mean(data_reg,3));colormap gray

%% get mask from ori tuning responses, then time-courses

if ~exist([fnpath '\mask&TCDir.mat'])
    error('need to get ROIs from direction tuning data')
else
    load([fnpath '\mask&TCDir.mat'])
    load([fnpath '\neuropil.mat'])
end

%get time-course
data_TC = stackGetTimeCourses(data_reg,mask_cell);
%% neuropil subtraction
nCells = size(data_TC,2);

buf = 4;
np = 6;

[npTC_iti_weighted, data_TCsubNP] = npSub(data_TC,mask_cell,data_reg,buf,np);
data_TC = data_TCsubNP;

%% thresholded (event extracted) time-courses
ROIcorrelations_threshold

%% Change parameters for downsampling
down = 10;
nON = double(input.nScansOn)./down;
nOFF = double(input.nScansOff)./down;
nStim = double(input.gratingDirectionStepN);
nRep = (size(data,3)./down)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
tOffInd = double(1:nON+nOFF:(size(data,3)./down));
tOnInd = double(1:nON+nOFF:(size(data,3)./down))+nOFF;

stimFramesInd = linspaceNDim(tOnInd,tOnInd+(frames2S/down)-1,(frames2S/down))';

%% Get and register stim frames

data_down = stackGroupProject(data,down);
data_stim = data_down(:,:,stimFramesInd(:));

data_stim_sub = subNegData(data_stim);
clear data_stim

% if exist([fnpath '\itiFrames.tif'])
%     data_stim_reg = double(readtiff([fnpath '\stimFrames.tif']));
% else
    [out data_stim_reg] = stackRegister(data_stim_sub, data_avg);
    clear data_iti_sub
    writetiff(data_stim_reg,'stimFrames');
% end
data_stim_reg = double(data_stim_reg);

% does the registered data look blurry?
figure;imagesc(mean(data_stim_reg,3));colormap gray

%get time-course
data_stim_TC = stackGetTimeCourses(data_stim_reg,mask_cell);

% subtract neuropil
[npTC_stim_weighted, data_TCsubNP neuropil_stim] = npSub(data_stim_TC,mask_cell,data_stim_reg,buf,np);
data_stim_TC = data_TCsubNP;
%% down-sample iti frames for rest of analysis
if divisible(frames2S,down)
    data_iti_reg = stackGroupProject(data_reg,down);
else
   error('fix this later') 
end
data_TC = stackGetTimeCourses(data_iti_reg,mask_cell);
[npTC_stim_weighted, data_TCsubNP neuropil_iti] = npSub(data_TC,mask_cell,data_iti_reg,buf,np);
data_TC = data_TCsubNP;

%% dF/F, global F
F = mean(data_TC,1);
dF = bsxfun(@minus,data_TC,F);
dFoverF = bsxfun(@rdivide,dF,F);

F = mean(data_stim_TC,1);
dF = bsxfun(@minus,data_stim_TC,F);
dFoverF_stim = bsxfun(@rdivide,dF,F);

%% remove motion trials
data_TC_trials = reshape(dFoverF,frames2S/down,nTrials,size(dFoverF,2));
data_stim_TC_trials = reshape(dFoverF_stim,frames2S/down,nTrials,size(dFoverF_stim,2));
iti_diff = max(diff(abs(mean(data_TC_trials,3)),1));
stim_diff = max(diff(abs(mean(data_stim_TC_trials,3)),1));

Fig1 = figure;
subplot(1,2,1)
hist(iti_diff)
title('iti motion')
subplot(1,2,2)
hist(stim_diff)
title('stim motion')
Fig1.PaperSize = [11 8.5];
print([fnpath '\itiMotionHist.pdf'],'-dpdf','-fillpage')

if expt(iexp).motionTestBimodal
    motionThreshold = expt(iexp).motionTestBimodal;
    ind_motion_iti = find(iti_diff>motionThreshold);
    ind_motion_stim = find(stim_diff>motionThreshold);
    
    dFoverF_temp = reshape(dFoverF,frames2S/down,nTrials,size(dFoverF,2));
    dFoverF_temp = reshape(dFoverF_temp(:,setdiff(1:nTrials,ind_motion_iti),:),(frames2S/down)*(nTrials-length(ind_motion_iti)),size(dFoverF,2));
    dFoverF = dFoverF_temp;
    
    dFoverF_temp = reshape(dFoverF_stim,frames2S/down,nTrials,size(dFoverF_stim,2));
    dFoverF_temp = reshape(dFoverF_temp(:,setdiff(1:nTrials,ind_motion_stim),:),(frames2S/down)*(nTrials-length(ind_motion_stim)),size(dFoverF_stim,2));
    dFoverF_stim = dFoverF_temp;
end


%% correlate time-courses
ROIcorrDown = corrcoef(dFoverF);
ROIstimCorrDown = corrcoef(dFoverF_stim);

ROIcorrDownChop = tril(ROIcorrDown,-1);
ROIstimCorrDownChop = tril(ROIstimCorrDown,-1);
%% make a histogram of correlations for iti and stim frames
itiCorrs = ROIcorrDownChop(ROIcorrDownChop ~=0);
stimCorrs = ROIstimCorrDownChop(ROIstimCorrDownChop ~=0);

bar_bins = -1:0.05:1;

iti_hist = hist(itiCorrs,bar_bins);
stim_hist = hist(stimCorrs,bar_bins);

figure(ROIcorrFig);
subplot(2,2,1)
bar(bar_bins,iti_hist)
title('iti frames correlations')
xlim([-1 1])
ylim([0 max([iti_hist stim_hist])])
xlabel('r')
ylabel('n pairs')
subplot(2,2,2)
bar(bar_bins,stim_hist)
title('stim frames correlations')
xlim([-1 1])
ylim([0 max([iti_hist stim_hist])])
xlabel('r')
ylabel('n pairs')
%% corr matrix fig
corrCutoffs = [0.25 0.5 0.75];
nPairs_str = cell(1,3);
nPairsStim_str = cell(1,3);
for iCut = 1:3
    nPairs_str{iCut} = [num2str(sum(sum(ROIcorrDownChop > corrCutoffs(iCut)))) '>' num2str(corrCutoffs(iCut))];
    nPairsStim_str{iCut} = [num2str(sum(sum(ROIstimCorrDownChop > corrCutoffs(iCut)))) '>' num2str(corrCutoffs(iCut))];
end
nPairs_str = [strjoin(nPairs_str,', ') ' (pairs>corr)'];
nPairsStim_str = [strjoin(nPairsStim_str,', ') ' (pairs>corr)'];


figure(ROIcorrFig);
colormap(brewermap([],'*RdBu'))
subplot(2,2,3)
imagesc(ROIcorrDown)
axis square
colorbar
caxis([-1 1])
title({'iti down-sampled'; nPairs_str})
subplot(2,2,4)
imagesc(ROIstimCorrDown)
axis square
colorbar
caxis([-1 1])
title({'stim frames down-sampled'; nPairsStim_str})

%% correlation between two halves of the same cell

sameCellCorr = zeros(1,nCells);
sameCellCorr_randPix = zeros(1,nCells);
for iCell = 1:nCells
    pix_ind = find(mask_cell == iCell);
    %halves
    half1_ind =  pix_ind(1:floor(length(pix_ind)/2));
    half2_ind = pix_ind(floor(length(pix_ind)/2)+1:end);
    
    mask1cell = zeros(size(mask_cell(:)));
    mask1cell(half1_ind) = 1;
    mask1cell = reshape(mask1cell,size(mask_cell));
    tc1 = stackGetTimeCourses(data_iti_reg,mask1cell);
    
    mask1cell = zeros(size(mask_cell(:)));
    mask1cell(half2_ind) = 1;
    mask1cell = reshape(mask1cell,size(mask_cell));
    tc2 = stackGetTimeCourses(data_iti_reg,mask1cell);
    
    c = corrcoef(tc1,tc2);
    sameCellCorr(iCell) = c(2,1);
    
    %rand pix
    half1_ind = randsample(pix_ind,floor(length(pix_ind)/2));
    half2_ind = setdiff(pix_ind, half1_ind);
    
    mask1cell = zeros(size(mask_cell(:)));
    mask1cell(half1_ind) = 1;
    mask1cell = reshape(mask1cell,size(mask_cell));
    tc1 = stackGetTimeCourses(data_iti_reg,mask1cell);
    
    mask1cell = zeros(size(mask_cell(:)));
    mask1cell(half2_ind) = 1;
    mask1cell = reshape(mask1cell,size(mask_cell));
    tc2 = stackGetTimeCourses(data_iti_reg,mask1cell);
    
    c = corrcoef(tc1,tc2);
    sameCellCorr_randPix(iCell) = c(2,1);
end
% compare rand pix vs. havles
figure(corrPropsFig);
subplot(3,2,1)
h = scatter(sameCellCorr,sameCellCorr_randPix,50,'k.');
hold on
plot(-1:0.1:1,-1:0.1:1,'k--')
xlim([-1 1])
ylim([-1 1])
axis square
xlabel('halves')
ylabel('rand pix')
title([exptName ', within cell correlations'])

% plot histogram of between halves vs between pairs with medians
figure(corrPropsFig);
subplot(3,2,2)
h = histogram(sameCellCorr,length(bar_bins)-1,'Normalization','probability');
hold on
h2 = histogram(itiCorrs,length(bar_bins)-1,'Normalization','probability');
h.EdgeColor = [1 1 1];
h.FaceColor = [0.25 0.25 1];
h2.FaceColor = 0.25*h.FaceColor;
h2.EdgeColor = [1 1 1];
hold on
l(1) = vline(median(sameCellCorr),'k--');
hold on 
l(2) = vline(median(itiCorrs),'r--');
ylabel('n (normalized)')
xlabel('correlation')
title('iti frames')
legend(l,{'median same cell'; 'median between cells'})
axis square



%% plot correlation vs distance bw cells
%get distance bw each ROI
ROIdist = ROIdistMat(mask_cell);
ROIdistChop = tril(ROIdist,-1);
ROIdistAll = ROIdistChop(ROIdistChop ~= 0);

%scatter r vs distance
figure(corrPropsFig)
subplot(3,2,3)
scatter(ROIdistAll,itiCorrs,50,'k.')
ylim([-1 1])
xlabel('distance between ROIs')
ylabel('correlation betweeen ROIs')
title('iti frames')
axis square

%% plot correlation vs. size
% get size of ROI in pixels
ROIsize = getROIsize(mask_cell);

ROIcorrNotSelf = reshape(ROIcorrDown(ROIcorrDown ~=1),nCells-1,nCells);
ROIcorrSortRows = sort(ROIcorrDown,1);
ROIcorrMaxMean = mean(ROIcorrSortRows(end-floor(nCells/20):end-1,:),1);
ROIcorrMax = ROIcorrSortRows(1,:);
ROIcorrMean = mean(ROIcorrNotSelf,1);

figure(corrPropsFig)
subplot(3,2,4)
scatter(ROIsize,ROIcorrMaxMean,50,'k.')
% ylim([-0.3 0.3])
xlabel('ROI size')
ylabel('mean corr with top 10% ROIs')
title('iti frames')
axis square


% plot within cell corr vs size
figure(corrPropsFig)
subplot(3,2,5)
scatter(ROIsize,sameCellCorr_randPix,50,'k.')
ylim([-1 1])
xlabel('ROI size')
ylabel('same cell corr')
title('iti frames')
axis square


%% plot dF/F vs size
ROIdfoverf = mean(dFoverF,1);

figure;
scatter(ROIsize,ROIdfoverf,50,'k.')
% ylim([-0.3 0.3])
xlabel('ROI size')
ylabel('mean dF/F')
title('iti frames')
axis square


%% test different correlation thresholds
for iCut = 1:3
%% find highly correlated pairs
    corrCutoff = corrCutoffs(iCut);

    nPairsDown = sum(sum(ROIcorrDownChop > corrCutoff));

    nPairs_stim = sum(sum(ROIstimCorrDownChop > corrCutoff));

    %only find pairs from down-sampled data
    [rowCell_iti,colCell_iti] = find(ROIcorrDownChop > corrCutoff);
    [rowCell_stim,colCell_stim] = find(ROIstimCorrDownChop > corrCutoff);

%% find if pair is correlated in stim frames too
    itiPairsInd = find(ROIcorrDownChop > corrCutoff);
    stimPairsInd = find(ROIstimCorrDownChop > corrCutoff);

    corrPairsInd = intersect(itiPairsInd,stimPairsInd);

    [rowCorrCell,colCorrCell] = ind2sub(size(ROIcorrDownChop),corrPairsInd);

%% ranking of correlations of cells that are in multiple pairs
sortEachCorrROI = sort(cat(1,rowCorrCell,colCorrCell));
multiMembers = unique(sortEachCorrROI(diff(sortEachCorrROI) == 0));

multiMemCorrs = cell(1,length(multiMembers));
for iCell = 1:length(multiMembers)
    roi = multiMembers(iCell);
   row_ind = rowCorrCell ==  roi;
   col_ind = colCorrCell == roi;
   pairs = cat(1,colCorrCell(row_ind),rowCorrCell(col_ind));
   multiMemCorrs{1,iCell} = ROIcorrDown(roi,pairs);
end
multiMax = cellfun(@max,multiMemCorrs);
multi2max = cellfun(@secondMax,multiMemCorrs);


    figure(maskFig);
    subplot(3,2,iCut*2)
    scatter(multiMax,multi2max,100,'k.')
    hold on
    plot(0:0.1:1,0:0.1:1,'k--')
    xlim([0 1])
    ylim([0 1])
    axis square
    xlabel('max corr')
    ylabel('2nd max corr')
    title('ROIs with multiple pairs')

%% plot time-course of individual cells
    allCorrCells = unique([rowCorrCell colCorrCell]);

    %%one cells correlation with several others
    cellTCFig1 = figure;
    nPlot = 5;
    nMatchCells = 10;
    colors = brewermap(nMatchCells+2,'*Greys');
    if length(allCorrCells) < nPlot
        nPlot = length(allCorrCells);
    end
    randCells = allCorrCells(randsample(length(allCorrCells),nPlot));

    for iplot = 1:nPlot
    cellCorrs = ROIcorrDown(randCells(iplot),:);
    [cellCorrsSort cellCorrsSortInd] = sort(cellCorrs);
    cellCorrsSort = fliplr(cellCorrsSort);
    cellCorrsSortInd = fliplr(cellCorrsSortInd);

    cellPlotInd = floor(linspace(2,length(cellCorrs),nMatchCells));
    cellPlotInd = [1 cellPlotInd(1:end-1)];

    figure(cellTCFig1)
    subplot(1,5,iplot)
    set(gca, 'ColorOrder', colors,'NextPlot', 'replacechildren');
    f = tcOffsetPlot(1:size(dFoverF,1),dFoverF(:,cellCorrsSortInd(cellPlotInd)),1);
    f(1).LineWidth = 2;
    title(['cell ' num2str(randCells(iplot))])

    text(repmat(size(dFoverF,1)+20,1,nMatchCells),0:nMatchCells-1,strread(num2str(chop(cellCorrsSort(cellPlotInd),2)),'%s'));
    xlim([0 size(dFoverF,1)+100])
    end
    cellTCFig1.PaperSize = [11 8.5];
    suptitle([exptName '-' num2str(corrCutoff)])
    print([fnpath '\exCorrCellTC1_' num2str(corrCutoff) '.pdf'],'-dpdf','-fillpage')

%% map correlated ROIs
    %draw lines bw correlated ROIs
        % find cell coordinates in image space, by indexing the first 
    xInd = NaN(1,length(allCorrCells));
    yInd = NaN(1,length(allCorrCells));
    xyInd = NaN(2,length(allCorrCells));
    for iCell= 1:length(allCorrCells)
        [yInd(1,iCell), xInd(1,iCell)] = find(mask_cell == allCorrCells(iCell),1);
    end

    xPair = NaN(2,1);
    yPair = NaN(2,1);
    figure(maskFig);
    subplot(3,2,iCut+(iCut-1))
    for iCell = 1:length(rowCorrCell)
        cell1 = rowCorrCell(iCell);
        cell2 = colCorrCell(iCell);

        xPair(1,1) = xInd(find(allCorrCells==cell1));
        yPair(1,1) = yInd(find(allCorrCells==cell1));

        xPair(2,1) = xInd(find(allCorrCells==cell2));
        yPair(2,1) = yInd(find(allCorrCells==cell2));

        hold on
        plot(xPair,yPair,'ko-')
    end
    xlim([1 size(mask_cell,2)])
    ylim([1 size(mask_cell,1)])

    ax = gca;
    ax.YDir = 'reverse';
    title('correlated cells connected')

%% combine ROIs
    allCorrCells = unique(cat(1,colCorrCell,rowCorrCell));
    comboROIs = [];
    for iCell = 1:length(allCorrCells)
        colPairs = find(colCorrCell == allCorrCells(iCell));
        rowPairs = find(rowCorrCell == allCorrCells(iCell));

        cGrp = cat(1,allCorrCells(iCell),rowCorrCell(colPairs));
        rGrp = cat(1,cGrp,colCorrCell(rowPairs));

        if iCell > 1
            grp_ind = find(cell2mat(cellfun(@(x) any(ismember(x,rGrp)),comboROIs,'unif',0)));
            if length(grp_ind) > 1
                cmbGrp = comboROIs(grp_ind);
                comboROIs{grp_ind(1)} = unique(cat(1,cell2mat(cmbGrp(:)),rGrp));
                comboROIs(grp_ind(2:end)) = {NaN};
                comboROIs{iCell} = NaN;
            elseif ~isempty(grp_ind)
                comboROIs{grp_ind} = unique(cat(1,comboROIs{grp_ind},rGrp));
                comboROIs{iCell} = NaN;
            else
                comboROIs{iCell} = rGrp;
            end

        else
            comboROIs{iCell} = rGrp;
        end

    end

    grp_ind = find(cell2mat(cellfun(@(x) sum(isnan(x)),comboROIs,'unif',0)) == 0);
    comboROIs = comboROIs(grp_ind);

    %change mask to be represented by smallest number cell in the combo ROI
    maskLinear = mask_cell(:);
    mask_cell_combo = maskLinear;
    neuropil_combo = zeros(size(neuropil_iti));
    discardFrames = [];
    for iCombo = 1:length(comboROIs)
        cellInd = find(ismember(maskLinear,comboROIs{iCombo}));
        mask_cell_combo(cellInd) = comboROIs{iCombo}(1);

        npNewFrame = sum(neuropil_iti(:,:,comboROIs{iCombo}),3);
        if length(unique(npNewFrame(:))) > 2
            npNewFrame(npNewFrame>1) = 0;
        end
        neuropil_combo(:,:,comboROIs{iCombo}(1)) = npNewFrame;
        discardFrames = cat(1,discardFrames,comboROIs{iCombo}(2:end));
    end

    mask_cell_combo = reshape(mask_cell_combo,size(mask_cell));
    neuropil_combo = neuropil_combo(:,:,setdiff(1:size(neuropil,3),unique(discardFrames)));
%% make spatial map of iti correlated pairs - each pair is a color
    cellImgLinear = zeros((size(mask_cell,1)*size(mask_cell,2)),3);

    % allCellsIndLinear = maskLinear > 0;

    %get linear index for each cell body
    corrCellColors = distinguishable_colors(length(comboROIs));

    for iCell = 1:length(comboROIs)
        grp_ind = find(ismember(maskLinear,comboROIs{iCell}));
        for iColor = 1:3
           cellImgLinear(grp_ind,iColor) = corrCellColors(iCell,iColor); 
        end    

    end
    cellImg = reshape(cellImgLinear,size(mask_cell,1),size(mask_cell,2),3);

%     figure(maskFig);
%     subplot(3,2,iCut*2)
%     image(cellImg)
%     title('correlated cells')
    
    imwrite(cellImg,fullfile(fnpath,[exptName '_corrMask_' num2str(corrCutoff) '.tif']));
end
%% save
save('maskROIcombo','mask_cell_combo','neuropil_combo','rowCorrCell','colCorrCell','corrCutoff')

figure(corrPropsFig)
suptitle(exptName)
print([fnpath '\' exptName 'corrPropsFig.pdf'],'-dpdf','-fillpage')
figure(ROIcorrFig)
suptitle(exptName)
print([fnpath '\' exptName 'ROIcorrFig.pdf'],'-dpdf','-fillpage')
% set params for figures
imgFig.PaperSize = [11 8.5];
figure(imgFig)
suptitle(exptName)
print([fnpath '\' exptName 'maxDFandMaskFig.pdf'],'-dpdf','-fillpage')
maskFig.PaperSize = [11 8.5];
figure(maskFig)
suptitle([exptName ' - corr cutoff = ' num2str(corrCutoff)])
print([fnpath '\' exptName 'ROImaskFig.pdf'],'-dpdf','-fillpage')

end