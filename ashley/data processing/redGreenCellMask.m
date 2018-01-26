clear all
close all
ds = 'FSAV_V1_GAD';
rc = behavConstsAV;
eval(ds)
slct_expt = 2:size(expt,2);
doMaxDFF = 1;
%%
iexp = 3

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;
redChannelOn = expt(iexp).greenredsimultaneous;
regImgInd = expt(iexp).regImgStartFrame;
down = 10;

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');
%%
if doMaxDFF == 0
%%
    load(fullfile(fnout,'corrImages.mat'))
    load(fullfile(fnout,'redChannelImage.mat'))

    [ypix,xpix] = size(redImage);

    corrImage = mean(cat(3,corrImage_bx,corrImage_tun),3);

    figure
    colormap gray
    imagesc(corrImage)

    %**enter vals here***
    xcrop = [1:15 786:xpix];
    ycrop = [1:10 261:ypix];

    corr_crop = corrImage;
    corr_crop(:,xcrop) = 0;
    corr_crop(ycrop,:) = 0;
    imagesc(corr_crop)

    red_crop = redImage;
    red_crop(:,xcrop) = 0;
    red_crop(ycrop,:) = 0;
    imagesc(red_crop)

    red_crop_norm = red_crop./(max(red_crop(:)));

    if expt(iexp).areaBorders    
        borders = readtiff(fullfile(fn,'FOVborders.tif'));
        redImage = mean(cat(3,red_crop_norm,borders),3);
        corrImage = mean(cat(3,corr_crop,borders),3);
    else
        redImage = red_crop_norm;
        corrImage = corr_crop;
    end

    figure
    suptitle([mouse '' expDate])
    colormap gray
    subplot 211
    imagesc(corrImage)
    title('behavior & tuning correlation image')
    subplot 212
    imagesc(redImage)
    title(sprintf('%s+ red channel',expt(iexp).redChannelLabel))

    %%
    redCellsClickLabel = imCellEditInteractive(redImage);
    activeCellsClickLabel = imCellEditInteractive(corrImage);

    redImageLeftovers = redImage;
    redBuffer = sum(imCellBuffer(bwlabel(redCellsClickLabel),2),3) + redCellsClickLabel;
    redImageLeftovers(redBuffer > 0) = 0;
    activeImageLeftovers = corrImage;
    activeBuffer = sum(imCellBuffer(bwlabel(corrImage),2),3) + activeCellsClickLabel;
    activeImageLeftovers(activeBuffer > 0) = 0;

    cellSelectQuads = {[1, floor(xpix/2), 1, floor(ypix/2)];...
        [floor(xpix/2)+1, xpix, 1, floor(ypix/2)];
        [1, floor(xpix/2), floor(ypix/2)+1, ypix];
        [floor(xpix/2)+1, xpix, floor(ypix/2)+1, ypix]};

    polySelectLabel = nan(ypix,xpix,4);
    for iquad = 1:4
        polySelectLabel(:,:,iquad) = imCellPolyEditInteractive(redImageLeftovers,[],cellSelectQuads{iquad});
    end
    redLeftoversLabel = sum(polySelectLabel,3);

    polySelectLabel = nan(ypix,xpix,4);
    for iquad = 1:4
        polySelectLabel(:,:,iquad) = imCellPolyEditInteractive(activeImageLeftovers,[],cellSelectQuads{iquad});
    end
    activeLeftoversLabel = sum(polySelectLabel,3);

    redCellsLabel = logical(redLeftoversLabel+redCellsClickLabel);
    activeCellsLabel = logical(activeLeftoversLabel+activeCellsClickLabel);

    %%
    redCellsMask = bwlabel(redCellsLabel);
    maskMeasurements = regionprops(redCellsLabel,'Centroid');
    redCellCenters = cell2mat({maskMeasurements.Centroid}');
    nRedCells = length(unique(redCellsMask(:)))-1;
    redCellTxtLabel = cellstr(num2str([1:nRedCells]'));
    redCellsLabel2 = double(redCellsLabel);
    redCellsLabel2(redCellsLabel2 > 0) = 2;

    activeCellsMask = bwlabel(activeCellsLabel);
    maskMeasurements = regionprops(activeCellsLabel,'Centroid');
    activeCellCenters = cell2mat({maskMeasurements.Centroid}');
    nActiveCells = length(unique(activeCellsMask(:)))-1;
    activeCellTxtLabel = cellstr(num2str([1:nActiveCells]'));
    activeCellsBuffers = imCellBuffer(activeCellsMask,1);
    activeBufferInd = sum(activeCellsBuffers,3) > 0;

    redAndActiveLabel = redCellsLabel2 + activeCellsLabel;

    figure
    imagesc(redAndActiveLabel)
    colorbar

    for icell = 1:nRedCells
        if icell == 1
            twoColorLabel = redAndActiveLabel > 0;
            twoColorID = zeros(ypix,xpix);
            redActiveCells = false(1,nRedCells);
        end
        activeCellPix = activeCellsMask(round(redCellCenters(icell,2)),round(redCellCenters(icell,1)));
        if activeCellPix > 0
            activeCell = activeCellsMask == activeCellPix;
            redCell = redCellsMask == icell;
            redEdge = ~activeCell & redCell;
            twoColorLabel(activeCell) = true;
            twoColorLabel(redEdge) = false;
            redActiveCells(icell) = true;
        elseif any(redCellsMask(:) == icell & activeCellsMask(:) > 0)
            activeBufferInd = sum(activeCellsBuffers,3) > 0;
            redCell = redCellsMask == icell;
            overlapInd = redCell & (activeCellsMask > 0 | activeBufferInd);
            twoColorLabel(overlapInd) = false;
        end
    end

    twoColorMask = bwlabel(twoColorLabel);
    figure;imagesc(twoColorMask)
    nCells = length(unique(twoColorMask(:)))-1;


    for icell = 1:nRedCells
        if icell == 1
            isLabeledCell = false(nCells,1);
            redCellFromActiveMask = false(nCells,1);
        end    
        cellPix = twoColorMask(round(redCellCenters(icell,2)),round(redCellCenters(icell,1)));
        if cellPix > 0
            isLabeledCell(cellPix) = true;
        end
        if redActiveCells(icell)
            redCellFromActiveMask(cellPix) = true;
        end
    end

    isFromActiveMask = true(nCells,1);
    isFromActiveMask(isLabeledCell & ~redCellFromActiveMask) = false;

    save(fullfile(fnout,'twoColorMask&CellLabels'),'twoColorMask',...
        'isLabeledCell','isFromActiveMask','redImage','corrImage')
elseif doMaxDFF == 1
    load(fullfile(fnout,'redChannelImage.mat'))
    load(fullfile(fnout,'bx_max_images.mat'))
    load(fullfile(fnout,'tun_max_images.mat'))
    
    [ypix,xpix] = size(redImage);

    bxDFFImage = mean(dFF_bxMax,3);
    tunDFFImage = mean(dFF_dirmax,3);

    figure
    colormap gray
    imagesc(bxDFFImage)

    %**enter vals here***
    xcrop = [1:15 794:xpix];
    ycrop = [1:10 261:ypix];

    mean_bx_crop = bxDFFImage;
    mean_bx_crop(:,xcrop) = 0;
    mean_bx_crop(ycrop,:) = 0;
    imagesc(mean_bx_crop)
    
    bx_crop = dFF_bxMax;
    bx_crop(:,xcrop) = 0;
    bx_crop(ycrop,:) = 0;
    imagesc(mean(bx_crop,3))
    
    tun_crop = dFF_dirmax;
    tun_crop(:,xcrop) = 0;
    tun_crop(ycrop,:) = 0;
    imagesc(mean(tun_crop,3))
    
    red_crop = redImage;
    red_crop(:,xcrop) = 0;
    red_crop(ycrop,:) = 0;
    imagesc(red_crop)

    red_crop_norm = red_crop./(max(red_crop(:)));

    if expt(iexp).areaBorders    
        borders = readtiff(fullfile(fn,'FOVborders.tif'));
        redImage = mean(cat(3,red_crop_norm,borders),3);
        bxImage = nan(size(bx_crop));
        tunImage = nan(size(tun_crop));
        for i = 1:size(tun_crop,3)
            if i <= 3
                bxImage(:,:,i) = mean(cat(3,bx_crop(:,:,i),borders),3);
            end
            tunImage(:,:,i) = mean(cat(3,tun_crop,borders),3);
        end
    else
        redImage = red_crop_norm;
        bxImage = bx_crop;
        tunImage = tun_crop;
    end

    figure
    suptitle([mouse '' expDate])
    colormap gray
    subplot 311
    imagesc(mean(bxImage,3))
    title('behavior max dF/F image')
    subplot 312
    imagesc(mean(tunImage,3))
    title('tuning max dF/F image')
    subplot 313
    imagesc(redImage)
    title(sprintf('%s+ red channel',expt(iexp).redChannelLabel))
    %%
    redCellsClickLabel = imCellEditInteractive(redImage);
    polySelectLabel = nan(ypix,xpix,4);
    redImageLeftovers = redImage;
    redBuffer = sum(imCellBuffer(bwlabel(redCellsClickLabel),2),3) + redCellsClickLabel;
    redImageLeftovers(redBuffer > 0) = 0;
    
    cellSelectQuads = {[1, floor(xpix/2), 1, floor(ypix/2)];...
        [floor(xpix/2)+1, xpix, 1, floor(ypix/2)];
        [1, floor(xpix/2), floor(ypix/2)+1, ypix];
        [floor(xpix/2)+1, xpix, floor(ypix/2)+1, ypix]};
    for iquad = 1:4
        polySelectLabel(:,:,iquad) = imCellPolyEditInteractive(redImageLeftovers,[],cellSelectQuads{iquad});
    end
    redPolyLabel = sum(polySelectLabel,3);
    redLabel = redCellsClickLabel+redPolyLabel;
    redCellsMask = bwlabel(redCellsClickLabel);
    maskMeasurements = regionprops(redCellsClickLabel,'Centroid');
    redCellCenters = cell2mat({maskMeasurements.Centroid}');
    redCellsLabel2 = double(redCellsClickLabel);
    redCellsLabel2(redCellsLabel2 > 0) = 2;
    nRedCells = length(unique(redCellsMask(:)))-1;
    bxCellsClickMask = maskFromMultiMaxDFFStack(cat(3,bxImage,tunImage));
    bxCellsClickLabel = bxCellsClickMask;
    bxCellsClickLabel(bxCellsClickLabel > 0) = 1;
    activeCellsBuffers = imCellBuffer(bxCellsClickMask,1);
    activeBufferInd = sum(activeCellsBuffers,3) > 0;
    
    redAndActiveLabel = redCellsLabel2 + bxCellsClickLabel;

    figure
    imagesc(redAndActiveLabel)
    colorbar

    for icell = 1:nRedCells
        if icell == 1
            twoColorLabel = redAndActiveLabel > 0;
            twoColorID = zeros(ypix,xpix);
            redActiveCells = false(1,nRedCells);
        end
        activeCellPix = bxCellsClickMask(round(redCellCenters(icell,2)),round(redCellCenters(icell,1)));
        if activeCellPix > 0
            activeCell = bxCellsClickMask == activeCellPix;
            redCell = redCellsMask == icell;
            redEdge = ~activeCell & redCell;
            twoColorLabel(activeCell) = true;
            twoColorLabel(redEdge) = false;
            redActiveCells(icell) = true;
        elseif any(redCellsMask(:) == icell & bxCellsClickMask(:) > 0)
            activeBufferInd = sum(activeCellsBuffers,3) > 0;
            redCell = redCellsMask == icell;
            overlapInd = redCell & (bxCellsClickMask > 0 | activeBufferInd);
            twoColorLabel(overlapInd) = false;
        end
    end
    twoColorMask = bwlabel(twoColorLabel);
    figure;imagesc(twoColorMask)
    nCells = length(unique(twoColorMask(:)))-1;
    
    for icell = 1:nRedCells
        if icell == 1
            isLabeledCell = false(nCells,1);
        end    
        cellPix = twoColorMask(round(redCellCenters(icell,2)),round(redCellCenters(icell,1)));
        if cellPix > 0
            isLabeledCell(cellPix) = true;
        end
    end
    save(fullfile(fnout,'twoColorMask&CellLabels'),'twoColorMask',...
        'isLabeledCell','redImage','bxImage','tunImage')
end
