close all
awFSAVdatasets
for iexp = 1:size(expt,2)
%     iexp = 5;

load(fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'mask&TCDir.mat'));
   
m = mask_cell(:);
nCells = length(unique(m))-1;

[nPixPerCellIncluding0 pixBinIx] = histc(m,unique(m));

nPixPerCell = nPixPerCellIncluding0(2:end); %get rid of 0 bin, which will be really high
binEdges = [4:4:400];
[pixBin b] = histc(nPixPerCell,binEdges);
%%
dendritesVsCells = figure;
%%
figure(dendritesVsCells);
subplot(2,2,1)
hist(nPixPerCell,100);
title([expt(iexp).mouse '-' expt(iexp).date]);


%% pixel threshold
pixCutoff = 50;
dendriteIx = find(nPixPerCell < pixCutoff);
cellIx = setdiff(1:nCells, dendriteIx);
figure(dendritesVsCells);
subplot(2,2,1)
hold on
vline(pixCutoff,'k')
%% show mask of cells vs dendrites
dendrites = find(ismember(m,dendriteIx));
dendritePix = zeros(size(m));
dendritePix(dendrites) = m(dendrites);
dendriteMask = reshape(dendritePix,(size(mask_cell)));
d1Mask = dendriteMask;
d1Mask(d1Mask >= 1) = 0.5;

cells = find(ismember(m,cellIx));
cellPix = zeros(size(m));
cellPix(cells) = m(cells);
cellMask = reshape(cellPix,(size(mask_cell)));
c1Mask = cellMask;
c1Mask(c1Mask >= 1) = 0.5;

%%
figure(dendritesVsCells);
subplot(2,2,2)
imagesc(dendriteMask)
title('dendrites')

figure(dendritesVsCells);
subplot(2,2,3)
imagesc(cellMask)
title('cells')

%% image of cells
try
img = readtiff(fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'DirectionTuning_V1.tif'));
imgtxt = 'max dF/F';
dF = bsxfun(@minus,img,(mean(img,3)));
dFoverF = bsxfun(@rdivide,dF,(mean(img,3)));
maxDFoverF = max(dFoverF,[],3);
croppedMaxDFoverF = maxDFoverF;
croppedMaxDFoverF(:,[1:50 750:end]) = 0;
croppedMaxDFoverF(1:50,:) = 0;
catch
    load(fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'regImg.mat'));
    img = data_avg;
    imgtxt = 'registration img';
end
shadeimg = imShade(croppedMaxDFoverF,d1Mask,c1Mask);
%%
figure(dendritesVsCells);
subplot(2,2,4)
imagesc(shadeimg)
title(imgtxt)

%% save
figure(dendritesVsCells)
set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
print(fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,'cellSizeHistImgDendriteCutoff'),'-dpdf')
%%
  end