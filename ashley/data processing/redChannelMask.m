ds = '_LMgad';
eval(['awFSAVdatasets' ds]);
rc = behavConstsAV;
iexp = 3;
%%
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
redFolder = expt(iexp).redChannelRun;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);

%% load red data
fName = [redFolder '_000_000'];
redData = loadsbx_choosepmt(2,mouse,expDate,redFolder,fName);
if size(redData,3) > 1000
    redData = redData(:,:,1:1000);
end
greenData = loadsbx_choosepmt(1,mouse,expDate,redFolder,fName);
if size(greenData,3) > 1000
    greenData = greenData(:,:,1:1000);
end
%% register
if expt(iexp).greenredsimultaneous
    load(fullfile(fn,'reg2RedOuts&Img.mat'))
    redDataDownSamp = stackGroupProject(redData,10);
    [~,redReg] = stackRegister(redDataDownSamp,image4Reg);
else
    load(fullfile(fn,'regOuts&Img.mat'))
    redDataDownSamp = stackGroupProject(redData,10);
    greenDataDownSamp = stackGroupProject(greenData,10);
    greenRegOuts = stackRegister(greenDataDownSamp,data_corr_img);
    [~,redReg] = stackRegister_MA(redDataDownSamp,[],[],greenRegOuts);
end
redImage = mean(redReg,3);
figure;colormap hot; imagesc(redImage)

%% image borders?
if expt(iexp).areaBorders
    borders = readtiff(fullfile(fn,'FOVborders.tif'));
    redImageNorm = (redImage - min(redImage(:)))./max(redImage(:));
    redImage4CellSelection = redImageNorm + borders;
    figure;colormap hot; imagesc(redImage4CellSelection)
else
    redImage4CellSelection = redImage;
end
%% save image
save(fullfile(fn,'redChannelImage'),'redImage')

%% select cells
bwout = imCellEditInteractive(redImage4CellSelection);
mask_redCell = bwlabel(bwout);
figure;imagesc(mask_redCell);
title(sprintf('%s cells',num2str(length(unique(mask_redCell(:)))-1)))

%% save mask
save(fullfile(fn,['redMask' ds]),'mask_redCell');