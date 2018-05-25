clear all
close all
ds = 'FSAV_V1_GAD';
rc = behavConstsAV;
eval(ds)
for iexp = 2:size(expt,2)
    
    
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    redFolder = expt(iexp).redChannelRun;
    redChannelOn = expt(iexp).greenredsimultaneous;
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,'data processing');


    load(fullfile(fn,'redChannelImage'),'redImage')
    
    [ypix,xpix] = size(redImage);
    redCellsClickLabel = imCellEditInteractive(redImage);
%     redImageLeftovers = redImage;
%     redBuffer = sum(imCellBuffer(bwlabel(redCellsClickLabel),2),3) + redCellsClickLabel;
%     redImageLeftovers(redBuffer > 0) = 0;
%     cellSelectQuads = {[1, floor(xpix/2), 1, floor(ypix/2)];...
%         [floor(xpix/2)+1, xpix, 1, floor(ypix/2)];
%         [1, floor(xpix/2), floor(ypix/2)+1, ypix];
%         [floor(xpix/2)+1, xpix, floor(ypix/2)+1, ypix]};
%     polySelectLabel = nan(ypix,xpix,4);
%     for iquad = 1:4
%         polySelectLabel(:,:,iquad) = imCellPolyEditInteractive(redImageLeftovers,[],cellSelectQuads{iquad});
%     end
%     redLeftoversLabel = sum(polySelectLabel,3);
%     redCellsLabel = logical(redLeftoversLabel+redCellsClickLabel);
    redCellsMask = bwlabel(redCellsClickLabel);
    
    save(fullfile(fn,'redChannelMask'),'redCellsMask')
end