function rgbFig = plotRGBImage(imageStack)
%     image should already be normalized
    if length(size(imageStack)) == 2
        nframes = 1;
    else
        [ypix,xpix,nframes] = size(imageStack);
    end
    rgbFig = figure;
    if nframes == 3
        imagesc(imageStack)
    elseif nframes < 3
        emptyRGBframes = zeros(ypix,xpix,3-nframes);
        paddedImageStack = cat(3,imageStack,emptyRGBframes);
        imagesc(paddedImageStack);
    elseif nframes > 3
        error('image stack has too many frames')
    end

end