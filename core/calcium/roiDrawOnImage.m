function imgH = roiDrawOnImage(bwMask, maskRGBColor)
%roiDrawOnImage (calcium): compose image with RGB roi mask
%
%$Id: roiDrawOnImage.m 104 2008-03-18 02:57:30Z histed $

%%% check opts
if nargin < 2 || isempty(maskRGBColor), maskRGBColor = [0 0 1]; end
assert(length(maskRGBColor) == 3);

[nRows nCols] = size(bwMask);

axH = newplot;
hold on;

maskImg = repmat(bwMask, [1 1 3]);

maskImg = maskImg .* repmat(permute(maskRGBColor(:),[2 3 1]), ...
                            [nRows nCols 1]);

imgH=image(maskImg);
set(imgH, 'Tag', 'draw_roi_on_image-ROI', ...
        'AlphaData', bwMask*0.5);

