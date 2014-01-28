function pixMat = get_tcoursesMH(stack, roiMask)
%GET_TCOURSESMH: return the pixel intensities inside a ROI
%
%   pixMat = get_tcourses(stack, roiMask)
%
%   stack is an image/stack, of size (nRows,nCols,nFrames)
%       nFrames can be 1 or 0 (i.e. stack can be 2-d)
%   roiMask is a logical matrix, size (nRows,nCols), as returned by
%       ROIPOLY.  nRoiPixels is sum(roiMask(:) == true);
%   
%   pixMat is a matrix of size (nRoiPixels, nFrames) with the intensity
%       values of stack where corresponding pixels in roiMask are true.
%       Each column corresponds to one frame.
%       It is the same numeric type (e.g. uint8, double) as stack.
%
%
%$Id: get_eachpix_tcoursesMH.m 59 2007-12-11 01:00:07Z histed $

if nargin < 2, error('Must specify both stack and roiMask'); end

[nRows,nCols,nFrames] = size(stack);
if ~(all(size(roiMask) == [nRows,nCols]))
    error('roiMask must be same size as stack(:,:,1)');
end
nRoiPixels = sum(roiMask(:) == 1);

% old way
% $$$ pixMat = zeros([nRoiPixels,nFrames], class(stack));  % match type
% $$$ for iF = 1:nFrames
% $$$     tFrame = stack(:,:,iF);
% $$$     pixMat(:,iF) = reshape(tFrame(roiMask == true), [nRoiPixels, 1]);
% $$$ end

% faster way
stack = reshape(stack, [nRows*nCols, nFrames]);  % convert images to cols
roiMask = reshape(roiMask, [nRows*nCols,1]);
pixMat = stack(roiMask,:);

