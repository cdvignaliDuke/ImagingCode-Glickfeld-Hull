function outS = stackDecimate(inS, nDec, outClass, doRescale)
%STACKDECIMATE (calcium): running avg groups of frames (like IJ Grouped ZProj)
%
%   outS = stackDecimate(inS, nDec, outClass, doRescale)
%   
%   outClass defaults to 'double'
%   doRescale (bool): true means map range after average to full range of
%   output class.  Default false
%
%$Id$

if nargin < 3, outClass = []; end

% reshape and average
[nRows,nCols,nFrames,nPlanes] = size(inS);
nOutFrames = floor(nFrames/nDec);
outS = mean(reshape(permute(inS(:,:,1:nOutFrames*nDec,:), ...
                            [1 2 4 3]), ...
                    [nRows nCols nPlanes nOutFrames nDec]), ...
            5);
outS = squeeze(permute(outS, [1 2 4 3]));

% convert to output type

% use cast and turn warning off here, instead of using floor, because it
% is likely faster
if ~isempty(outClass)   %(else leave as double)
    ws = warning('off', 'MATLAB:intConvertNonIntVal');
    outS = cast(outS,outClass);
    warning(ws);
end



