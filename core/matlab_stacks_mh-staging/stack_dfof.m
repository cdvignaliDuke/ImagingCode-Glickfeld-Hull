function outStack = stack_dfof(stack, F0, nGaussPix, pctOffset)
%
%$Id: stack_dfof.m 58 2007-11-15 21:38:14Z histed $

if nargin < 4, pctOffset = 15; end

[nRows,nCols,nInFrames] = size(stack);

% smooth with gauss
F0 = double(F0);
smFr = smooth2(F0, 'gauss', [2 2], nGaussPix);

smOffset = range(smFr(:)) * pctOffset/100;
fprintf(1, 'Smoothed frame: max %5.2f, min %5.2f, offset by %5.2f\n', ...
        max(smFr(:)), min(smFr(:)), smOffset);
smFr = smFr + smOffset;


%figure;
%imshow(smFr./max(smFr(:)));

outStack = double(stack);
for iF = 1:nInFrames
    % do this inplace
    tFr = (double(stack(:,:,iF)) - F0) ./ smFr;
    outStack(:,:,iF) = tFr;
end

% rescale to full frame
maxV = max(outStack(:));
minV = min(outStack(:));

outStack = uint8(floor( ((outStack - minV) ./ (maxV-minV)) * 255));    



