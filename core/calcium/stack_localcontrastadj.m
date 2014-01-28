function outStack = stack_localcontrastadj(stack, nGaussPix, pctOffset, ...
    fovImg, topQuantToDiscard)
%
%$Id: stack_localcontrastadj.m 93 2008-03-15 20:30:19Z histed $

if nargin < 3, pctOffset = 8; end
if nargin < 4, fovImg = []; end
if nargin < 5, topQuantToDiscard = 0.9999; end

[nRows,nCols,nInFrames,nColors] = size(stack);

assert(isempty(fovImg) ...
       || (ndims(fovImg) == 2 ...
           && all(size(fovImg) == [nRows nCols])), ...
       'input fov image must be same size as frames of stack');



%% fov
if nInFrames == 1
    fov = double(stack);
elseif isempty(fovImg) 
    fprintf(1, '%s: averaging to get FOV image... ', mfilename);
    fov = mean(stack, 3);
    fprintf(1, 'done.\n');
else
    fov = double(fovImg);
end


% smooth with gauss
smFr = smooth2(fov, 'gauss', [2 2], nGaussPix);
smOffset = range(smFr(:)) * pctOffset/100;
fprintf(1, 'Smoothed frame: max %5.2f, min %5.2f, offset by %5.2f\n', ...
        max(smFr(:)), min(smFr(:)), smOffset);
smFr = smFr + smOffset;

%figure;
%imshow(smFr./max(smFr(:)));

stack = single(stack);
for iF = 1:nInFrames
    % do this inplace
    tFr = stack(:,:,iF) ./ smFr;
    stack(:,:,iF) = tFr;
end

%% rescale to full frame

% first find max and min vals of stack
nToSkip = floor(nInFrames / 100);
if nInFrames <= 100
    subNs = 1:nInFrames;
else
    subNs = 1:nToSkip:nInFrames;  % downsample to 100 fr
end
minV = min(reshape(stack(:,:,subNs),[],1));
if isempty(topQuantToDiscard)
    maxV = max(reshape(stack(:,:,subNs),[],1));
else
    maxV = quantile(reshape(stack(:,:,subNs), [], 1), topQuantToDiscard);
end

% do frame by frame to save mem, rely on JIT to make fast
for iF = 1:nInFrames
    stack(:,:,iF,:) = floor( ((stack(:,:,iF,:) - minV) ./ (maxV-minV)) * 255 );    
end
% no warn on overflow
wS = warning('off', 'MATLAB:intConvertOverflow');
stack = uint8(stack);
warning(wS);

% note; nearly 2x as long to do uint8(min(stack,255)) than just uint8(stack)
% for 512x512x1 stack:          = 0.012s                        = 0.007s

%stack = uint8(floor( ((stack - minV) ./ (maxV-minV)) * 255));    
%stack = stack_rescale(stack);

outStack = stack;


