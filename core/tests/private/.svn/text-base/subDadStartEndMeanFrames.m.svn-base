function outStack = subDadStartEndMeanFrames(stack, framesPerDadFile)
%% make a stack with averages of first 30 and last 30 frames in each dad file
%$Id$

nFrToAvg = 30;

[nRows nCols nFrames] = size(stack);

firstFrNs = cumsum([1 framesPerDadFile]);
assert(firstFrNs(end) == (nFrames+1));
firstFrNs = firstFrNs(1:end-1); % remove last entry, same as nFrames
nDadFiles = length(framesPerDadFile);


% loop over dad files and compute
outFr = repmat(NaN, [nRows, nCols, nDadFiles*2]);
for iF = 1:nDadFiles
    tFirst = firstFrNs(iF);
    tOutN = (iF-1)*2+1;    
    
    % check for errors
    if framesPerDadFile(iF) < (nFrToAvg*2) %any(tFRange >= tLRange)
        fprintf(1, 'Too few frames (%d) in dad file (file %d)\n', ...
                framesPerDadFile(iF), iF);
        outFr(:,:,tOutN:tOutN+1) = 0;
        continue
    end
    
    % do average, beginning frames
    tFRange = tFirst:(tFirst+nFrToAvg-1);
    outFr(:,:,tOutN) = mean(stack(:,:,tFRange),3);
    
    % do average, end frames
    if iF == nDadFiles
        tLast = nFrames;
    else
        tLast = firstFrNs(iF+1)-1;
    end
    tLRange = (tLast-nFrToAvg+1):tLast;
    outFr(:,:,tOutN+1) = mean(stack(:,:,tLRange),3);
    
end

% rescale double to integers
outStack = uint16(floor(imScale(outFr, [], [0 65535])));
