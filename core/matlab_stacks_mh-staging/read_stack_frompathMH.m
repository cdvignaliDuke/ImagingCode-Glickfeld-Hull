function [stack] = read_stack_frompathMH(fnameFmt, frameNums, doCaching)
%read_stack_frompathMH: read a set of tiff files into a 3d matrix
%
%   [stack] = read_stack_frompathMH(formatStr, frameNums, doCaching)
%
%$Id: read_stack_frompathMH.m 80 2008-01-18 00:08:56Z histed $

    
% process arguments
if nargin < 3, doCaching = false; end 

if doCaching
    %%% check mem cache %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cacheKey = {fnameFmt,frameNums};
    outdatedKey = {}; 
    cDat = memory_cache('get', cacheKey, 'allowMissing', outdatedKey);
    if ~isempty(cDat), [stack] = deal(cDat{:}); return; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% test if dir exists
[fPath fName fExt] = fileparts(fnameFmt);
if ~exist(fPath)
    error('Path does not exist: %s', fPath);
end

% test if last image is there
firstImg = sprintf(fnameFmt, frameNums(1));
if ~exist(firstImg, 'file')
    error('First file (%d) not found, is fmt str correct? (name: %s)', ...
          1, firstImg);
end
lastImg = sprintf(fnameFmt, frameNums(end));
if ~exist(firstImg, 'file')
    error('First file (%d) not found, is numbering correct? (name: %s)', ...
          frameNums(end), lastImg);
end



nImages = length(frameNums);

startT = double(tic);
lastT = startT;
for iI=1:nImages
    tInd = frameNums(iI);
    tImage = imread(sprintf(fnameFmt, tInd));

    if iI == 1
        % set up array
        assert(isa(tImage, 'uint8')||isa(tImage, 'uint16'));
        stack = repmat(tImage*0, [1 1 nImages]);
        fprintf(1, '%s: reading %d images\n', mfilename, nImages);
    end
    
    % copy image in
    stack(:,:,iI) = tImage;
    
    % disp
    elapsedS = (double(tic)-lastT)/1e9;
    totalElS = (double(tic)-startT)/1e9;
    if elapsedS > 5 || iI == nImages
        fprintf(1,'%6d/%d  %5.1fs elapsed\n', iI, nImages, totalElS);
        lastT = double(tic);  % restart counter
    end
end

% display only at end
if iI == nImages
    nB = whos('stack');
    fprintf(1, 'done.  %8gMB total, %8gMB/sec, %8g images/sec\n', ...
            chop(nB.bytes ./ 1e6, 3), ...
            chop(nB.bytes ./ 1e6 / totalElS, 3), ...
            chop(nImages / totalElS, 2));
end


if doCaching
    %%% store data in cache %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    memory_cache('set', {stack}, cacheKey,outdatedKey,'allowDups');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end





