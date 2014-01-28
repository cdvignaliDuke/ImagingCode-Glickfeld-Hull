function [decStack outDir outName] ...
        = frReadWithDecimation(exptName, seriesName, ...
                               channelName, ...
                               doWrite, ...
                               timeDecFrames)
%frReadWithDecimation(fastrig): read stack, avg frames, write dec. stack
%
%   [decStack outDir outName] ...
%       = frReadWithDecimation(imageDir, doWrite, timeDecFrames)
%       
%       directory is: dirs.images/exptName/seriesName/seriesName_channelName
%       imageDir: 'exptName/seriesName/
%       doWrite: boolean {true}: output decimated stack
%       if output argument is present, return decimated stack
%
%   output stacks are written to 
%       dirs.images/exptName/seriesName/binned_channelName
%
%   Mainly useful for online analysis
%
%$Id$

dirs = frGetDirs;

assert(any(ismember({'red', 'green'}, channelName)), ...
       'channelName must be red or green');
if nargin < 4, doWrite = true; end
if nargin < 5, timeDecFrames = 10; end

imDir = fullfileforw(dirs.images, exptName, seriesName, ...
                     [seriesName '_' channelName]);
outDir = fullfileforw(dirs.images, exptName, seriesName, ...
                      ['binned' '_' channelName]);
if ~exist(outDir)
    mkdir(outDir);
end

% read one tiff file at a time
dirS = dir(fullfileforw(imDir, '*.tif*'));
fNames = {dirS.name};
nFiles = length(fNames);

% read first file always
tStack = readtiff(imDir, 1, [], false);
[nRows nCols nFrames] = size(tStack);
for iF = 1:nFiles
    if iF ~= nFiles  % can't do for last
        nextStack = readtiff(imDir, iF+1, [], false);
    end
    
    nDecFrames = ceil(nFrames./timeDecFrames);
    nTotalFr = nDecFrames * timeDecFrames;
    nNeededFromNext = nTotalFr - nFrames;
    
    if iF == nFiles
        % last, can't take extra frames
        tDecStack = stackGroupProject(tStack, timeDecFrames);
    elseif nNeededFromNext == 0
        tDecStack = stackGroupProject(tStack, timeDecFrames);
        tStack = nextStack;
    elseif nNeededFromNext > 0
        tDecStack = stackGroupProject(cat(3,tStack,...
                                          nextStack(:,:,1:nNeededFromNext)), ...
                                      timeDecFrames);
        tStack = nextStack(:,:,nNeededFromNext+1:end);
    else
        error('bug');
    end
    [nRows nCols nFrames] = size(tStack);    
    
    % save to cell array
    if nargout > 0
        decC{iF} = tDecStack;
    end
end

decStack = cat(3, decC{:});

% write output
if doWrite
    outName = ['binned_' seriesName '.tif'];
    writetiff(decStack, fullfileforw(outDir,outName));
end
    


if nargout == 0
    clear decStack outDir outName  % return nothing
end

    



