function [gZ rZ gStack rStack] ...
    = frVisZStack(exptName, seriesName, doOverwrite, splitPlanes);
%[gZ rZ gStack rStack] ...
%   = frVisZStack(exptName, seriesName, doOverwrite, splitPlanes);
%
%   splitPlanes: bool {false}|true:  if true, make separate stacks for
%   each color channel; advantage is that separate channels can each be
%   32-bit while RGB is limited to 8bits/channel
%
%$Id$

if nargin < 1, exptName = 'calib081230'; end
if nargin < 2, seriesName = 'zstack3A'; end
if nargin < 3, doOverwrite = false; end
if nargin < 4, splitPlanes = false; end

%% consts
%exptName = 'calib081230';
%seriesName = 'zstack3A';

%% compute relevant directories
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);

%% data is copied to j:\data\filename
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false)
  

%% get info
sl = frParseScanLog(dadDir);

%% read files
imDir = fullfileforw(dirs.images, exptName, seriesName);
gStack = readtiff(fullfileforw(imDir, [seriesName '_green']));
rStack = readtiff(fullfileforw(imDir, [seriesName '_red']));


%% average
[nRows nCols nFrames] = size(gStack);
expectedNFrames = sl.cycle_NSteps * sl.cycle_NFramesPrStep;
if nFrames < expectedNFrames
    warning('Not enough frames found: if acq not stopped early, likely bug');
elseif nFrames > expectedNFrames
    error('Too many frames?  bug somewhere');
end
gZ = stackGroupProject(gStack, sl.cycle_NFramesPrStep);
rZ = stackGroupProject(rStack, sl.cycle_NFramesPrStep);

if splitPlanes
    redI = ijarray2plus(rZ, [], [sl.seriesName ' red']);
    greenI = ijarray2plus(gZ, [], [sl.seriesName ' green']);    
    ijStart;
    redI.show;
    greenI.show;
else
    % make RGB image
    rgbZ = cat(4, rZ, gZ, gZ*0);
    rgbZ8 = imScale(rgbZ, [], [0 255], 'uint8');

    %% open in imagej
    rI = ijarray2plus(rgbZ8, [], sl.seriesName);
    ijStart;
    rI.show;
end

%% open stack
%rsIJ = ijarray2plus(rStack);
%rsIJ.show;
