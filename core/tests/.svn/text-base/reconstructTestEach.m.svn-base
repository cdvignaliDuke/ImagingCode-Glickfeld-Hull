function reconstructTestEach(dadDataRoot, outRoot, tExptName, tSeriesName, ...
                             ijd2tOpts)
%
% Reconstruct a data series, run some computations and averages, write
% output figures and images, return
%
% 080619 histed: v1
% 
% $Id$

%% consts
dataSampleFreq = 3.373e6;    % check with vincent

if nargin < 5, ...
        ijd2tOpts = { 'FilesN', [], ... % 20, ...%1:5:20, ...
                      'Channels', [1 1], ...
                      'AverageN', 1, ...
                      'Limits', [2048 4096], ...
                      'FixDelay', [], ...
                      'Depth', 16, ...
                      'RecalculateDelay', true, ...
                      'ContinuousData', true };
end


%%%%%%%%%%%%%%%%

%% fill in consts
sourceName = fullfileMH(dadDataRoot, tExptName, tSeriesName);
outName = fullfileMH(outRoot, sprintf('%s-%s', tExptName, tSeriesName));
outFigDir = fullfileMH(outRoot, 'figures');
if exist(outFigDir,'dir')~=7, mkdir(outFigDir); end
outFigBasename = fullfileMH(outFigDir, sprintf('%s-%s', tExptName, tSeriesName));
outMatName = [outFigBasename '-statoutput.mat'];

%% do reconstruct
if exist(outName, 'dir')==7, deltree(outName); end% clear output dir
warning off;
[periods stats] = ijdad2tiffseq(sourceName, outName, ...
                                ijd2tOpts{:});

% remove holes in periods
periods = periods(periods ~= 0);
assert(length(periods) == stats.nUsedFiles);

% save periods, stats
save(outMatName, 'stats', 'periods');

load(outMatName, 'stats', 'periods');

%%%%%%%%%%%%%%%%


%% figures
figH = [];

% period progression
figH(end+1) = figure;
periodsUs = (periods / dataSampleFreq) / 1e-6;
%periodsUs = periodsUs(periodsUs~
plot(periodsUs, 'x-');
xlabel('File #');
ylabel('Fast mirror periods (\mus)');
yLim = get(gca, 'YLim');
yLim = [min(yLim(1), 249.9), max(yLim(2), 250.1)];  % set minimum y limits
set(gca, 'YLim', yLim);
title(sprintf('%s-%s: Dad file calculated periods', tExptName, tSeriesName));

% delay progression
figH(end+1) = figure;
plot(stats.delayForFile, 'x-');
xlabel('File #');
ylabel('Delay (pixels)');
%yLim = get(gca, 'YLim');
%yLim = [min(yLim(1), 18), max(yLim(2), 22)];  % set minimum y limits
%set(gca, 'YLim', yLim);
title(sprintf('%s-%s: Dad file delays used in reconstruction', tExptName, tSeriesName));

% save figures to disk
set(figH, 'Visible', 'off');
nFigs = length(figH);
for iF = 1:nFigs
    tOutName = [outFigBasename, sprintf('-fig%03d', iF)];
    exportfigPrint(figH(iF), tOutName, ...
                   'FileFormat', 'png', ...
                   'Size', 8 * [1 0.75]);
    close(figH(iF));
end

%%%%%%%%%%%%%%%%

%%% analysis from tiffs on disk

%% get whole stacks
gStack = readtiff(stats.targetGreen);
[nRows nCols nFrames] = size(gStack);
rStack = readtiff(stats.targetRed);
assert(all(size(rStack) == size(gStack)));

%% whole mean of red
gFov = mean(gStack,3);
gFovU8 = uint8(floor(imScale(gFov, [], [0 255])));
imwrite(gFovU8, [outFigBasename, '-fovGreen.png'], 'png');
rFov = mean(rStack,3);
rFovU8 = uint8(floor(imScale(rFov, [], [0 255])));
imwrite(rFovU8, [outFigBasename, '-fovRed.png'], 'png');

%% pull out first/last few frames in dad
dadGStack = subDadStartEndMeanFrames(gStack, stats.framesPerDadFile);
writetiff(dadGStack, [outFigBasename '-dadAvgFirstLastGreen']);
dadRStack = subDadStartEndMeanFrames(rStack, stats.framesPerDadFile);
writetiff(dadRStack, [outFigBasename '-dadAvgFirstLastRed']);
