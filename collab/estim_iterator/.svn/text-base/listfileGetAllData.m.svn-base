function [os] = listfileGetAllData(exptName, seriesName, doCaching, ...
    redoCache, forcePhaseSkip)
%
%  get listfile entry and stack data
%
%$Id: listfileGetStack.m 199 2008-04-29 21:11:01Z histed $

if nargin < 3, doCaching = true; end
if nargin < 4, redoCache = false; end
if nargin < 5, forcePhaseSkip = false; end

exptName = lower(exptName);  % so don't need to worry about case
seriesName = lower(seriesName);

%% get all metadata
le = listfileGetEntry(exptName, seriesName);

%% first check caches
if ~le.isFastscan
    % disk cache: for Leica only
    cacheKey = {exptName, seriesName}; outdatedKey = {}; cClass = mfilename;

    if doCaching && ~redoCache
        cDat = disk_cache('get', cacheKey, 'allowMissing', cClass, outdatedKey);
        if ~isempty(cDat), 
            [os] = deal(cDat{:});
            fprintf(1, '%s: read stack from disk cache %s %s\n', ...
                mfilename, exptName, seriesName);
            return;
        end
    end
else
    % memory cache for fastscan
    % only allow one file in cache now
    cacheKey = {mfilename}; outdatedKey = {exptName,seriesName};
    cDat = memory_cache('get', cacheKey, 'allowMissing', outdatedKey);
    if ~isempty(cDat)
        [os] = deal(cDat{:});

        fprintf(1, '%s: got fastscan stack from memory cache %s %s\n', ...
            mfilename, exptName, seriesName);
        return;
    end
end

    
%%% else not found in cache, continue

           
%% read data from disk

if le.isFastscan
    %% read fastscan here
    % things to add:
    %    break to a separate function
    %    subset red stack
    %    detect both kinds of reconstruction directories
    %    don't check cache for fastscan?

    % hack it here; should be a separate function
    dirs = directories;
    
    %rootDir = fullfile(dirs.fastscanDataRoot, lower(exptName), lower(seriesName));
    % default to registered
    dirs = frGetDirs;
    if exist(fullfile(dirs.registered, lower(exptName)), 'dir')
        rootDir = fullfile(dirs.registered, lower(exptName), lower(seriesName));
    elseif exist(fullfile(dirs.images, lower(exptName)), 'dir')
        rootDir = fullfile(dirs.images, lower(exptName), lower(seriesName));
    else
        error('reading fast scan images: none found');
    end
    ds.greenStack = readtiff(fullfile(rootDir, 'green'), []);
    
    % read red stack
    %    subsampTotalN =
    rStack = readtiff(fullfile(rootDir, 'red'), []);
    ds.redImg = mean(rStack,3);

else
    %% read leica
    nAvgFrames = le.PhaseNAvgFrames;
    if isnan(nAvgFrames), nAvgFrames = 1; end
    assert(nAvgFrames>=1 && nAvgFrames < 10000, 'check listfile field');
    doPhase = ~(forcePhaseSkip==1) && ~(le.SkipPhaseCorrect==1)
    [ds] ...
        = readStackComplete('ExptName', exptName, ...
                            'SeriesName', seriesName, ...
                            'IsLeicaData', le.IsLeicaData, ...
                            'DataPathOverride', le.DataPathOverride, ...
                            'DoRed', true, ...
                            'DoPhaseAdjust', doPhase, ...
                            'PhaseAdjustNAvgFrames', nAvgFrames, ...
                            'DoPhaseFigureSave', true, ...
                            'ForceRedImageName', le.RedImageOverride, ...
                            'ForceRedFromSeries', le.RedFromSeries, ...
                            'DoCaching', false);  % done in this function
end

%% see if green stack needs truncating
[nRows nCols nFrames] = size(ds.greenStack);
if isnan(le.DoTruncate) ...   % if missing, try it anyway
        || le.DoTruncate == true  % and if specified, definitely do it
    maxStims = floor((nFrames-(le.StimEvery-1)) ./ le.StimEvery);
    maxReps = floor(maxStims ./ le.NStimsInSer);
    assert(maxReps >= 1, ...
           'Too few frames for even 1 rep');
    nStimsToKeep = maxReps*le.NStimsInSer;
    nFramesToKeep = le.StimEvery*nStimsToKeep+(le.StimEvery-1);
    assert(nFramesToKeep <= nFrames, 'bug');
    if nFramesToKeep < nFrames;
        assert(nFramesToKeep > nFrames/2, ...
               'Truncated by a factor of 2?');
        fprintf(1, '\n    *** Too few frames for stim: Truncated green stack to %d frs from %d\n', ...
                nFramesToKeep, nFrames);
        ds.greenStack = ds.greenStack(:,:,1:nFramesToKeep);
    end
else
    %fprintf(1, '\n    ***   Skipping truncate\n');
end


os = struct_union(le, ds);

%% compute an fov
os.fov = mean(os.greenStack,3);
os.fovL = stack_localcontrastadj(os.fov, 31, 5, [], 0.999);
if ~isempty(os.redImg)
    os.redL = stack_localcontrastadj(os.redImg, 31, 5, [], 0.9999);
    os.fovRGB = cat(3, os.redL, os.fovL, 0*os.redL);
    os.fovN_RGB = imScale(cat(3, os.redImg, os.fov, 0*os.fov));
else
    % no red image available
    os.redL = [];
    os.fovRGB = os.fov;
    os.fovN_RGB = os.fovL;
end


%% store data in cache
if doCaching
    
    if le.isFastscan
        %% store in memory cache
        memory_cache('set', {os}, cacheKey,outdatedKey,'allowDups');
    else
        %% store data in disk cache 
        disk_cache('set', {os}, cacheKey,outdatedKey,'allowDups',cClass);
    end

end

    
