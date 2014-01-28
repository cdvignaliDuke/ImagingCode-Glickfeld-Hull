function [os] = readStackComplete(varargin)
% output is a structure
%
%   warning: lots of data gets saved to disk cache, but it's better than
%   recomputing all the phase shifts etc.
%
%$Id$

defs = { 'ExptName', [], ...
         'SeriesName', [], ...
         'IsLeicaData', true, ...
         'DataPathOverride', [], ...
         'DoRed', true, ...
         'ForceRedImageName', [], ...
         'ForceRedFromSeries', [], ...
         'DoCaching', true, ...
         'DoPhaseAdjust', true, ...
         'DoPhaseFigureSave', true, ...
         'PhaseAdjustNAvgFrames', [], ...
         };

uo = stropt2struct(stropt_defaults(defs, varargin));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% first check disk cache 
if uo.DoCaching
    cacheKey = {lower(uo.ExptName), lower(uo.SeriesName)}; % ignore case
    outdatedKey = struct_trim(uo, ...
                              {'DataPathOverride', 'DoRed', ...
                               'ForceRedImageName', ...
                               'ForceRedFromSeries', 'DoPhaseAdjust'});
    cClass = mfilename; 
    cDat = disk_cache('get', cacheKey, 'allowMissing', cClass, outdatedKey);
    if ~isempty(cDat), [os] = deal(cDat{:}); 
        fprintf(1, '%s: read stack from disk cache %s %s\n', ...
                mfilename, uo.ExptName, uo.SeriesName);
        return; 
    end
end
%%% else not found in cache, continue

%% find data location
[dataPath, isRegisteredTif, origPath] ...
    = subGetDataPath(uo.ExptName, uo.SeriesName, uo.IsLeicaData, ...
                     uo.DataPathOverride);
                      

if uo.IsLeicaData
    
    %% read green tiffs 
    if isRegisteredTif
        % path hardcoded in subGetDataPath
        os.greenStack = readtiff(dataPath);
    else
        % leica raw files on disk
        
        %% autodetect which channel to use
        [chanNo redChanNo] = subLeicaGuessChanNos(origPath, uo.SeriesName);
        
        restrictRE = sprintf('%s_t[0-9]+_ch%02d', uo.SeriesName, chanNo);
        os.greenStack = readtiff(dataPath, [], restrictRE, true);
    end
    [nRows nCols nFrames] = size(os.greenStack);

    %% do red if requested
    if uo.DoRed
        % which series to get from?
        if ~isempty(uo.ForceRedFromSeries)
            redSeriesName = uo.ForceRedFromSeries;
        else
            redSeriesName = uo.SeriesName;
        end
        
        [os.redImg lagsUsed phaseFigH] = subReadRedImage('ExptName', uo.ExptName, ...
            'SeriesName', redSeriesName, ...
            'IsLeicaData', uo.IsLeicaData, ...
            'DataPathOverride', uo.DataPathOverride, ...
            'ForceRedImageName', uo.ForceRedImageName, ...
            'DoPhaseAdjust', uo.DoPhaseAdjust, ...
            'PhaseAdjustNAvgFrames', uo.PhaseAdjustNAvgFrames, ...
            'DoPhaseFigureSave', uo.DoPhaseFigureSave);

    end
    
    %% do phase adjust, after red
    if uo.DoPhaseAdjust && ~isRegisteredTif 
        if ~isempty(lagsUsed)  % means that we did phase correction in red channel
            % adjust using red shifts
            os.greenStack = imBidirShift(os.greenStack, lagsUsed);
        else
            % compute on green
            [os.greenStack,crap,crap,crap,phaseFigH] ...
                = stackFixBidirPhase(os.greenStack, [], ...
                                     uo.DoPhaseFigureSave, ...
                                     uo.PhaseAdjustNAvgFrames);
        end

        if uo.DoPhaseFigureSave
            %% save phase correct figure (either from red or green above)
            set(phaseFigH, 'Visible', 'Off');
            dirs = directories;
            fprintf(1, 'Saving phase shift figure (from stackBidirPhaseCorrect.m)\n');
            set(phaseFigH, 'Visible', 'off');        
            outName = fullfile(dirs.phaseFigOut, ...
                               sprintf('%s-%s_phase_correct_shift', ...
                                       uo.ExptName, uo.SeriesName));
            exportfigPrint(phaseFigH, outName, ...
                           'FileFormat', 'pdf', ...
                           'Size', 10*[1 0.75]);
            close(phaseFigH);
        end
    end
        
    
else  % if isLeicaData
    error('only leica implemented now');
end
    

if uo.DoCaching
    %%% store data in cache and return
    disk_cache('set', {os}, cacheKey,outdatedKey,'allowDups',cClass);
end
    

