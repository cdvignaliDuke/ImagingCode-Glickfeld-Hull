function iterateStacksQuick(funcHandle, listIndices)
%
%  This reads SkipPhaseCorrect DataPathOverride
%
%  and uses dirs.rawdataRoot and dirs.iterQuickOut and
%           dirs.listfile
%   (see directories.m)
%
%$Id: iterateStacksQuick.m 348 2008-10-03 20:00:17Z histed $

if nargin < 2, listIndices = []; end

dirs = directories;

% start saving a diary
outDiaryName = fullfile(dirs.iterQuickOut, ...
                        sprintf('diary-%s-start.txt', ...
                                datestr(now, 'yymmdd-HHMM')));
diary(outDiaryName);


% read listfile
fs = readListfile(dirs.listfile);

% check for this fn's necessary field names
fNames = fieldnames(fs);
subCheckForField(fNames, 'SkipPhaseCorrect');
subCheckForField(fNames, 'DataPathOverride');

% restrict by indices
indicesToDo = 1:fs.nRows;
if ~isempty(listIndices)
    misIndices = listIndices(~ismember(listIndices, indicesToDo));
    if ~isempty(misIndices)
        error('specified some indices not in listfile: %s', mat2str(misIndices));
    end
    indicesToDo = listIndices;
end


%%% iterate through listfile, 
for iR = indicesToDo
    
    [tExptName tSeriesName] = listfileGetExptId(fs, iR);

    if isempty(tExptName) || isempty(tSeriesName)
        error('Expt/series name is empty; missing rows or end of file?');
    end

    % message user
    fprintf(1, '\n-----------------------------------------------\n');
    fprintf(1, '%s: reading stack %d/%d: %s %s\n', ...
            mfilename, iR, fs.nRows, ...
            tExptName, tSeriesName);
    
    %% get data
    [le] = listfileGetAllData(tExptName, tSeriesName);

    %% define output dir   
    outDir = fullfile(dirs.iterQuickOut, [tExptName '_' tSeriesName]);
    if ~exist(outDir, 'dir')
        mkdir(outDir);  % create if missing
        fprintf(1, 'Created output directory %s\n', outDir);
    else
        rmdir(outDir, 's'); % delete it and its contents
        mkdir(outDir); % recreate
        fprintf(1, 'Cleared and recreated output directory %s\n', outDir); 
    end

% $$$     % delete the fields from ls that we don't want to pass
% $$$     args = subExtractArguments(fs, iR, ...
% $$$                                { 'ExptName', ...
% $$$                                  'SeriesName', ...
% $$$                                  'DataPathOverride', ...
% $$$                                  'SeriesEstimNum', ...
% $$$                                  'IsLeicaData', ...
% $$$                                  'SkipPhaseCorrect' });

    % call function with listfile entry
    feval(funcHandle, ...
          le, ...
          outDir);

          
    % done, iterate
end

diary('off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subCheckForField(fNames, desName);
if ~ismember(desName, fNames)
    error('%s: missing listfile field: %s', mfilename, desName);
end


