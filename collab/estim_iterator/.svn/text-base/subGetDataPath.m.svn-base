function [dataPath, isRegisteredTif, origPath] = subGetDataPath(exptName, ...
    seriesName, isLeicaData, dataPathOverride);
%
%$Id$


dirs = directories;

%% what is the correct path to the data?
if ~exist(dirs.rawdataRoot)
    error('rawdata root dir does not exist (%s)', dirs.rawdataRoot);
end

isRegisteredTif = false;
if (~isempty(dataPathOverride) ...
    && ~any(isnan(dataPathOverride)))
    % override path set, use it
    origPath = fullfileMH(dirs.rawdataRoot, dataPathOverride);
    dataPath = origPath;
else
    if isLeicaData
        origPath = fullfileMH(dirs.rawdataRoot, exptName, ...
                              [exptName '_' seriesName]);
        
        % check for registered
        testPath = fullfileMH(dirs.rawdataRoot, 'registered', ...
                              [exptName '_' seriesName '.tif']);
        if exist(testPath, 'file')
            dataPath = testPath;
            fprintf(1, 'Using registered stack for %s %s\n', ...
                exptName, seriesName);
            isRegisteredTif = true;
        else
            % use normal leica path
            dataPath = origPath;
        end

    else
        error('only leica / registered tif stacks from leica implemented now');
    end
end

%% check for errors
foundFile = exist(dataPath, 'file');

if ~exist(dataPath, 'dir') && ~foundFile
    error('Data path does not exist: %s', dataPath);
end
if isLeicaData && ~isRegisteredTif && ~foundFile
    %if ~exist(fullfile(dataPath, [exptName '_' seriesName '.lei']),
    %'file');
    if isempty(ls(fullfile(dataPath, '*.lei')))
        error('Target dir does not look like Leica data: %s', dataPath);
    end
end


