function sl = frParseScanLog(dirName)
%FRPARSESCANLOG
% SL = FRPARSESCANLOG(DIRNAME)
%
% dirName is of the form 'exptName/seriesName' and dad files usually at
%  j:/data/exptName/seriesName/seriesName.dad 
%  etc.
%
%
%
%$Id$

dirs = frGetDirs;

if exist(dirName, 'dir')
    sl.logPath = fullfileforw(dirName, '00fastscan-data-log.txt');
    % old format:
    %sl.logPath = fullfileforw(dirName, 'fastscan-data-log.txt');    
elseif exist(dirName, 'file')
    % full path passed in, do nothing
    [tPath tName tExt] = fileparts(dirName);
    assert(strcmp([tName tExt], '00fastscan-data-log.txt'), ...
        'Passed in an existing file but not a scan log file?');
    sl.logPath = dirName;
else
    % add data dir to it
    sl.logPath = fullfileforw(dirs.data, dirName, '00fastscan-data-log.txt');
end

if ~exist(sl.logPath, 'file');
    % if can't find the log, return empty
    sl = []; 
    return
end

% parse expt name
[tempPath tFN] = fileparts(sl.logPath);
[tempPath sl.seriesName] = fileparts(tempPath);
[tempPath sl.animalName] = fileparts(tempPath);

% read the whole file into a cell array
fid = fopen(sl.logPath, 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
C = C{1};

%% pull out fields 
% first structures
tokC = regexp(C, 'fs\.(.*): (.*)$', 'tokens');
noMatchIx = cellfun(@isempty, tokC);
toks1 = cat(1, tokC{:});
toks = cat(1, toks1{:});
fNames = toks(:,1);
fValStrs = toks(:,2);
nFields = length(fNames);

% then bare names
tokC = regexp(C, '^([A-Za-z0-9_]*): (.*)$', 'tokens'); 
       %will not match fs. lines because of dot
noMatchIx = cellfun(@isempty, tokC);
toks1 = cat(1, tokC{:});
toks = cat(1, toks1{:});
if ~isempty(toks)
    fNames = cat(1, fNames, toks(:,1));
    fValStrs = cat(1, fValStrs, toks(:,2));
    nFields = length(fNames);
end

%assert(length(unique(fNames)) == length(fNames), ...
%       'bug: duplicate field names?');

[cr1 keepIx] = unique(fNames);
fNames = fNames(keepIx);
fValStrs = fValStrs(keepIx);
nFields = length(fNames);

%fprintf(1, 'debug! fix dup field names\n');

% make into a struct
for iF = 1:nFields
    tFName = strrep(fNames{iF}, '.', '_');
    sl.(tFName) = str2num(fValStrs{iF});
end

%% extract time stamps

% start time
tag = 'Current time: ';
ind = strmatch(tag,C);
if length(ind)
    l = C{ind};
    sl.start_date = l(findstr(l,tag)+length(tag):end);
end

%
tag = 'Stopped at: ';
ind = strmatch(tag,C);
if length(ind)
    l = C{ind};
    sl.stop_date = l(findstr(l,tag)+length(tag):end);
end

% return

    
% old log format:    
%sl.DAQ_dad_size_in_frames = 150;
