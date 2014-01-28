
function [filelist, nFiles] = GetMatFileNames(pathname, fileind, textInName, textNotInName)

%%%%%%%% extracted from readtiff.m

if nargin < 2, fileind = []; end
if nargin < 3, textInName = ''; end
if nargin < 4, textNotInName = ''; end

if exist(pathname)  == 7 % directory
    %% list directory
    pathstr = pathname;
    list = dir(pathname); 
    list(1:2)=[];  % remove '.' and '..'
    filelist = {list.name};
    if isempty(filelist)
        error('No files found at path %s', pathname);
    end

    %% first restrict by indices passed in
    if ~isempty(fileind)
        filelist = filelist(fileind);
    end
    
    %% remove non-MPDs
    rMatchNs = regexpi(filelist, '\.(mat)$');
    filelist = filelist(~cellfun(@isempty, rMatchNs));
    
    %% restrict list by text str if necessary
    if ~isempty(textInName)
        rMatchNs = regexpi(filelist, textInName);
        filelist = filelist(~cellfun(@isempty, rMatchNs));
    end
    if ~isempty(textNotInName)
        rMatchNs = regexpi(filelist, textNotInName);
        filelist = filelist(cellfun(@isempty, rMatchNs));
    end
        
    %---  must have full list by the time we hit here, sorting below

    nFiles = length(filelist);

    %% check number resulting
    if nFiles<1
        error('No files remaining after restrictions applied');
    else
        fprintf(1,'Found %i files \n', nFiles);
    end

else  % could not find directory/file
    error('pathname does not exist: %s', pathname);
end