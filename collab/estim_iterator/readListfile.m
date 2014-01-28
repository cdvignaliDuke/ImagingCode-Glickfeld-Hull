function listStruct = readListfile(fileName)
%readListfile (calcium): read a csv listfile from disk
%   listStruct = readListfile(fileName)
%
%   listStruct is a single struct of length 1
%      First line in the csv file is a header line of strings identifying
%      each column.  
%      Fields of listStruct have the same *name* as the string in the
%      corresponding column of the header line
%      Field values are vectors; numeric if possible and cell arrays otherwise
%
%   histed 04/25/08: created
%
%$Id: readListfile.m 394 2008-11-21 07:02:34Z histed $

if nargin < 1
    dirs = directories;
    fileName = dirs.listfile;
end

% this is here to allow future smarter file reading (i.e. of text files)
% and to provide documentation


ls = csvread_textornum(fileName, true); %caching handled here
                                                %+let mixed fields stay cell
ls.nRows = length(ls.ExptName);

% change this to boolean right up (NaNs get converted to false)
%ls.IsRegisteredTif = ls.IsRegisteredTif == true;

%% extra estim stuff
% convert estimNum field to seriesName
estimName = eval([ '{' sprintf('''estim%d'',', ls.SeriesEstimNum) '}'])';
% add here!  check to make sure only one of SeriesEstimNum or SeriesName is
% set

% copy in seriesNames
mergedNames = ls.SeriesName;
nameBlankIx = cellfun(@isempty, ls.SeriesName);
mergedNames(nameBlankIx) = estimName(nameBlankIx);
ls.specifiedSeriesName = ls.SeriesName;
ls.SeriesName = mergedNames;


% kludge this- if current is a string mat, eval it to a double mat
for iR = 1:ls.nRows
    tC = ls.Current{iR};
    if ischar(tC) && ~isempty(tC)
        ls.Current{iR} = eval(tC);
    end
end

listStruct = ls;
