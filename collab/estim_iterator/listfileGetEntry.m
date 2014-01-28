function [lOut] = listfileGetEntry(exptName, seriesName)
%
%  read listfile data and return
%
%$Id: listfileGetEntry.m 394 2008-11-21 07:02:34Z histed $

lf = readListfile;

tN = find(strcmpi(lf.ExptName, exptName) & strcmpi(lf.SeriesName, seriesName));
if isempty(tN)
    error('Listfile entry not found: %s %s', exptName, seriesName);
elseif length(tN) > 1
    error('More than one listfile match for: %s %s', exptName, ...
          seriesName);
end
lOut.listfileIndexN = tN;




fNames = fieldnames(lf);
nFields = length(fNames);
for iF = 1:nFields
    tFN = fNames{iF};
    tFV = lf.(tFN);
    if length(tFV) == 1
        assert(strcmp(tFN, 'nRows'));
        continue;
    end
    
    if iscell(tFV)
        tVal = tFV{tN};
    else
        tVal = tFV(tN);
    end
    lOut.(tFN) = tVal;
end


%% handle fastscan sorting here

%% determine what type of data
switch lower(lOut.LeicaOrFastscan)
    case 'fastscan'
        lOut.isFastscan = true;
    case {'', 'leica'}
        lOut.isFastscan = false;
    otherwise
        error('Invalid value: %s ', lOut.LeicaOrFastscan);
end

%% fastscan details
if lOut.isFastscan
    assert(isempty(lOut.FrameTimeMs) || isnan(lOut.FrameTimeMs));
    lOut.FrameTimeMs = 31.25;  % 4000/128
end

%% one-off force for cat data
%assert(strcmpi(exptName, 'cat061012') && strcmpi(seriesName, 'estim19'));
%lOut.DataPathOverride = 'cat061012/cat061006_estim19';

%keyboard
    
