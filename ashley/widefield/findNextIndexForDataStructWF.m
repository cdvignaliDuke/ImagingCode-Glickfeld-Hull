function [xd indexRowN] = findNextIndexForDataStructWF(rc, dS); %(varargin)
%  read subj xls, find a record
% 
% histed 121003

%% arg processing
%userDefs = {};
%uo = stropt2struct(stropt_defaults(userDefs, varargin));

if ~exist(rc.indexFilename, 'file')
    error('XLS file not found, create it with no non-header rows to start anew');
end

xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols); 
indexRowN = []; % this is returned if no data found

%compile all dates previously in structure
dates = [];
if isfield(dS, 'date')
    for id = 1:length(dS)
        dates = [dates dS(id).date];
    end
end

%search for unadded dates
for iD = 1:xd.nRows
    dateStr = xd.DateStr{iD};
    if or(sum(strcmp(dates, dateStr))==0, isempty(dates))
        % found an index field with no fit
        indexRowN = iD;
        break;
    else
        % data is already there
        continue
    end
end
