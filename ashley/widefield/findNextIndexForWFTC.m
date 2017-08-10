function [xd indexRowN] = findNextIndexForWFTC(rc); %(varargin)
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
for iD = 1:xd.nRows
    subjNum = xd.Subject(iD);
    dateStr = xd.DateStr{iD};
    if ~exist(fullfile(rc.structOutput, ['i' num2str(subjNum) '_' dateStr '_tc.mat']),'file')
        % found an index field with no fit
        indexRowN = iD;
        break;
    else
        % data is already there
        continue
    end
end
