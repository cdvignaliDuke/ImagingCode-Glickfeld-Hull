function writeXLSfilelist (s, xlsfname)

% write file list to xls file (Excel 97 or 2003 format)
% s: structure array
% xls file will be written as follows:
%
% name1     name2   name3   name4
% xxx       xxx     xxx     xxx
% xxx       xxx     xxx     xxx
%
% xxx could be either numbers of letters
%
% 7/17/08 Kenichi Ohki

data=struct2cell(s)';
names=fieldnames(s)';

xlswrite(xlsfname,[names;data]);