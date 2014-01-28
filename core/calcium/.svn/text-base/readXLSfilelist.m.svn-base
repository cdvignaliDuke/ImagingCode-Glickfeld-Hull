function s = readXLSfilelist (xlsfname)

% read file list from xls file (Excel 97 or 2003 format)
% xls file should be written as follows:
%
% name1     name2   name3   name4
% xxx       xxx     xxx     xxx
% xxx       xxx     xxx     xxx
%
% xxx could be either numbers of letters
%
% 7/17/08 Kenichi Ohki

[num,text,raw]=xlsread(xlsfname);
[Nlines, Nparams]=size(raw)
s=struct([]);
for param=1:Nparams
    fields{param}=char(raw(1,param));
end
s=cell2struct(raw(2:Nlines,:), fields, 2);

% if the first column is empty, delete the lines
for lines=1:Nlines-1
    if isnan(getfield(s,{lines},char(raw(1,1))))
        s(lines:Nlines-1,:)=[];
        break;
    end
end

% convert text to number
N=length(s);
for i=1:N
    for param=1:Nparams
        temp=getfield(s,{i},fields{param});
        if ischar(temp)
            if ~isempty(str2num(temp))
                s=setfield(s,{i},fields{param},str2num(temp));
            end
        end
    end
end

