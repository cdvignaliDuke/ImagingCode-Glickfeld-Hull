function s=readPresentationLog(filename)
% READPRESENTATIONLOG Reads and parses presentation log file 
% S=READPRESENTATIONLOG(FILENAME) where S is log file elements as STRUCT
%

fid = fopen(filename);

% discard first 3 lines
tline = fgetl(fid); % 1st line
tline = fgetl(fid); % 2nd
tline = fgetl(fid);

% disambiguate format
tline = fgetl(fid);
c = strtokc(tline);
ntoks = length(c);

switch ntoks
    case 10 % default log format
        fmt = '%d32%s%d32%d32%d32%d32%d32%d32%d32%s';
        fields = {'Trial','EventType','Code','Time',...
                  'TTime','TimeUncertainty','Duration','DurationUncertainty','ReqTime','ReqDur'};
        
    case 13
        fmt = '%d32%s%s%d32%d32%d32%d32%d32%d32%d32%d32%d32%s';
        fields = {'Trial','EventType','Code','iTrial','iStim','iFrame','Time',...
          'TTime','TimeUncertainty','Duration','DurationUncertainty','ReqTime','ReqDur'};
      
    case 14
        fmt = '%d32%s%s%d32%d32%d32%d32%d32%d32%d32%d32%d32%d32%s';
        fields = {'Trial','EventType','Code','iTrial','iStim','bBlank','iFrame','Time',...
          'TTime','TimeUncertainty','Duration','DurationUncertainty','ReqTime','ReqDur'};
      
    otherwise
        fprintf(1,'log file has %i columns\n',ntoks);
        error('unexpected format');
end
    
tline = fgetl(fid);

blocksiz = 25000;
blockinc = 10000;
tlog = cell(blocksiz,ntoks);
nlines = 1;

w = warning('off', 'MATLAB:intConvertNaN');

% process lines until blank line
while true
    line = fgetl(fid);
	if isempty(line);break;end;

    if nlines > blocksiz
        tlog(end+1:end+blockinc,:)=cell(blockinc,ntoks);
        blocksiz = blocksiz + blockinc;
    end
    
    c = textscan(line,fmt,1,'delimiter','\t');

    tlog(nlines,1:length(c))=c;
    nlines = nlines +1;
end

warning(w);

tlog(nlines:end,:)=[];nlines = nlines -1;

s = cell2struct(tlog,fields,2);

fclose(fid);

return;
