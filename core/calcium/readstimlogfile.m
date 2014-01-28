function slog=readstimlogfile(fn)

fid = fopen(fn);

% discard first five lines
for iline = 1:5
    tline = fgetl(fid);
%    fprintf(1,[tline '\n']);
end

blocksiz = 25000;
blockinc = 10000;
tlog = cell(blocksiz,10);
nlines = 1;

% process lines until blank line
while true
    line = fgetl(fid);
	if isempty(line);break;end;

    if nlines > blocksiz
        tlog(end+1:end+blockinc,:)=cell(blockinc,10);
        blocksiz = blocksiz + blockinc;
    end
    
    tlog(nlines,:)=textscan(line,'%d32%s%d32%d32%d32%d32%d32%d32%d32%s',2,'delimiter','\t');
    nlines = nlines +1;
end

tlog(nlines:end,:)=[];nlines = nlines -1;

fields = {'Trial','EventType','Code','Time',...
          'TTime','TimeUncertainty','Duration','DurationUncertainty','ReqTime','ReqDur'};

slog = cell2struct(tlog,fields,2);

fclose(fid);

return;
