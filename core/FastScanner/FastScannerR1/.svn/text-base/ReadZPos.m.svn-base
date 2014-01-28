function pos=ReadZPos
pos=[];
global fs;
if isempty(fs.stage.comport)
    return
end
set(fs.handles.lblStatus,'String',['ESP300 ReadZPos: reading position']);


%flash input
%%
n=get(fs.stage.comport,'BytesAvailable');
if  n > 0
    temp=fread(fs.stage.comport,n);
end

fprintf(fs.stage.comport, '1TP');          
pause(0.1);
n=get(fs.stage.comport,'BytesAvailable');
if n>0
a=fread(fs.stage.comport,n)';
a(find(a>57))=32;
%pos=str2num(char(a))/1000;
%pos=str2num(char(a));
pos=round(str2num(char(a))*1000); %if pos is in um
end
set(fs.handles.lblStatus,'String',['ESP300 ReadZPos: position=' num2str(pos)]);
