function pos=SetPositionZ(po)
global fs;

if isempty(fs.stage.comport)
    return;
end
%%

%pos=ReadZPos;
%pause(0.05);

set(fs.handles.lblStatus,'String',['ESP300 SetPositionZ: setting position ' num2str(po)]);

%flash input
n=get(fs.stage.comport,'BytesAvailable');
if  n > 0
    temp=fread(fs.stage.comport,n);
end
pause(0.1);

fwrite(fs.stage.comport, '1PA'); 
%fwrite(fs.stage.comport, num2str(po*1000)); 
%fwrite(fs.stage.comport, num2str(po));
fwrite(fs.stage.comport, num2str(po/1000)); %if po is in um
fwrite(fs.stage.comport, [13]); 

pause(0.1);

pos=ReadZPos;

n=0;
while abs(pos-po)>0.002
    pause(0.05);
    pos=ReadZPos;
    if ~isempty(pos)
        set(fs.handles.lblCurrentZ,'String',num2str(round(pos)));
    end
    n=n+1;
    if n>40
        set(fs.handles.lblStatus,'String',['SetPositionZ: was not able get to '  num2str(round(po)) ' um']);
        break
    end
end

