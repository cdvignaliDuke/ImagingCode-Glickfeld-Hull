function pos=SetPositionZ_no_wait(po)
global fs;

if isempty(fs.stage.comport)
    return;
end

set(fs.handles.lblStatus,'String',['ESP300 SetPositionZ_no_wait: setting position ' num2str(po)]);


%flash input
n=get(fs.stage.comport,'BytesAvailable');
if  n > 0
    temp=fread(fs.stage.comport,n);
end

fwrite(fs.stage.comport, '1PA'); 
%fwrite(fs.stage.comport, num2str(po*1000)); 
%fwrite(fs.stage.comport, num2str(po));
fwrite(fs.stage.comport, num2str(po/1000)); %if po is in um
fwrite(fs.stage.comport, [13]); 

pause(0.1);
%{
pos=ReadZPos;
    if ~isempty(pos)
        set(fs.handles.lblCurrentZ,'String',num2str(round(pos)));
%        set(fs.handles.sliderZPos,'Value',pos);
    end
%}