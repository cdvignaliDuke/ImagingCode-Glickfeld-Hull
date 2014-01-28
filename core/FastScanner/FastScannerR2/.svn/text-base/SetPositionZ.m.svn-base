function pos=SetPositionZ(po)
global fs;

if isempty(fs.stage.comport)
    return;
end


pos=ReadZPos;
pause(0.05);
fwrite(fs.stage.comport, '1PA'); 
%fwrite(fs.stage.comport, num2str(po*1000)); 
fwrite(fs.stage.comport, num2str(po)); 
fwrite(fs.stage.comport, [13]); 

pause(0.05);

pos=ReadZPos;

n=0;
while abs(pos-po)>0.002
    pause(0.02);
    pos=ReadZPos;
    if ~isempty(pos)
        set(fs.handles.lblCurrentZ,'String',num2str(pos));
    end
    n=n+1;
    if n>50
        set(fs.handles.lblStatus,'String',['SetPositionZ: was not able get to '  num2str(po) ' mm']);
        break
    end
end

