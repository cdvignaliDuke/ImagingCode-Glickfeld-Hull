function pos=SetPositionZ_no_wait(po)
global fs;

if isempty(fs.stage.comport)
    return;
end


%pos=ReadZPos;
%pause(0.05);
fwrite(fs.stage.comport, '1PA'); 
%fwrite(fs.stage.comport, num2str(po*1000)); 
fwrite(fs.stage.comport, num2str(po)); 
fwrite(fs.stage.comport, [13]); 

pause(0.1);

pos=ReadZPos;
%po
%pos
    if ~isempty(pos)
        set(fs.handles.lblCurrentZ,'String',num2str(pos));
%        set(fs.handles.sliderZPos,'Value',pos);
    end

