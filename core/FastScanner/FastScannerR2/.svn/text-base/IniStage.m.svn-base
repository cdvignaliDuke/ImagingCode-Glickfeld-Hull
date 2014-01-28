function p=IniStage
global fs;
p=0;
fs.stage.ready=0;
% check communication with stage, read current position
if ~isempty(fs.stage.comport)
    fclose(fs.stage.comport);
    pause(0.1);
end


fs.stage.comport = serial(fs.stage.comportname);
set(fs.stage.comport, 'Parity', 'none' , 'Terminator', 'CR', 'StopBits', 1);
set(fs.stage.comport,'DataBits',8); 
set(fs.stage.comport,'FlowControl','hardware');


% open and check status
fopen(fs.stage.comport);
stat=get(fs.stage.comport, 'Status');
if ~strcmp(stat, 'open')
    set(fs.handles.lblStatus,'String',[' MM3000Config: trouble opening port; cannot to proceed']);
    return;
else
    set(fs.handles.lblStatus,'String',[' MM3000Config: port opend']);
end

pause(0.05);
fwrite(fs.stage.comport, '1US1um'); 
fwrite(fs.stage.comport, [13]); 

pause(0.05);

set(fs.handles.lblStatus,'String',[' Stage resolution: 1um']);

pos=ReadZPos;
set(fs.handles.lblCurrentZ,'String',num2str(pos));
fs.position.z=pos;
set(fs.handles.txtZPos,'String',num2str(fs.position.z));

if fs.position.z<get(fs.handles.sliderZPos,'Min')
    set(fs.handles.sliderZPos,'Min',fs.position.z);
end
if fs.position.z>get(fs.handles.sliderZPos,'Max')
    set(fs.handles.sliderZPos,'Max',fs.position.z);
end
set(fs.handles.sliderZPos,'Value',fs.position.z);

fs.stage.ready=1;

p=1;