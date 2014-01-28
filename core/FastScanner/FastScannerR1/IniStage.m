function p=IniStage
global fs;
p=0;
fs.stage.ready=0;

%%

%% close old ports
% ip = instrfind;
% if ~isempty(ip)
%     fclose(ip);
% end


%% init
if isempty(fs.stage.comport)
 set(fs.handles.lblStatus,'String',['ESP300 IniStage: start creating port']);
drawnow   
pause(0.1);
fs.stage.comport = serial(fs.stage.comportname);
set(fs.stage.comport, 'Parity', 'none' , 'Terminator', 'CR', 'StopBits', 1);
set(fs.stage.comport,'DataBits',8); 
set(fs.stage.comport,'FlowControl','hardware');
set(fs.stage.comport,'BaudRate',19200);
else
    set(fs.handles.lblStatus,'String',['ESP300 IniStage: port already exists']);
    drawnow
    pause(0.1);
end


stat=get(fs.stage.comport, 'Status');
if ~strcmp(stat, 'open')
        set(fs.handles.lblStatus,'String',['ESP300 IniStage: oppening port']);
    drawnow
    pause(0.1);
% open and check status
fopen(fs.stage.comport);
stat=get(fs.stage.comport, 'Status');
end

%%
if ~strcmp(stat, 'open')
    set(fs.handles.lblStatus,'String',['ESP300 IniStage: trouble opening port; cannot to proceed']);
    return;
else
    set(fs.handles.lblStatus,'String',['ESP300 IniStage: port opend']);
end
drawnow
pause(0.1);
%%
%{
pause(0.05);
fwrite(fs.stage.comport, '1US1um'); 
fwrite(fs.stage.comport, [13]); 

pause(0.05);

set(fs.handles.lblStatus,'String',[' Stage resolution: 1um']);
%}

%%

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
set(fs.handles.txtStartZ, 'String', num2str(fs.position.z));

fs.stage.ready=1;

if fs.stage.JS==1
    %transfer control to JS
    fwrite(fs.stage.comport, 'EX JOYSTICK ON');
    fwrite(fs.stage.comport, [13]);
else
    %disable JS
    fwrite(fs.stage.comport, 'EX JOYSTICK OFF');
    fwrite(fs.stage.comport, [13]);
end



p=1;