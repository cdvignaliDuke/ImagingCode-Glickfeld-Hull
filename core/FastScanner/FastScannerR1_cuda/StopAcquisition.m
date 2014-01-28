function p = StopAcquisition
p =0;
global fs;

set(fs.handles.lblProgress,'String','Stopping...');
%refresh;

StopDAQ;
StopStream;

if fs.DAQ.adjusting==0;
    closeShutter;
end

set(fs.handles.lblProgress,'String','Idle');
set(fs.handles.btnStopAquisition,'BackgroundColor',.9*[1 1 1]);
set(fs.handles.btnStream,'BackgroundColor',get(fs.handles.uipanel1,'BackgroundColor'));
set(fs.handles.btnScan,'BackgroundColor',get(fs.handles.uipanel1,'BackgroundColor'));
set(fs.handles.btnStream,'Enable','on');
set(fs.handles.btnScan,'Enable','on');
set(fs.handles.btnStartCycles,'Enable','on');


p = 1;

return;
