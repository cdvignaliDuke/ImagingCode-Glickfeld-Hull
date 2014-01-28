function p=GetScannerParameters
%function will start aquisition on hardware trigger and measure
%phase and period
p=0;

global fs;
fs.cycle.Do=0;
fs.DAQ.DoNotSave=1;
fs.DAQ.Running=1;
RunAcquisition;

set(fs.handles.lblProgress,'String','Getting Scanner Parameters');
set(fs.handles.btnStopAquisition,'BackgroundColor',get(fs.handles.uipanel1,'BackgroundColor'));

p=1;