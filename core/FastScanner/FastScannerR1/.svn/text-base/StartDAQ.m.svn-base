function p=StartDAQ
%%
p=0;
global fs;

%fs.DAQ.framecounter = 0; %moved to RunAcquisition
fs.DAQ.framecounter_thisacq = 0;

%% load waveform to PCell channel
SetPCell;
%PreparePCellForWaveform prepares and loads waveformes for PCell and "TTL" puls
PreparePCellForWaveform;
if fs.DAQ.BlockEnds==1    
	set(fs.handles.sliderPCellMin,'Enable','off');
    set(fs.handles.sliderPCellMax,'Enable','off'); 
	set(fs.handles.txtPCellMin,'Enable','off');
    set(fs.handles.txtPCellMax,'Enable','off');     
end

fs.DAQ.pCellCounter=0;
    
%% acquisition object triggered by aoAcqTrigger
start(fs.DAQ.ai);

%% start A0 - Pcell modulation and trigger for AI
start(fs.DAQ.aoPockels);

fs.DAQ.started = 1;
p = 1;
