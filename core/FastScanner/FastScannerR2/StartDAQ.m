function p=StartDAQ
%%
p=0;
global fs;

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
%next three lines are in SetupDAQ.m, here for information only
%set(fs.DAQ.ai, 'TriggerType', 'HwDigital');  % trigged from PFI0, source AO PCI-6117
%set(fs.DAQ.ai, 'TriggerRepeat', 0);
%set(fs.DAQ.ai, 'TriggerFcn', {'StartStream'});

%% start A0 - Pcell modulation(does not work yet) and trigger for AI
start(fs.DAQ.aoPockels);

fs.DAQ.started = 1;
p = 1;
