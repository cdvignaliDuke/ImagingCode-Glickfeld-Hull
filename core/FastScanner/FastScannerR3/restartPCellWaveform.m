function restartPCellWaveform
global fs;
%warning off all
% events = fs.DAQ.aoPockels.EventLog;
% {events.Type}
% trigdata = events(1).Data
% trigdata = events(end).Data
% get(fs.DAQ.aoPockels,'TriggerType')
% get(fs.DAQ.aoPockels,'TriggersExecuted')
% set(fs.DAQ.aoPockels,'TriggerType')
%%showdaqevents(fs.DAQ.aoPockels)
% to block following warning:
%Warning: NI-DAQ: Samples put to device are less then fifo size.
%Continuous operation may be impossible.

warning off daq:daqmex:propertyConfigurationError 
%stop(fs.DAQ.aoPockels);
putdata(fs.DAQ.aoPockels,[fs.pCell.waveform fs.pCell.TTLout]);
%s=warning('query','last');
% %s
%warning on all
%set(fs.DAQ.aoPockels, 'TriggerType', 'HwDigital'); % % trigged from PFI6, source AO PCI-6117
%set(fs.DAQ.aoPockels, 'TriggerType', 'Manual'); 

start(fs.DAQ.aoPockels);
fs.DAQ.pCellCounter=fs.DAQ.pCellCounter+50;

