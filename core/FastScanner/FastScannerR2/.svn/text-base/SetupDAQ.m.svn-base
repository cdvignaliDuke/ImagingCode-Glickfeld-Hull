function p=SetupDAQ
% create all channels for in and out

global fs;
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCI-6115
%% analog inputs 
fs.DAQ.ai = analoginput('nidaq',fs.DAQ.acquisitionBoardIndex); % NI PCI-6115 Board
set(fs.DAQ.ai, 'SampleRate', fs.DAQ.inputRate);
fs.DAQ.realinputRate=get(fs.DAQ.ai, 'SampleRate');
set(fs.DAQ.ai, 'TriggerType', 'HwDigital');  % trigged from PFI0, source AO PCI-6117
set(fs.DAQ.ai, 'TriggerRepeat', 0);
set(fs.DAQ.ai, 'TriggerFcn', {'StartStream'});
rate=fs.DAQ.realinputRate

%dma=get(fs.DAQ.ai, 'TransferMode') %SingleDMA
fs.ParamMap = fopen(fs.DAQ.ParamMap,'w');
fwrite(fs.ParamMap,fs.DAQ.realinputRate,'double');
fclose(fs.ParamMap);

% get input ranges and load them into list box
DAQinfo= daqhwinfo(fs.DAQ.ai);
fs.DAQ.InputRanges=DAQinfo.InputRanges;
set(fs.handles.lstRanges,'String',num2str(fs.DAQ.InputRanges));

%two PMT channels
fs.DAQ.aiChanPMT1=addchannel(fs.DAQ.ai,fs.DAQ.PMTChannelIndex(1));
fs.DAQ.aiChanPMT2=addchannel(fs.DAQ.ai,fs.DAQ.PMTChannelIndex(2));

%fast mirror - X coordinate channel
fs.DAQ.aiChanXPos=addchannel(fs.DAQ.ai,fs.DAQ.XcoordChannelIndex);
fs.DAQ.aiChanYPos=addchannel(fs.DAQ.ai,fs.DAQ.YcoordChannelIndex);
fs.DAQ.aiChanXPos.Coupling='AC';
fs.DAQ.aiChanXPos.Coupling='DC';

% acceptable ranges +/- 5.0 V, 2.0 V, 1.0 V, 0.5 V, and 0.2 V
fs.DAQ.ai.Channel(1:end).InputRange=[-5 5]; 
fs.DAQ.ai.Channel(1:2).InputRange=fs.DAQ.InputRanges(1,:);

%specify DAQ paramters (number of points to collect)
fs.DAQ.InPointsToCollect=round((256-16)*(fs.DAQ.realinputRate*fs.timing.FastMirrorPeriod/2000));

frame_time=fs.DAQ.InPointsToCollect/fs.DAQ.realinputRate %should directly correspond to MirrorPeriod

set(fs.DAQ.ai,'BufferingConfig',[fs.DAQ.InPointsToCollect 196]);
set(fs.DAQ.ai, 'SamplesPerTrigger', inf);
set(fs.DAQ.ai, 'Timeout', 10);
% set(fs.DAQ.ai,'DataMissedFcn','HandleDAQError');
% set(fs.DAQ.ai,'RuntimeErrorFcn','HandleDAQError');


%{
%% external trigger, this signals to stimulator that acquisition is ongoing
%% 6115
fs.DAQ.aoExtTrigger = analogoutput('nidaq',fs.DAQ.acquisitionBoardIndex);
set(fs.DAQ.aoExtTrigger, 'SampleRate', fs.DAQ.outputRate);
fs.DAQ.aoChanExtTrigger= addchannel(fs.DAQ.aoExtTrigger,...
                                    fs.DAQ.aoExtTriggerChannelIndex,'ExternalTrigger');
set(fs.DAQ.aoExtTrigger, 'TriggerType', 'HwDigital');
%set(fs.DAQ.aoExtTrigger, 'HwDigitalTriggerSource','PFI6');
set(fs.DAQ.aoExtTrigger, 'TriggerFcn', {'FastScannerCallback'});
% set(fs.DAQ.aoExtTrigger,'RuntimeErrorFcn','HandleDAQError');
putsample(fs.DAQ.aoExtTrigger,[5]); % receiver triggers on negative edge
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCI-6711
%% pockell cell and acuisition trigger (combined in one AO in V6)
fs.DAQ.aoPockels = analogoutput('nidaq',fs.DAQ.secondaryBoardIndex);
set(fs.DAQ.aoPockels, 'SampleRate', fs.DAQ.outputRate);
fs.DAQ.realoutputRate=get(fs.DAQ.aoPockels, 'SampleRate');
fs.DAQ.aoChanPockels = addchannel(fs.DAQ.aoPockels, fs.DAQ.pockelsChannelIndex);
fs.DAQ.aoChanAcqTrigger= addchannel(fs.DAQ.aoPockels, fs.DAQ.aoAcqTriggerChannelIndex,'AcquisitionTrigger');
set(fs.DAQ.aoPockels, 'TriggerType', 'HwDigital'); % default PFI6
%set(fs.DAQ.aoPockels, 'TriggerFcn', {'FastScannerCallback'});
set(fs.DAQ.aoPockels,'BufferingConfig',fs.DAQ.aoAcqTriggerBufferConfig);
fs.DAQ.aoCurrentTriggerValue=5;
putsample(fs.DAQ.aoPockels,[0 fs.DAQ.aoCurrentTriggerValue]); % receiver triggers on negative edge, set trigger channel to 0V

% set(fs.DAQ.aoPockels,'RuntimeErrorFcn','HandleDAQError');
% syncNIDAQBoards(fs.DAQ.aoExtTrigger, fs.DAQ.aoPockels);	%  sync the 2 board clocks.

%% X mirror, Y mirror
% analog output, mirror control voltages
fs.DAQ.aoXYMirrors = analogoutput('nidaq',fs.DAQ.secondaryBoardIndex);
set(fs.DAQ.aoXYMirrors, 'SampleRate', fs.DAQ.outputRate);
fs.DAQ.aoChanXMirror = addchannel(fs.DAQ.aoXYMirrors, fs.DAQ.XMirrorChannelIndex);
fs.DAQ.aoChanYMirror= addchannel(fs.DAQ.aoXYMirrors, fs.DAQ.YMirrorChannelIndex);
% set(fs.DAQ.aoXYMirrors,'RuntimeErrorFcn','HandleDAQError');

%{
%% acquisition trigger, this signal triggers start of acquisition % 6711
fs.DAQ.aoAcqTrigger = analogoutput('nidaq',fs.DAQ.secondaryBoardIndex);
set(fs.DAQ.aoAcqTrigger, 'SampleRate', fs.DAQ.outputRate);
fs.DAQ.aoChanAcqTrigger= addchannel(fs.DAQ.aoAcqTrigger, ...
    fs.DAQ.aoAcqTriggerChannelIndex,'AcquisitionTrigger');
set(fs.DAQ.aoAcqTrigger, 'TriggerType', 'HwDigital'); % default PFI6
% set(fs.DAQ.aoAcqTrigger, 'HwDigitalTriggerSource','PFI0');
set(fs.DAQ.aoAcqTrigger, 'TriggerFcn', {'FastScannerCallback'});
set(fs.DAQ.aoAcqTrigger,'BufferingConfig',fs.DAQ.aoAcqTriggerBufferConfig);
% set(fs.DAQ.aoAcqTrigger,'RuntimeErrorFcn','HandleDAQError');
putsample(fs.DAQ.aoAcqTrigger,[5]); % receiver triggers on negative edge
%}
%% digital outputs
fs.DAQ.dio = digitalio('nidaq', fs.DAQ.triggerBoardIndex);
% trigger line (deprecated)
fs.DAQ.triggerLine = addline(fs.DAQ.dio, fs.DAQ.triggerLineIndex, 'out');
putvalue(fs.DAQ.triggerLine, 0); 
% shutter digital out line
fs.DAQ.shutterLine = addline(fs.DAQ.dio, fs.DAQ.shutterLineIndex, 'out');
closeShutter;
fs.DAQ.stopLine = addline(fs.DAQ.dio, fs.DAQ.stopLineIndex,'in');

%%
fs.DAQ.created=1;