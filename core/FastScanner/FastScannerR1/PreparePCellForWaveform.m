function PreparePCellForWaveform
global fs

if isempty(fs.DAQ.linestart) || isempty(fs.DAQ.period)
    fs.DAQ.linestart=10;
    fs.DAQ.period=800;
end

if fs.iniDone
    if isempty(fs.pCell.LUT)
       voltMax=fs.pCell.max;
       voltMin=fs.pCell.min;
    else
       voltMax=Power2PockelsVoltage(fs.pCell.max);
       voltMin=Power2PockelsVoltage(fs.pCell.min);
    end
end

linestartOut=fs.DAQ.linestart*fs.DAQ.realoutputRate/fs.DAQ.realinputRate;
periodOut=fs.DAQ.period*fs.DAQ.realoutputRate/fs.DAQ.realinputRate;
%oneFrameLe=round((120+128*50)*periodOut); %50 frames without last lightback
oneFrameLe=round((120+128)*periodOut); %1 frames without last lightback

blockFraction=0.1;

outInd=1:oneFrameLe;
outInd=outInd-(linestartOut-blockFraction*periodOut/2);


fs.pCell.waveform=ones(oneFrameLe,1)*voltMax;

%if no blocking required then waveform is flat voltMax
if fs.DAQ.BlockEnds==1
    fs.pCell.waveform(mod(outInd,periodOut/2)<blockFraction*periodOut)=voltMin;
    fs.pCell.waveform(end+1)=voltMin;
else
    fs.pCell.waveform(end+1)=voltMax;
end

%fs.pCell.TTLout(1:fs.DAQ.realoutputRate/1000)=5; %1 ms TTL puls
%initial state of TTLout is 5v, here we just set it to zero
%fs.pCell.TTLout=zeros(oneFrameLe+1,1);

if floor(oneFrameLe/2)==oneFrameLe/2    
%     % trigger is a half-a-frame long pulse
%     fs.pCell.TTLout=[fs.DAQ.aoAcqTriggerStreamValue*ones(floor(oneFrameLe/2),1);...
%                      fs.DAQ.aoAcqTriggerIdleValue*ones(floor(oneFrameLe/2)+1,1)]; 
                 
    % trigger is a level transition 
    fs.pCell.TTLout=fs.DAQ.aoAcqTriggerStreamValue*...
        [ones(floor(oneFrameLe/2),1);ones(floor(oneFrameLe/2)+1,1)]; 
else
%     % trigger is a half-a-frame long pulse
%     fs.pCell.TTLout=[fs.DAQ.aoAcqTriggerStreamValue*ones(floor(oneFrameLe/2),1);...
%                      fs.DAQ.aoAcqTriggerIdleValue*ones(floor(oneFrameLe/2)+2,1)]; 
                 
    % trigger is a level transition 
    fs.pCell.TTLout=fs.DAQ.aoAcqTriggerStreamValue*...
        [ones(floor(oneFrameLe/2),1);ones(floor(oneFrameLe/2)+2,1)]; 
end

set(fs.DAQ.aoPockels,'BufferingConfig',[200 800]);

putdata(fs.DAQ.aoPockels,[fs.pCell.waveform fs.pCell.TTLout]);

%fs.pCell.TTLout(:)=0; %reset TTL out to create single TTL puls

%%next to lines - attempt to restart waveform in 'SamplesOutputFcn'
%set(fs.DAQ.aoPockels,'SamplesOutputFcn','restartPCellWaveform');
%set(fs.DAQ.aoPockels,'SamplesOutputFcnCount',oneFrameLe+1);



