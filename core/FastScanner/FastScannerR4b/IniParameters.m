function IniParameters(process)
        global fd;
        global fs;
if nargin < 1
    process = 'FastScanner';
end
inputRate=7000000;
switch process
    case 'FastScanner';

        fs.DAQ.MemMapFile='C:\MemMap.dat';
        fs.DAQ.ParamMap='C:\ParamMap.dat';
        fs.pCell.min=0;
        fs.pCell.max=1;
        fs.pCell.LUTfile='C:\Documents and Settings\user\My Documents\logs\laser\PockelsCellLUT.txt';
        fs.pCell.LUT=[];
        fs.position.x=0;
        fs.position.y=0;
        fs.position.z=0;
        fs.position.zoom=2.5;
        fs.position.scale=5;
        fs.position.XtoYscale=2.5;
        fs.position.YXratio=1;

        fs.cycle.XStart=0;
        fs.cycle.YStart=0;
        fs.cycle.ZStart=0;
        fs.cycle.XStep=0;
        fs.cycle.YStep=0;
        fs.cycle.ZStep=0;
        fs.cycle.NSteps=1;
        fs.cycle.NReps=1;
        fs.cycle.NFramesPrStep=0;
        fs.cycle.Avg=0;
        
        % pci-6115, analog input/outputs
        fs.DAQ.acquisitionBoardIndex=1;
        fs.DAQ.inputRate=inputRate;
        fs.DAQ.PMTChannelIndex(1)=0;
        fs.DAQ.PMTChannelIndex(2)=1;
        fs.DAQ.XcoordChannelIndex=2;
        fs.DAQ.YcoordChannelIndex=3;
        fs.DAQ.aoExtTriggerChannelIndex=0;
        fs.DAQ.outputRate=100000;
        
        % pci-6711 analog outputs only
        fs.DAQ.secondaryBoardIndex=2;
        fs.DAQ.pockelsChannelIndex=0;
        fs.DAQ.XMirrorChannelIndex=1;
        fs.DAQ.YMirrorChannelIndex=2;
        fs.DAQ.aoAcqTriggerChannelIndex = 3;
        fs.DAQ.aoAcqTriggerBufferConfig = [2048 2];

        % digital inputs / outputs
        fs.DAQ.triggerBoardIndex=1;
        fs.DAQ.triggerLineIndex=1;
        fs.DAQ.shutterLineIndex=2;
        fs.DAQ.stopLineIndex=3;
        fs.DAQ.MainshutterLineIndex=4;

        
        fs.DAQ.aiTriggerPolarity = 'PositiveEdge';        
        fs.DAQ.aoAcqTriggerIdleValue = 0;
        fs.DAQ.aoAcqTriggerStreamValue = 6;
        fs.DAQ.aoCurrentTriggerValue=fs.DAQ.aoAcqTriggerIdleValue;
               
        fs.DAQ.BaseFileName='E:\Data\tst';

        fs.DAQ.created=0;
        fs.DAQ.started=0;
        fs.DAQ.TestMode=0;
        fs.DAQ.adjusting=0;
        fs.DAQ.DoNotSave=0;
        
        fs.DAQ.BlockEnds=0;

        fs.shutter.closedLevel=1;
        fs.shutter.MainclosedLevel=1;

        fs.DAQ.framecounter=0;
        fs.DAQ.Running=0;
        fs.DAQ.Streaming=0;
        
        fs.DAQ.linestart=[];
        fs.DAQ.period=[];
        fs.DAQ.delay=20;
        fs.DAQ.calculating_delay=0;

        fs.stage.comport=[];
        fs.stage.comportname='COM1';
        fs.stage.ready=0;
        
        fs.stage.JS=0;

        fs.timing.FastMirrorPeriod=0.2545;

        UpdateFastScannerByParam;
    case 'FastDisplay';

        fd.DAQ.MemMapFile='C:\MemMap.dat';
        fd.DAQ.ParamMap='C:\ParamMap.dat';

        fd.img.Ch(1).min=0;
        fd.img.Ch(1).max=4096;
        fd.img.Ch(2).min=0;
        fd.img.Ch(2).max=4096;

        fd.img.Ch(1).use=1;
        fd.img.Ch(1).h=[];
        fd.img.Ch(1).CData=[];
        fd.img.Ch(1).CData1D=[];

        fd.img.Ch(2).use=0;
        fd.img.Ch(2).h=[];
        fd.img.Ch(2).CData=[];
        fd.img.Ch(2).CData1D=[];

        fd.img.XSizePix=256;
        fd.img.YSizePix=256;

        fd.timing.FastMirrorTriggerDelay=0;

        fd.img.Ch(1).ax=[];
        fd.img.Ch(2).ax=[];
        fd.img.AvgOnDisplay=2;

        fd.Display.DoDisplay=0;
        fd.DisplayTimer=[];
        
        fd.stats.show=0;
        fd.stats.region.x1=1;
        fd.stats.region.y1=1;
        if inputRate<4000000
            fd.stats.region.x2=256;
        else
            fd.stats.region.x2=512;
        end
        fd.stats.region.x2=256;
        fd.stats.region.y2=240;

        UpdateFastDisplayByParam;
end

