function p=ReadIniForFS(IniFName,gui)
global fs;
global fd;
p=0;

if nargin < 2 
    gui = 'FastScanner';
end

ReadIniCommon(IniFName);

switch gui
    case 'FastScanner'
        ReadIniFastScanner(IniFName);
        UpdateFastScannerByParam;
    case 'FastDisplay'
        ReadIniFastDisplay(IniFName);
        UpdateFastDisplayByParam;
end

return;

%%
function ReadIniCommon(IniFName);
global fs;
global fd;

retString=GetProfileString(IniFName,'DAQ','MemMapFile');
if isempty(retString)
    fs.DAQ.MemMapFile='C:\MemMap.dat';
    fd.DAQ.MemMapFile='C:\MemMap.dat';
else
    fs.DAQ.MemMapFile=retString;
    fd.DAQ.MemMapFile=retString;
end

return;

%%
function ReadIniFastScanner(IniFName);
global fs;
retString=GetProfileString(IniFName,'pCell','min');
if isempty(retString)
    fs.pCell.min=0;
else
    fs.pCell.min=str2num(retString);
end

retString=GetProfileString(IniFName,'pCell','max');
if isempty(retString)
    fs.pCell.max=0;
else
    fs.pCell.max=str2num(retString);
end


retString=GetProfileString(IniFName,'pCell','LUTfile');
if isempty(retString)
    fs.pCell.LUTfile='C:\projects\fastscanner\PockelsCellLUT.txt';
else
    fs.pCell.LUTfile=retString;
end

retString=GetProfileString(IniFName,'pCell','LUTdir');
if isempty(retString)
    fs.pCell.LUTdir='C:\projects\fastscanner';
else
    fs.pCell.LUTdir=retString;
end
files = dir(fs.pCell.LUTdir);
ind=1;
for i=1:length(files)
     if strcmp(fs.pCell.LUTfile, fullfile(fs.pCell.LUTdir,files(i).name))
        ind=i;
    end    
end
set(fs.handles.lstboxCalibFiles,'String',{files(3:end).name});
if ind>2
    set(fs.handles.lstboxCalibFiles,'Value',ind-2);
else
    set(fs.handles.lstboxCalibFiles,'Value',1);
    contents = get(hObject,'String');
    fs.pCell.LUTfile=fullfile(fs.pCell.LUTdir,contents{get(hObject,'Value')});
end



retString=GetProfileString(IniFName,'position','zoom');
if isempty(retString)
    fs.position.zoom=1;
else
    fs.position.zoom=str2num(retString);
end

retString=GetProfileString(IniFName,'stage','comportname');
if isempty(retString)
    fs.stage.comportname='COM1';
else
    fs.stage.comportname=retString;
end


fs.position.scale=1;


%{
% TODO - read from ini folowing parameters
fs.DAQ.acquisitionBoardIndex=2;
fs.DAQ.inputRate=10000000;
fs.DAQ.PMTChannelIndex(1)=0;
fs.DAQ.PMTChannelIndex(2)=1;
fs.DAQ.XcoordChannelIndex=2;

fs.DAQ.mirrorOutputBoardIndex=2;
fs.DAQ.outputRate=40000;
fs.DAQ.XMirrorChannelIndex=0;
fs.DAQ.YMirrorChannelIndex=1;

fs.DAQ.secondaryBoardIndex=1;
fs.DAQ.pockelsChannelIndex=0;
fs.DAQ.aoTriggerChannelIndex1;

fs.DAQ.triggerBoardIndex=2;
fs.DAQ.triggerLineIndex=1;
fs.DAQ.shutterLineIndex=2;

fs.shutter.closedLevel=1;
%}

return;

%%
function ReadIniFastDisplay(IniFName);

global fd;

retString=GetProfileString(IniFName,'img','Ch1Use');
if isempty(retString)
    fd.img.Ch(1).use=1;
else
    fd.img.Ch(1).use=str2num(retString);
end

retString=GetProfileString(IniFName,'img','Ch2Use');
if isempty(retString)
    fd.img.Ch(2).use=1;
else
    fd.img.Ch(2).use=str2num(retString);
end

retString=GetProfileString(IniFName,'img','XSizePix');
if isempty(retString)
    fd.img.XSizePix=256;
else
    fd.img.XSizePix=str2num(retString);
end

retString=GetProfileString(IniFName,'img','YSizePix');
if isempty(retString)
    fd.img.YSizePix=220;
else
    fd.img.YSizePix=str2num(retString);
end

retString=GetProfileString(IniFName,'timing','FastMirrorTriggerDelay');
if isempty(retString)
    fd.timing.FastMirrorTriggerDelay=0;
else
    fd.timing.FastMirrorTriggerDelay=str2num(retString);
end


return;
