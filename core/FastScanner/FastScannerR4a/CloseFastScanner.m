function p=CloseFastScanner

global fs;

closeShutter

try
    if ~isempty(fs.img.Ch(2).h)
        close(fs.img.Ch(2).h);
        fs.img.Ch(2).h=[];
    end
    if ~isempty(fs.img.Ch(1).h)
        close(fs.img.Ch(1).h);
        fs.img.Ch(1).h=[];
    end
catch
end

putsample(fs.DAQ.aoXYMirrors,[0 0]); %zero scan

% in case objects were not stopped graciously
if fs.DAQ.started
    stop(fs.DAQ.ai);
    stop(fs.DAQ.aoAcqTrigger);
    stop(fs.DAQ.aoPockels);
end

DestroyDAQ;

if ~isempty(fs.stage.comport)
    fclose(fs.stage.comport);
    fs.stage.compor=[];
end

close all;
clear all;
button = questdlg('Do you want to close Matlab?','Exit','Yes','No','No');
if strcmp(button,'Yes')
    exit
end
return;
