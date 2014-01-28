function [p msg]= RunAcquisition

global fs;
p=0;
msg='';
%% where to save the data
if ~fs.DAQ.DoNotSave
    if fs.DAQ.DoNotSave==0 
        % create base directory
        [datPath baseDir] = fileparts(fs.DAQ.BaseFileName);

        % data path with dad file name but without dad or number extension
        fs.DAQ.BaseFileNameLong = fullfile(fs.DAQ.BaseFileName, baseDir);
        fs.DAQ.BaseFileNameLongGreen = [fullfile(fs.DAQ.BaseFileName, baseDir) '_green\' baseDir];
        fs.DAQ.BaseFileNameLongRed = [fullfile(fs.DAQ.BaseFileName, baseDir) '_red\' baseDir];
        
        if isempty(baseDir)
            errordlg('Missing last directory, check data path');
            return
        else
            [success, message] = mkdir([fullfile(fs.DAQ.BaseFileName, baseDir) '_green\']);
            if ~success
                errordlg('Error creating data directory');
                return
            end
            [success, message] = mkdir([fullfile(fs.DAQ.BaseFileName, baseDir) '_red\']);
            if ~success
                errordlg('Error creating data directory');
                return
            end            
        end        
        
        % check if any dad files are there, prompt and then delete
        fnm=[fs.DAQ.BaseFileNameLongGreen '*.tif'];
        dfiles=dir(fnm);
        if ~isempty(dfiles)
            button = questdlg('Data file already exists. Do you want to continue','Attention','Yes','No','Yes');
            switch button
                case 'Yes'
                    w = warning('off', 'MATLAB:DELETE:Permission');  % some may be open from last run, too hard to deal, ignore
                    delete(fullfile(fs.DAQ.BaseFileName, '*.tif'));  % remove all dad files
                    warning(w);
                case 'No'
                    return
            end
        end
    end % if fs.DAQ.DoNotSave==0
end

fs.fid_curr=[];
fs.fid_prev=[];

fs.timing.runStartClock = clock;

%% initialize counters
fs.DAQ.framecounter=0;
fs.DAQ.filecounter=0;
fs.DAQ.clockcounter=0;

%% initialize z stack counters
fs.DAQ.zIsWaiting = false;
fs.DAQ.zNextStartTime = NaN;
fs.DAQ.zClosePrevDadTime = NaN;
fs.DAQ.framecounter_thisacq = 0;
fs.DAQ.dad_size_in_frames = 150; 
%% initialize data buffer and images

fs.datasize1=240;

if fs.DAQ.realinputRate<4000000
    fs.datasize2=256;
else
    fs.datasize2=512;    
end


fs.DAQ.data_leftover=[];
fs.DAQ.fgreen=zeros(fs.datasize1,fs.datasize2,'uint16')*2048;
fs.DAQ.fred=zeros(fs.datasize1,fs.datasize2,'uint16')*2048;
fs.DAQ.InPointsToCollect=round(1.5*(256)*(fs.DAQ.realinputRate*fs.timing.FastMirrorPeriod/2000));

%% debug information
fs.tictoc=[];
fs.sAcquired=[];
fs.sAvailable=[];
fs.totaltime=[];
fs.totaldelta=[];
fs.sAvailableFrames=[];
fs.phase=[];

%% initialize stage if needed
if fs.cycle.Do
    if fs.stage.ready
        if fs.cycle.NFramesPrStep < 50
            errordlg('Frames Per Step must be > 50 during Z stack; avoids too frequent acq restart');
            return
        end
        if fs.stage.JS==1
            errordlg('Disable JS before start');
            return
        end
        
        %assume that start Z position is current position
        pos=ReadZPos;
        fs.position.z=pos;
        fs.cycle.ZStart=pos;
        set(fs.handles.txtZPos,'String',num2str(fs.position.z));
        set(fs.handles.txtStartZ,'String',num2str(fs.cycle.ZStart));
        if fs.position.z<get(fs.handles.sliderZPos,'Min')
            set(fs.handles.sliderZPos,'Min',fs.position.z);
        end
        if fs.position.z>get(fs.handles.sliderZPos,'Max')
            set(fs.handles.sliderZPos,'Max',fs.position.z);
        end
        set(fs.handles.sliderZPos,'Value',fs.position.z);
    end
%    fs.DAQ.dad_size_in_frames= fs.cycle.NFramesPrStep;
end

%% check if visual stimulator ready
if getvalue(fs.DAQ.stopLine) == 0
    disp('Stop line low. Automatic stopping enabled.');
    fs.DAQ.stopLineEnabled = 1; % should be low upon startup
else
    fs.DAQ.stopLineEnabled = 0; % line not attached or stimulator not ready
	disp('Stop line high. Automatic stopping disabled.');
end

%% write log file
fs.DAQ.LogName = fullfile(fs.DAQ.BaseFileName, '00fastscan-data-log.txt');
if ~fs.DAQ.DoNotSave
    fs.DAQ.fdLog = fopen(fs.DAQ.LogName, 'wt');
    fd1 = fs.DAQ.fdLog;
    fprintf(fd1, 'FastScanner data log - 00fastscan-data-log.txt\n');
    fprintf(fd1, 'Running from dir: %s\n', pwd);
    fprintf(fd1, 'Current time: %s\n', datestr(now));
    fprintf(fd1, 'fs.DAQ.realinputRate: %d\n', round(fs.DAQ.realinputRate));
    fprintf(fd1, 'fs.DAQ.baseFileName: %s\n', fs.DAQ.BaseFileName);
    fprintf(fd1, 'fs.cycle.ZStart: %d\n', fs.cycle.ZStart);
    fprintf(fd1, 'fs.cycle.ZStep: %d\n', fs.cycle.ZStep);
    fprintf(fd1, 'fs.cycle.NSteps: %d\n', fs.cycle.NSteps);
    fprintf(fd1, 'fs.cycle.NReps: %d\n', fs.cycle.NReps);
    fprintf(fd1, 'fs.cycle.NFramesPrStep: %d\n', fs.cycle.NFramesPrStep);
    fprintf(fd1, 'fs.cycle.Do: %d\n', fs.cycle.Do);
    fprintf(fd1, 'fs.position.zoom: %5.2f\n', fs.position.zoom); 
    fprintf(fd1, 'fs.DAQ.dad_size_in_frames: %d\n', fs.DAQ.dad_size_in_frames);
end

%% start acq
openShutter; pause(.125);
StartDAQ;
if fs.DAQ.DoNotSave
fs.DAQ.modeString='Preview';
else
fs.DAQ.modeString='Stream';
end

if fs.cycle.Do
fs.DAQ.modeString='Stream-stack';
end
p=1

return;