function [p msg]= RunAcquisition

global fs;
p=0;
msg='';
%% where to save the data
if ~fs.DAQ.DoNotSave
    %% check for existing data file
    if fs.DAQ.DoNotSave==0 
        k = strfind(fs.DAQ.BaseFileName, '\');
        if isempty(k)
            msg='Check data file name';
            display('Check data file name');
            return
        else
            [s,mess,messid]=mkdir(fs.DAQ.BaseFileName);
            if s~=1                
                msg='Error creating data directory';
                display(mess);
                return;
            end

            fs.DAQ.BaseFileNameLong=[fs.DAQ.BaseFileName fs.DAQ.BaseFileName(k(end):end)];
        end
        fnm=[fs.DAQ.BaseFileNameLong '*.dad'];
        dfiles=dir(fnm);
        if ~isempty(dfiles)
            button = questdlg('Data file already exists. Do you want to continue','Attention','Yes','No','Yes');
            switch button
                case 'Yes'
                case 'No'
                    return
            end
        end
    end % if fs.DAQ.DoNotSave==0
end

fs.fid_curr=[];
fs.fid_prev=[];

%% initialize counters
fs.DAQ.framecounter=0;
fs.DAQ.filecounter=0;
fs.DAQ.clockcounter=0;

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
end

%% check if visual stimulator ready
if getvalue(fs.DAQ.stopLine) == 0
    disp('Stop line low. Automatic stopping enabled.');
    fs.DAQ.stopLineEnabled = 1; % should be low upon startup
else
    fs.DAQ.stopLineEnabled = 0; % line not attached or stimulator not ready
	disp('Stop line high. Automatic stopping disabled.');
end

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