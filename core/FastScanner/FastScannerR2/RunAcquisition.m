function p = RunAcquisition

global fs;

%% where to save the data
if ~fs.DAQ.DoNotSave
    %% check for existing data file
    if fs.DAQ.DoNotSave==0 
        k = strfind(fs.DAQ.BaseFileName, '\');
        if isempty(k)
            display('Check data file name');
            return
        else
            [s,mess,messid]=mkdir(fs.DAQ.BaseFileName);
            if s~=1
                display('Error creating data directory');
                return;
            end

            fs.DAQ.BaseFileNameLong=[fs.DAQ.BaseFileName fs.DAQ.BaseFileName(k(end):end)];
        end
        fnm=[fs.DAQ.BaseFileNameLong '*.tif']; 
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
fs.sAcquired=zeros(100000,1);
fs.sAvailable=zeros(100000,1);
fs.totaltime=zeros(100000,1);
fs.tictoc=zeros(100000,1);
fs.sAvailableFrames=zeros(100000,1);
fs.totaldelta=zeros(100000,1);
fs.phase=zeros(100000,1);


%% initilize some matrixes and start values
fs.DAQ.ch(1).CData1D=zeros(256,256);
fs.DAQ.ch(2).CData1D=zeros(256,256);

%fs.DAQ.ch(1).CData1D=zeros(240,256); %im size
%fs.DAQ.ch(2).CData1D=zeros(240,256);

fs.datasize=0;
fs.dataTail=[];
%fs.DAQ.InPointsToCollect - specified in SetupDAQ

%fs.timing.FastMirrorTriggerDelay=0; %just for test

fs.tifStackSize=150;

%% create cos LUT
cos_in=0:4000000;
fs.cos_table=uint16((1-cos(cos_in/4000))*255/2+1);
%fs.cos_table=uint16((1-cos(cos_in*128*2*pi/4000000))*255/2+1);


%% prepare RAM file for FastDisplay
fs.map = fopen(fs.DAQ.MemMapFile,'W');

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

%% open shutter
openShutter; pause(.125);

if fs.DAQ.DoNotSave
    fs.DAQ.modeString='Preview';
else
    fs.DAQ.modeString='Stream';
end

if fs.cycle.Do
    fs.DAQ.modeString='Stream-stack';
end

%% start DAQ
StartDAQ;
