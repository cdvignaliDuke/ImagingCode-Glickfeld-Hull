function varargout = FastScanner(varargin)
% FASTSCANNER M-file for FastScanner.fig
%      FASTSCANNER, by itself, creates a new FASTSCANNER or raises the existing
%      singleton*.
%
%      H = FASTSCANNER returns the handle to a new FASTSCANNER or the handle to
%      the existing singleton*.
%
%      FASTSCANNER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FASTSCANNER.M with the given input arguments.
%
%      FASTSCANNER('Property','Value',...) creates a new FASTSCANNER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FastScanner_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FastScanner_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FastScanner

% Last Modified by GUIDE v2.5 12-Jun-2009 15:08:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FastScanner_OpeningFcn, ...
                   'gui_OutputFcn',  @FastScanner_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FastScanner is made visible.
function FastScanner_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FastScanner (see VARARGIN)

global fs;

fs.iniDone=0;
% Choose default command line output for FastScanner
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
fs.handles=handles;
set(handles.btnStopAquisition,'Background',.9*[1 1 1]);
OpenFastScanner;

% --- Outputs from this function are returned to the command line.
function varargout = FastScanner_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliderXPos_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global fs;
fs.position.x=get(hObject,'Value');
set(handles.txtXPos,'String',num2str(fs.position.x));


% --- Executes during object creation, after setting all properties.
function sliderXPos_CreateFcn(hObject, eventdata, handles)
global fs;
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function txtXPos_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of txtXPos as text
%        str2double(get(hObject,'String')) returns contents of txtXPos as a double
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.position.x=str2num(get(hObject,'String'));
end
if fs.position.x<get(handles.sliderXPos,'Min')
    fs.position.x=get(handles.sliderXPos,'Min');
end
if fs.position.x>get(handles.sliderXPos,'Max')
    fs.position.x=get(handles.sliderXPos,'Max');
end
set(handles.sliderXPos,'Value',fs.position.x);
set(hObject,'String',fs.position.x);

% --- Executes during object creation, after setting all properties.
function txtXPos_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliderYPos_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global fs;
fs.position.y=get(hObject,'Value');
set(handles.txtYPos,'String',num2str(fs.position.y));

% --- Executes during object creation, after setting all properties.
function sliderYPos_CreateFcn(hObject, eventdata, handles)
global fs;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function txtYPos_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of txtYPos as text
%        str2double(get(hObject,'String')) returns contents of txtYPos as a double
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.position.y=str2num(get(hObject,'String'));
end
if fs.position.y<get(handles.sliderYPos,'Min')
    fs.position.y=get(handles.sliderYPos,'Min');
end
if fs.position.y>get(handles.sliderYPos,'Max')
    fs.position.y=get(handles.sliderYPos,'Max');
end
set(handles.sliderYPos,'Value',fs.position.y);
set(hObject,'String',fs.position.y);

% --- Executes during object creation, after setting all properties.
function txtYPos_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderZPos_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global fs;
fs.position.z=round(get(hObject,'Value'));
set(handles.txtZPos,'String',num2str(fs.position.z));
set(handles.lblCurrentZ,'String',num2str(fs.position.z));
SetPositionZ_no_wait(fs.position.z);
fs.cycle.ZStart=fs.position.z;
set(fs.handles.txtStartZ,'String',fs.cycle.ZStart);


% --- Executes during object creation, after setting all properties.
function sliderZPos_CreateFcn(hObject, eventdata, handles)
global fs;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function txtZPos_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of txtZPos as text
%        str2double(get(hObject,'String')) returns contents of txtZPos as a double
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.position.z=str2num(get(hObject,'String'));
    if fs.position.z<get(handles.sliderZPos,'Min')
        fs.position.z=get(handles.sliderZPos,'Min');
    end
    if fs.position.z>get(handles.sliderZPos,'Max')
        fs.position.z=get(handles.sliderZPos,'Max');
    end
    set(handles.sliderZPos,'Value',fs.position.z);
    set(hObject,'String',fs.position.z);

    SetPositionZ(fs.position.z);
    fs.cycle.ZStart=fs.position.z;
    set(fs.handles.txtStartZ,'String',fs.cycle.ZStart);

end

% --- Executes during object creation, after setting all properties.
function txtZPos_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderZoom_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global fs;
fs.DAQ.adjusting=1;
fs.position.zoom=get(handles.sliderZoom,'Value');
set(handles.txtZoom,'String',num2str(fs.position.zoom));
%set voltage, that controls X mirror amplitude here
if fs.iniDone
%    putsample(fs.DAQ.aoXYMirrors,[fs.position.scale/fs.position.zoom fs.position.scale*fs.position.XtoYscale/fs.position.zoom]); %here we assume that max voltage is 5V
    putsample(fs.DAQ.aoXYMirrors,[fs.position.zoom fs.position.XtoYscale*fs.position.zoom*fs.position.YXratio]); %here we assume that max voltage is 5V
end

ranges = [.2 .5 1 2 5];

if fs.DAQ.created==1
    % driving signal is 0 to 5 V. 
    % X and Y signals always slightly below driving signal
    temp = ranges - fs.position.zoom;     
    bestrange = min(ranges(find(temp>=0)));

    r=fs.DAQ.ai.Running;

    if strcmpi(r,'On')
        btnStopAquisition_Callback([],[], handles);
    while strcmpi(fs.DAQ.ai.Running,'On')
        pause(0.1);
    end
    end

    fs.DAQ.ai.Channel(3:4).InputRange=[-bestrange bestrange];
    
% 	disp([-bestrange bestrange]);
    
    if strcmpi(r,'On')
    if fs.DAQ.DoNotSave
        btnScan_Callback([],[], handles);
    else
        btnStream_Callback([],[], handles);
    end
    end
end
fs.DAQ.adjusting=0;



% --- Executes during object creation, after setting all properties.
function sliderZoom_CreateFcn(hObject, eventdata, handles)
global fs;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function txtZoom_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of txtZoom as text
%        str2double(get(hObject,'String')) returns contents of txtZoom as a double
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.position.zoom=str2num(get(hObject,'String'));
end
if fs.position.zoom<get(handles.sliderZoom,'Min')
    fs.position.zoom=get(handles.sliderZoom,'Min');
end
if fs.position.zoom>get(handles.sliderZoom,'Max')
    fs.position.zoom=get(handles.sliderZoom,'Max');
end
set(handles.sliderZoom,'Value',fs.position.zoom);
set(hObject,'String',fs.position.zoom);
sliderZoom_Callback([],[],handles);

% --- Executes during object creation, after setting all properties.
function txtZoom_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in btnStopAquisition.
function btnStopAquisition_Callback(hObject, eventdata, handles)
global fs;

fs.DAQ.Running=0;

if fs.DAQ.started && ~fs.DAQ.Streaming
    StopAcquisition;
end

set(handles.btnStream,'Enable','on');
set(handles.btnScan,'Enable','on');
set(handles.btnStartCycles,'Enable','on');


% --- Executes on button press in btnStream.
function btnStream_Callback(hObject, eventdata, handles)
global fs;
fs.cycle.Do=0;
fs.DAQ.DoNotSave=0;
fs.DAQ.Running=1;
[p msg]=RunAcquisition;
if p
set(handles.lblProgress,'String','Streaming to disk');
set(handles.btnStopAquisition,'BackgroundColor',get(handles.uipanel1,'BackgroundColor'));
set(handles.btnStream,'BackgroundColor',.9*[1 1 1]);
set(handles.btnScan,'BackgroundColor',get(handles.uipanel1,'BackgroundColor'));
set(handles.btnStream,'Enable','off');
set(handles.btnScan,'Enable','off');
set(handles.btnStartCycles,'Enable','off');
else
set(handles.lblProgress,'String',msg);    
end


function txtStartX_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.XStart=str2num(get(hObject,'String'));
end
set(hObject,'String',fs.cycle.XStart);

% --- Executes during object creation, after setting all properties.
function txtStartX_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txtStartY_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.YStart=str2num(get(hObject,'String'));
end
set(hObject,'String',fs.cycle.YStart);


% --- Executes during object creation, after setting all properties.
function txtStartY_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txtStartZ_Callback(hObject, eventdata, handles)
%if ~isempty(str2num(get(hObject,'String')))
%    fs.cycle.ZStart=str2num(get(hObject,'String'));
%end
%set(hObject,'String',fs.cycle.ZStart);

% --- Executes during object creation, after setting all properties.
function txtStartZ_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtStepX_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.XStep=str2num(get(hObject,'String'));
end
set(hObject,'String',fs.cycle.XStep);


% --- Executes during object creation, after setting all properties.
function txtStepX_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtStepY_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.YStep=str2num(get(hObject,'String'));
end
set(hObject,'String',fs.cycle.YStep);


% --- Executes during object creation, after setting all properties.
function txtStepY_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txtStepZ_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.ZStep=str2num(get(hObject,'String'));
end
set(hObject,'String',fs.cycle.ZStep);


% --- Executes during object creation, after setting all properties.
function txtStepZ_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtNofSteps_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.NSteps=str2num(get(hObject,'String'));
end
set(hObject,'String',fs.cycle.NSteps);


% --- Executes during object creation, after setting all properties.
function txtNofSteps_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtNofRepetitions_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.NReps=str2num(get(hObject,'String'));
end
set(hObject,'String',fs.cycle.NReps);


% --- Executes during object creation, after setting all properties.
function txtNofRepetitions_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkAvg.
function chkAvg_Callback(hObject, eventdata, handles)
global fs;
fs.cycle.Avg=get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of chkAvg



function txtFramesPerStep_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.NFramesPrStep=str2num(get(hObject,'String'));
end
set(hObject,'String',fs.cycle.NFramesPrStep);
% Hints: get(hObject,'String') returns contents of txtFramesPerStep as text
%        str2double(get(hObject,'String')) returns contents of txtFramesPerStep as a double


% --- Executes during object creation, after setting all properties.
function txtFramesPerStep_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnStartCycles.
function btnStartCycles_Callback(hObject, eventdata, handles)
global fs;
%if ~isempty(str2num(get(hObject,'String')))
    fs.cycle.ZStart=ReadZPos;
%end
set(fs.handles.txtStartZ,'String',fs.cycle.ZStart);
fs.cycle.Do=1;
fs.DAQ.DoNotSave=0;
fs.DAQ.Running=1;

if fs.stage.JS==1
    btnJS_Callback([],[], handles);
end

[p msg]=RunAcquisition;
if p
set(handles.lblProgress,'String','Streaming to disk - stack');
set(handles.btnStopAquisition,'BackgroundColor',get(handles.uipanel1,'BackgroundColor'));
set(handles.btnStream,'BackgroundColor',.9*[1 1 1]);
set(handles.btnScan,'BackgroundColor',get(handles.uipanel1,'BackgroundColor'));

set(handles.btnStream,'Enable','off');
set(handles.btnScan,'Enable','off');
set(handles.btnStartCycles,'Enable','off');
else
set(handles.lblProgress,'String',msg);    
end

% --- Executes on button press in btnStopCycles.
function btnStopCycles_Callback(hObject, eventdata, handles)
global fs;
btnStopAquisition_Callback([],[], handles);
set(handles.btnStream,'Enable','on');
set(handles.btnScan,'Enable','on');
set(handles.btnStartCycles,'Enable','on');


function txtPCellMax_Callback(hObject, eventdata, handles)
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.pCell.max=str2num(get(hObject,'String'));
end
if fs.pCell.max<get(handles.sliderPCellMax,'Min')
    fs.pCell.max=get(handles.sliderPCellMax,'Min');
end
if fs.pCell.max>get(handles.sliderPCellMax,'Max')
    fs.pCell.max=get(handles.sliderPCellMax,'Max');
end
set(handles.sliderPCellMax,'Value',fs.pCell.max);
set(hObject,'String',fs.pCell.max);
SetPCell;
% Hints: get(hObject,'String') returns contents of txtPCellMax as text
%        str2double(get(hObject,'String')) returns contents of txtPCellMax as a double


% --- Executes during object creation, after setting all properties.
function txtPCellMax_CreateFcn(hObject, eventdata, handles)
global fs;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function sliderPCellMax_Callback(hObject, eventdata, handles)
global fs;
fs.pCell.max=get(hObject,'Value');
fs.pCell.max=round(fs.pCell.max);
if fs.pCell.max<get(handles.sliderPCellMax,'Min')
    fs.pCell.max=get(handles.sliderPCellMax,'Min');
end
if fs.pCell.max>get(handles.sliderPCellMax,'Max')
    fs.pCell.max=get(handles.sliderPCellMax,'Max');
end
set(handles.sliderPCellMax,'Value',fs.pCell.max);

set(handles.txtPCellMax,'String',num2str(fs.pCell.max));
SetPCell;

% --- Executes during object creation, after setting all properties.
function sliderPCellMax_CreateFcn(hObject, eventdata, handles)
global fs;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnIniStage.
function btnIniStage_Callback(hObject, eventdata, handles)
global fs;
if IniStage
    set(handles.txtZPos,'Enable','on');
    set(handles.sliderZPos,'Enable','on');
%    set(handles.txtStartZ,'Enable','on');
    set(handles.txtStepZ,'Enable','on');
end


% --- Executes on button press in btnIniScanner.
function btnIniScanner_Callback(hObject, eventdata, handles)
global fs;
GetScannerParameters

function txtXSize_Callback(hObject, eventdata, handles)
% fs.img.XSizePix is not in use in FastScanner
% global fs;
% if ~isempty(str2num(get(hObject,'String')))
%     fs.img.XSizePix=str2num(get(hObject,'String'));
% end
% set(hObject,'String',fs.img.XSizePix);

% Hints: get(hObject,'String') returns contents of txtXSize as text
%        str2double(get(hObject,'String')) returns contents of txtXSize as a double


% --- Executes during object creation, after setting all properties.
function txtXSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtYSize_Callback(hObject, eventdata, handles)
% fs.img.YSizePix is not in use in FastScanner

% global fs;
% if ~isempty(str2num(get(hObject,'String')))
%     fs.img.YSizePix=str2num(get(hObject,'String'));
% end
% set(hObject,'String',fs.img.YSizePix);
% Hints: get(hObject,'String') returns contents of txtYSize as text
%        str2double(get(hObject,'String')) returns contents of txtYSize as a double


% --- Executes during object creation, after setting all properties.
function txtYSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
delete(hObject);
CloseFastScanner
return;




% --- Executes on selection change in lstRanges.
function lstRanges_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns lstRanges contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstRanges
global fs;
fs.DAQ.adjusting=1;
r=fs.DAQ.ai.Running;

if strcmpi(r,'On')
    btnStopAquisition_Callback([],[], handles);
    while strcmpi(fs.DAQ.ai.Running,'On')
        pause(0.1);
    end
end
index_selected = get(hObject,'Value');
fs.DAQ.ai.Channel(1:2).InputRange=fs.DAQ.InputRanges(index_selected,:);

if strcmpi(r,'On')
    if fs.DAQ.DoNotSave
        btnScan_Callback([],[], handles);
    else
        btnStream_Callback([],[], handles);
    end
end
fs.DAQ.adjusting=0;

% --- Executes during object creation, after setting all properties.
function lstRanges_CreateFcn(hObject, eventdata, handles)
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtBaseFileName_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of txtBaseFileName as text
%        str2double(get(hObject,'String')) returns contents of txtBaseFileName as a double
global fs;
fs.DAQ.BaseFileName=get(hObject,'String');

% --- Executes during object creation, after setting all properties.
function txtBaseFileName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function lblCurrentZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblCurrentZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on button press in chkTestMode.
function chkTestMode_Callback(hObject, eventdata, handles)
% hObject    handle to chkTestMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fs;
fs.DAQ.TestMode=get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of chkTestMode


% --- Executes on button press in btnScan.
function btnScan_Callback(hObject, eventdata, handles)
% hObject    handle to btnScan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global fs;
fs.cycle.Do=0;
fs.DAQ.DoNotSave=1;
fs.DAQ.Running=1;
RunAcquisition;

set(handles.lblProgress,'String','Previewing');
set(handles.btnStopAquisition,'BackgroundColor',get(handles.uipanel1,'BackgroundColor'));
set(handles.btnStream,'BackgroundColor',get(handles.uipanel1,'BackgroundColor'));
set(handles.btnScan,'BackgroundColor',.9*[1 1 1]);
set(handles.btnStream,'Enable','off');
set(handles.btnScan,'Enable','off');
set(handles.btnStartCycles,'Enable','off');



function txtPCellMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtPCellMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPCellMin as text
%        str2double(get(hObject,'String')) returns contents of txtPCellMin as a double
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.pCell.min=str2num(get(hObject,'String'));
end
if fs.pCell.min<get(handles.sliderPCellMin,'Min')
    fs.pCell.min=get(handles.sliderPCellMin,'Min');
end
if fs.pCell.min>get(handles.sliderPCellMin,'Max')
    fs.pCell.min=get(handles.sliderPCellMin,'Max');
end
set(handles.sliderPCellMin,'Value',fs.pCell.min);
set(hObject,'String',fs.pCell.min);

% --- Executes during object creation, after setting all properties.
function txtPCellMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPCellMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderPCellMin_Callback(hObject, eventdata, handles)
% hObject    handle to sliderPCellMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global fs;
fs.pCell.min=get(hObject,'Value');
fs.pCell.min=round(fs.pCell.min);
if fs.pCell.min<get(handles.sliderPCellMin,'Min')
    fs.pCell.min=get(handles.sliderPCellMin,'Min');
end
if fs.pCell.min>get(handles.sliderPCellMin,'Max')
    fs.pCell.min=get(handles.sliderPCellMin,'Max');
end
set(handles.sliderPCellMin,'Value',fs.pCell.min);

set(handles.txtPCellMin,'String',num2str(fs.pCell.min));


% --- Executes during object creation, after setting all properties.
function sliderPCellMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderPCellMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in chkBlockEnds.
function chkBlockEnds_Callback(hObject, eventdata, handles)
global fs;
fs.DAQ.BlockEnds=get(hObject,'Value');




% --- Executes on key press with focus on figure1 and no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function GenericKeyPress
%just placeholder needed for correct GUI operation



% --- Executes on selection change in lstboxCalibFiles.
function lstboxCalibFiles_Callback(hObject, eventdata, handles)
% hObject    handle to lstboxCalibFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lstboxCalibFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstboxCalibFiles

global fs;

contents = get(hObject,'String');
fs.pCell.LUTfile=fullfile(fs.pCell.LUTdir,contents{get(hObject,'Value')});
%fs.pCell.LUTfile=contents{get(hObject,'Value')};
fs.pCell.LUT=ReadPockelsVoltageLUT(fs.pCell.LUTfile);
UpdateFastScannerByParam;

SetPCell;



% --- Executes during object creation, after setting all properties.
function lstboxCalibFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstboxCalibFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkMainShutter.
function chkMainShutter_Callback(hObject, eventdata, handles)
% hObject    handle to chkMainShutter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkMainShutter

if get(hObject,'Value')
    openMainShutter;
else
    closeMainShutter;
end



function txtYXratio_Callback(hObject, eventdata, handles)
% hObject    handle to txtYXratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtYXratio as text
%        str2double(get(hObject,'String')) returns contents of txtYXratio as a double
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.position.YXratio=str2num(get(hObject,'String'));
end
if fs.iniDone
%    putsample(fs.DAQ.aoXYMirrors,[fs.position.scale/fs.position.zoom fs.position.scale*fs.position.XtoYscale/fs.position.zoom]); %here we assume that max voltage is 5V
    putsample(fs.DAQ.aoXYMirrors,[fs.position.zoom fs.position.XtoYscale*fs.position.zoom*fs.position.YXratio]); %here we assume that max voltage is 5V
end

% --- Executes during object creation, after setting all properties.
function txtYXratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtYXratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btnJS.
function btnJS_Callback(hObject, eventdata, handles)
% hObject    handle to btnJS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fs;
if fs.stage.ready==1
if fs.stage.JS==0
    %transfer control to JS
    fwrite(fs.stage.comport, 'EX JOYSTICK ON');
    fwrite(fs.stage.comport, [13]);
    set(handles.btnJS,'String','JS is ON');
    set(handles.btnJS,'ForegroundColor',[1 0 0]);
    set(handles.txtZPos,'Enable','off');
    set(handles.sliderZPos,'Enable','off');
    set(handles.txtStepZ,'Enable','off');
    fs.stage.JS=1;
else
    %disable JS
    fwrite(fs.stage.comport, 'EX JOYSTICK OFF');
    fwrite(fs.stage.comport, [13]);
    set(handles.btnJS,'String','JS is OFF');
    set(handles.btnJS,'ForegroundColor',[0 0 0]);
    set(handles.txtZPos,'Enable','on');
    set(handles.sliderZPos,'Enable','on');
    set(handles.txtStepZ,'Enable','on');
    fs.position.z=ReadZPos;
    if fs.position.z<get(handles.sliderZPos,'Min')
        fs.position.z=get(handles.sliderZPos,'Min');
    end
    if fs.position.z>get(handles.sliderZPos,'Max')
        fs.position.z=get(handles.sliderZPos,'Max');
    end
    set(handles.sliderZPos,'Value',fs.position.z);
    set(handles.txtZPos,'String',num2str(fs.position.z));
    set(handles.lblCurrentZ,'String',num2str(fs.position.z));
    fs.cycle.ZStart=fs.position.z;
    set(fs.handles.txtStartZ,'String',fs.cycle.ZStart);
    fs.stage.JS=0;
end
end


function txtDelay_Callback(hObject, eventdata, handles)
% hObject    handle to txtDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtDelay as text
%        str2double(get(hObject,'String')) returns contents of txtDelay as a double
global fs;
if ~isempty(str2num(get(hObject,'String')))
    fs.DAQ.delay=str2num(get(hObject,'String'));
end


% --- Executes during object creation, after setting all properties.
function txtDelay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btnCalculateDelay.
function btnCalculateDelay_Callback(hObject, eventdata, handles)
% hObject    handle to btnCalculateDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fs;
fs.DAQ.calculating_delay=1;
