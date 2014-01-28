function varargout = FastDisplay(varargin)
% FastDisplay M-file for FastDisplay.fig
%      FastDisplay, by itself, creates a new FastDisplay or raises the existing
%      singleton*.
%
%      H = FastDisplay returns the handle to a new FastDisplay or the handle to
%      the existing singleton*.
%
%      FastDisplay('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FastDisplay.M with the given input arguments.
%
%      FastDisplay('Property','Value',...) creates a new FastDisplay or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FastDisplay_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FastDisplay_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FastDisplay

% Last Modified by GUIDE v2.5 14-Sep-2007 11:34:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FastDisplay_OpeningFcn, ...
                   'gui_OutputFcn',  @FastDisplay_OutputFcn, ...
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


% --- Executes just before FastDisplay is made visible.
function FastDisplay_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FastDisplay (see VARARGIN)

global fd;

fd.iniDone=0;
% Choose default command line output for FastDisplay
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
fd.handles=handles;

OpenFastDisplay;

% --- Outputs from this function are returned to the command line.
function varargout = FastDisplay_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% Channels
% --- Executes on button press in chkChannel1.
function chkChannel1_Callback(hObject, eventdata, handles)
global fd;
if get(hObject,'Value')
    set(handles.txtMinChannel1,'Enable','on');
    set(handles.txtMaxChannel1,'Enable','on');
    set(handles.sliderMinChannel1,'Enable','on');
    set(handles.sliderMaxChannel1,'Enable','on');
    set(handles.btnAutoscaleChannel1,'Enable','on');
    if isempty(fd.img.Ch(1).h)
        fd.img.Ch(1).h=figure('NumberTitle','off','Name','Channels','WindowButtonDownFcn','SelectRectOnFig');
        fd.img.Ch(1).imH=image;
        fd.img.Ch(1).ax=gca;
        colormap('gray');
        set(fd.img.Ch(1).ax,'XLim',[1 256]);
        set(fd.img.Ch(1).ax,'YLim',[1 256]);
        set(fd.img.Ch(1).ax, 'YDir', 'reverse', 'XDir', 'reverse'); % match oculars
        set(fd.img.Ch(1).imH, 'EraseMode', 'none', 'CDataMapping','scaled');  
        set(fd.img.Ch(1).imH, 'CData', zeros(256,256,'uint16'));  
        truesize(fd.img.Ch(1).h,[256 256]);
        pos = get(fd.img.Ch(1).h,'position');
        pos(1)=600;pos(2)=600;
        set(fd.img.Ch(1).h,'position',pos);
%        set(fd.img.Ch(1).ax,'CLim',[fd.img.Ch(1).min+1 fd.img.Ch(1).max]);
    end
    fd.img.Ch(1).use=get(hObject,'Value');
else
    fd.img.Ch(1).use=get(hObject,'Value');
    set(handles.txtMinChannel1,'Enable','off');
    set(handles.txtMaxChannel1,'Enable','off');
    set(handles.sliderMinChannel1,'Enable','off');
    set(handles.sliderMaxChannel1,'Enable','off');
    set(handles.btnAutoscaleChannel1,'Enable','off');
%     if ~isempty(fd.img.Ch(1).h)
%         close(fd.img.Ch(1).h);
%         fd.img.Ch(1).h=[];
%     end
end

function txtMinChannel1_Callback(hObject, eventdata, handles)
global fd;

if ~isempty(str2num(get(handles.txtMinChannel1,'String')))
    fd.img.Ch(1).min=str2num(get(handles.txtMinChannel1,'String'));
end
if fd.img.Ch(1).min<get(handles.sliderMinChannel1,'Min')
    fd.img.Ch(1).min=get(handles.sliderMinChannel1,'Min');
end
if fd.img.Ch(1).min>get(handles.sliderMinChannel1,'Max')
    fd.img.Ch(1).min=get(handles.sliderMinChannel1,'Max');
end
set(handles.sliderMinChannel1,'Value',fd.img.Ch(1).min);
set(hObject,'String',fd.img.Ch(1).min);
% if ~isempty(fd.img.Ch(1).ax)
%     set(fd.img.Ch(1).ax,'CLim',[fd.img.Ch(1).min+1 fd.img.Ch(1).max]);
% end


% --- Executes during object creation, after setting all properties.
function txtMinChannel1_CreateFcn(hObject, eventdata, handles)
global fd;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderMinChannel1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMinChannel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fd;
fd.img.Ch(1).min=round(get(handles.sliderMinChannel1,'Value'));
if fd.img.Ch(1).max<fd.img.Ch(1).min
    fd.img.Ch(1).min=fd.img.Ch(1).max-1;
end
set(handles.txtMinChannel1,'String',num2str(fd.img.Ch(1).min));
% if ~isempty(fd.img.Ch(1).ax)
%     set(fd.img.Ch(1).ax,'CLim',[fd.img.Ch(1).min fd.img.Ch(1).max]);
% end
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderMinChannel1_CreateFcn(hObject, eventdata, handles)
global fd;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function txtMaxChannel1_Callback(hObject, eventdata, handles)
global fd;
if ~isempty(str2num(get(handles.txtMaxChannel1,'String')))
    fd.img.Ch(1).max=str2num(get(handles.txtMaxChannel1,'String'));
end
if fd.img.Ch(1).max<get(handles.sliderMaxChannel1,'Min')
    fd.img.Ch(1).max=get(handles.sliderMaxChannel1,'Min');
end
if fd.img.Ch(1).max>get(handles.sliderMaxChannel1,'Max')
    fd.img.Ch(1).max=get(handles.sliderMaxChannel1,'Max');
end

set(handles.sliderMaxChannel1,'Value',fd.img.Ch(1).max);
set(hObject,'String',fd.img.Ch(1).max);

% if ~isempty(fd.img.Ch(1).ax)
%     set(fd.img.Ch(1).ax,'CLim',[fd.img.Ch(1).min fd.img.Ch(1).max]);
% end
% Hints: get(hObject,'String') returns contents of txtMaxChannel1 as text
%        str2double(get(hObject,'String')) returns contents of txtMaxChannel1 as a double


% --- Executes during object creation, after setting all properties.
function txtMaxChannel1_CreateFcn(hObject, eventdata, handles)
global fd;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderMaxChannel1_Callback(hObject, eventdata, handles)
global fd;
fd.img.Ch(1).max=round(get(handles.sliderMaxChannel1,'Value'));
if fd.img.Ch(1).max<fd.img.Ch(1).min
    fd.img.Ch(1).max=fd.img.Ch(1).min+1;
end
set(handles.txtMaxChannel1,'String',num2str(fd.img.Ch(1).max));
% if ~isempty(fd.img.Ch(1).ax)
%     set(fd.img.Ch(1).ax,'CLim',[fd.img.Ch(1).min fd.img.Ch(1).max]);
% end



% --- Executes during object creation, after setting all properties.
function sliderMaxChannel1_CreateFcn(hObject, eventdata, handles)
global fd;
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
 

% --- Executes on button press in btnAutoscaleChannel1.
function btnAutoscaleChannel1_Callback(hObject, eventdata, handles)
global fd;
xoffset=round((256-fd.img.XSizePix)/2);
yoffset=round((256-fd.img.YSizePix)/2);
xend=min(xoffset+fd.img.XSizePix,256);
yend=min(yoffset+fd.img.YSizePix,256);
fd.img.Ch(1).CData = reshape(fd.img.Ch(1).CDataV2, 256,256);
set(handles.txtMinChannel1,'String',num2str(min(min(fd.img.Ch(1).CData(yoffset+1:yend,xoffset+1:xend)))));
set(handles.txtMaxChannel1,'String',num2str(max(max(fd.img.Ch(1).CData(yoffset+1:yend,xoffset+1:xend)))));
txtMinChannel1_Callback([],[],handles);
txtMaxChannel1_Callback([],[],handles);

% --- Executes on button press in chkChannel2.
function chkChannel2_Callback(hObject, eventdata, handles)
global fd;
if get(hObject,'Value')
    set(handles.txtMinChannel2,'Enable','on');
    set(handles.txtMaxChannel2,'Enable','on');
    set(handles.sliderMinChannel2,'Enable','on');
    set(handles.sliderMaxChannel2,'Enable','on');
    set(handles.btnAutoscaleChannel2,'Enable','on');
%     if isempty(fd.img.Ch(2).h)
%         fd.img.Ch(2).h=figure('NumberTitle','off','Name','Channel 2');
%         fd.img.Ch(2).imH=image;
%         fd.img.Ch(2).ax=gca;
%         colormap('gray');
%         set(fd.img.Ch(2).ax,'XLim',[1 256]);
%         set(fd.img.Ch(2).ax,'YLim',[1 256]);
%         set(fd.img.Ch(2).imH, 'EraseMode', 'none', 'CDataMapping','scaled');
%         set(fd.img.Ch(2).ax,'CLim',[fd.img.Ch(2).min fd.img.Ch(2).max]);
%     end
    fd.img.Ch(2).use=get(hObject,'Value');
else
    fd.img.Ch(2).use=get(hObject,'Value');
    set(handles.txtMinChannel2,'Enable','off');
    set(handles.txtMaxChannel2,'Enable','off');
    set(handles.sliderMinChannel2,'Enable','off');
    set(handles.sliderMaxChannel2,'Enable','off');
    set(handles.btnAutoscaleChannel2,'Enable','off');
    if ~isempty(fd.img.Ch(2).h)
        close(fd.img.Ch(2).h);
        fd.img.Ch(2).h=[];
    end
end


function txtMinChannel2_Callback(hObject, eventdata, handles)
global fd;
if ~isempty(str2num(get(handles.txtMinChannel2,'String')))
    fd.img.Ch(2).min=str2num(get(handles.txtMinChannel2,'String'));
end

if fd.img.Ch(2).min<get(handles.sliderMinChannel2,'Min')
    fd.img.Ch(2).min=get(handles.sliderMinChannel2,'Min');
end
if fd.img.Ch(2).min>get(handles.sliderMinChannel2,'Max')
    fd.img.Ch(2).min=get(handles.sliderMinChannel2,'Max');
end

set(handles.sliderMinChannel2,'Value',fd.img.Ch(2).min);
set(hObject,'String',fd.img.Ch(2).min);
% if ~isempty(fd.img.Ch(2).ax)
%     set(fd.img.Ch(2).ax,'CLim',[fd.img.Ch(2).min fd.img.Ch(2).max]);
% end

% --- Executes during object creation, after setting all properties.
function txtMinChannel2_CreateFcn(hObject, eventdata, handles)
global fd;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMaxChannel2_Callback(hObject, eventdata, handles)
global fd;
if ~isempty(str2num(get(handles.txtMaxChannel2,'String')))
    fd.img.Ch(2).max=str2num(get(handles.txtMaxChannel2,'String'));
end
if fd.img.Ch(2).max<get(handles.sliderMaxChannel2,'Min')
    fd.img.Ch(2).max=get(handles.sliderMaxChannel2,'Min');
end
if fd.img.Ch(2).max>get(handles.sliderMaxChannel2,'Max')
    fd.img.Ch(2).max=get(handles.sliderMaxChannel2,'Max');
end
set(handles.sliderMaxChannel2,'Value',fd.img.Ch(2).max);
set(hObject,'String',fd.img.Ch(2).max);
% if ~isempty(fd.img.Ch(2).ax)
%     set(fd.img.Ch(2).ax,'CLim',[fd.img.Ch(2).min fd.img.Ch(2).max]);
% end


% --- Executes during object creation, after setting all properties.
function txtMaxChannel2_CreateFcn(hObject, eventdata, handles)
global fd;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderMinChannel2_Callback(hObject, eventdata, handles)
global fd;
fd.img.Ch(2).min=round(get(handles.sliderMinChannel2,'Value'));
if fd.img.Ch(2).max<fd.img.Ch(2).min
    fd.img.Ch(2).min=fd.img.Ch(2).max-1;
end
set(handles.txtMinChannel2,'String',num2str(fd.img.Ch(2).min));
% if ~isempty(fd.img.Ch(2).ax)
%     set(fd.img.Ch(2).ax,'CLim',[fd.img.Ch(2).min fd.img.Ch(2).max]);
% end


% --- Executes during object creation, after setting all properties.
function sliderMinChannel2_CreateFcn(hObject, eventdata, handles)
global fd;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderMaxChannel2_Callback(hObject, eventdata, handles)
global fd;
fd.img.Ch(2).max=round(get(handles.sliderMaxChannel2,'Value'));
if fd.img.Ch(2).max<fd.img.Ch(2).min
    fd.img.Ch(2).max=fd.img.Ch(2).min+1;
end
set(handles.txtMaxChannel2,'String',num2str(fd.img.Ch(2).max));
% if ~isempty(fd.img.Ch(2).ax)
%     set(fd.img.Ch(2).ax,'CLim',[fd.img.Ch(2).min fd.img.Ch(2).max]);
% end


% --- Executes during object creation, after setting all properties.
function sliderMaxChannel2_CreateFcn(hObject, eventdata, handles)
global fd;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnAutoscaleChannel2.
function btnAutoscaleChannel2_Callback(hObject, eventdata, handles)
global fd;
xoffset=round((256-fd.img.XSizePix)/2);
yoffset=round((256-fd.img.YSizePix)/2);
xend=min(xoffset+fd.img.XSizePix,256);
yend=min(yoffset+fd.img.YSizePix,256);
fd.img.Ch(2).CData = reshape(fd.img.Ch(2).CDataV2, 256,256);
set(handles.txtMinChannel2,'String',num2str(min(min(fd.img.Ch(2).CData(yoffset+1:yend,xoffset+1:xend)))));
set(handles.txtMaxChannel2,'String',num2str(max(max(fd.img.Ch(2).CData(yoffset+1:yend,xoffset+1:xend)))));
txtMinChannel2_Callback([],[],handles);
txtMaxChannel2_Callback([],[],handles);

function txtXSize_Callback(hObject, eventdata, handles)
global fd;
if ~isempty(str2num(get(hObject,'String')))
    fd.img.XSizePix=str2num(get(hObject,'String'));
    fd.img.XSizePix=min(fd.img.XSizePix,256);
end
set(hObject,'String',fd.img.XSizePix);
% Hints: get(hObject,'String') returns contents of txtXSize as text
%        str2double(get(hObject,'String')) returns contents of txtXSize as a double


% --- Executes during object creation, after setting all properties.
function txtXSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txtYSize_Callback(hObject, eventdata, handles)
global fd;
if ~isempty(str2num(get(hObject,'String')))
    fd.img.YSizePix=str2num(get(hObject,'String'));
    fd.img.YSizePix=min(fd.img.YSizePix,256);
end
set(hObject,'String',fd.img.YSizePix);
% Hints: get(hObject,'String') returns contents of txtYSize as text
%        str2double(get(hObject,'String')) returns contents of txtYSize as a double

% --- Executes during object creation, after setting all properties.
function txtYSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
global fd;
try
    if ~isempty(fd.img.Ch(2).h)
        close(fd.img.Ch(2).h);
        fd.img.Ch(2).h=[];
    end
    if ~isempty(fd.img.Ch(1).h)
        close(fd.img.Ch(1).h);
        fd.img.Ch(1).h=[];
    end


% Hint: delete(hObject) closes the figure
delete(hObject);
close all;
clear all
button = questdlg('Do you want to close Matlab?','Exit','Yes','No','No');
if strcmp(button,'Yes')
    exit
end
catch
end

function txtFastMirrorTriggerDelay_Callback(hObject, eventdata, handles)
global fd;
r=fd.Display.DoDisplay;
btnStopDisplay_Callback([],[], handles);

if ~isempty(str2num(get(hObject,'String')))
    fd.timing.FastMirrorTriggerDelay=str2num(get(hObject,'String'));
end
set(hObject,'String',fd.timing.FastMirrorTriggerDelay);

if r==1
    btnStartDisplay_Callback([],[], handles);
end


% --- Executes during object creation, after setting all properties.
function txtFastMirrorTriggerDelay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtAvgOnDisplay_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of txtAvgOnDisplay as text
%        str2double(get(hObject,'String')) returns contents of txtAvgOnDisplay as a double
global fd;

r=fd.Display.DoDisplay;
btnStopDisplay_Callback([],[],handles);
if ~isempty(str2num(get(hObject,'String')))
    fd.img.AvgOnDisplay=str2num(get(hObject,'String'));
end
if r==1
    btnStartDisplay_Callback([],[],handles);
end



% --- Executes during object creation, after setting all properties.
function txtAvgOnDisplay_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in btnStartDisplay.
function btnStartDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to btnStartDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fd;
set(handles.lblStatus,'String','Running Display');
set(handles.lblStatus,'ForegroundColor',[1 0 0]);
fd.Display.DoDisplay=1;
set(handles.btnStartDisplay,'Enable','off');
DoDisplay

% --- Executes on button press in btnStopDisplay.
function btnStopDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to btnStopDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fd;
set(handles.lblStatus,'String','Display Stopped');
set(handles.lblStatus,'ForegroundColor',[0 1 0]);

if ~isempty(fd.DisplayTimer)

    stop(fd.DisplayTimer);
 %   pause(0.5);
    delete(fd.DisplayTimer); 
    fd.DisplayTimer = [];
end
fd.Display.DoDisplay=0;
set(handles.btnStartDisplay,'Enable','on');


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GenericKeyPress_FastDisplay



% --- Executes on button press in btnCalculateDelay.
function btnCalculateDelay_Callback(hObject, eventdata, handles)
% hObject    handle to btnCalculateDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fd;
CalculateDelay;

% --- Executes on button press in zoomToggle.
function zoomToggle_Callback(hObject, eventdata, handles)
% hObject    handle to zoomToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoomToggle

global fd;

if get(hObject,'Value');
    truesize(fd.img.Ch(1).h, [512,512]);
else
    truesize(fd.img.Ch(1).h, [256,256]);
end


% --- Executes on button press in statsToggleButton.
function statsToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to statsToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of statsToggleButton

global fd;

fd.stats.show = get(hObject,'Value');
r=fd.Display.DoDisplay;
btnStopDisplay_Callback([],[], handles);

if fd.stats.show
    fd.stats.nsamples = 1000;
    fd.stats.h=figure('NumberTitle','off','Name','Stats');
    fd.stats.Ch1 = nan*zeros(fd.stats.nsamples,3);
    fd.stats.Ch2 = nan*zeros(fd.stats.nsamples,3);
    fd.stats.plot1 = plot(fd.stats.Ch1,'g');
    hold on;
    fd.stats.plot2 = plot(fd.stats.Ch2,'r');    
    title('Stats');
    fd.stats.ax = gca;
    set(fd.stats.plot1(1),'linestyle','-');
    set(fd.stats.plot1(2),'linewidth',2);
    set(fd.stats.plot1(3),'linestyle','-');
    set(fd.stats.plot2(1),'linestyle','-');
    set(fd.stats.plot2(2),'linewidth',2);
    set(fd.stats.plot2(3),'linestyle','-');
    set(fd.stats.ax,'ylimmode','manual','ylim', [1900 4096]);
    set(fd.stats.ax,'xlimmode','manual','xlim', [0 1000]);
    set(hObject,'String','Stats - close');
else
    delete(fd.stats.h);
    set(hObject,'String','Stats - show');
end

if r==1
    btnStartDisplay_Callback([],[], handles);
end