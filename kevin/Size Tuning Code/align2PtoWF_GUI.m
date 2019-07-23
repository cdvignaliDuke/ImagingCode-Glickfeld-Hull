function varargout = align2PtoWF_GUI(varargin)
% ALIGN2PTOWF_GUI MATLAB code for align2PtoWF_GUI.fig
%      ALIGN2PTOWF_GUI, by itself, creates a new ALIGN2PTOWF_GUI or raises the existing
%      singleton*.
%
%      H = ALIGN2PTOWF_GUI returns the handle to a new ALIGN2PTOWF_GUI or the handle to
%      the existing singleton*.
%
%      ALIGN2PTOWF_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALIGN2PTOWF_GUI.M with the given input arguments.
%
%      ALIGN2PTOWF_GUI('Property','Value',...) creates a new ALIGN2PTOWF_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before align2PtoWF_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to align2PtoWF_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help align2PtoWF_GUI

% Last Modified by GUIDE v2.5 03-Jul-2019 17:02:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @align2PtoWF_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @align2PtoWF_GUI_OutputFcn, ...
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


% --- Executes just before align2PtoWF_GUI is made visible.
function align2PtoWF_GUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for align2PtoWF_GUI
handles.output = hObject;

% set defaults
handles.fixedim = varargin{1};
handles.movingim = varargin{2};

handles.Xpos = 120;
handles.Ypos = 180;
handles.rotation = 0;
handles.scale = 20;
handles.fixedim_resize = imresize(handles.fixedim,handles.scale);

set(handles.slider1,'Min',0)
set(handles.slider1,'Value',handles.Xpos)
set(handles.slider1,'Max',round(size(handles.fixedim,2)))
set(handles.slider2,'Min',0)
set(handles.slider2,'Value',handles.Ypos)
set(handles.slider2,'Max',round(size(handles.fixedim,1)))
set(handles.slider3,'Min',-360)
set(handles.slider3,'Value',handles.rotation)
set(handles.slider3,'Max',360)
set(handles.slider4,'Min',10)
set(handles.slider4,'Value',handles.scale)
set(handles.slider4,'Max',30)

guidata(hObject,handles)

update_axes(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);


function update_axes(hObject, eventdata, handles)
set(handles.text8,'String','Calculating...');
% guidata(hObject,handles)
moving_overlay = 0*handles.fixedim_resize;
moving = handles.movingim';
% if handles.scale==1
%     moving_scaled = handles.movingim;
% else
%     moving_scaled = imresize(handles.movingim,handles.scale);
% end
% if handles.rotation == 0
%     moving_rot = moving_scaled;
% else
%     moving_rot = imrotate(moving_scaled,handles.rotation);
% end
% [r, c] = ind2sub(size(moving_rot),1:numel(moving_rot));
% rshift = 10*handles.Ypos-size(moving_rot,1)/2;
% cshift = 10*handles.Xpos-size(moving_rot,2)/2;
% r_new = r+rshift; c_new = c+cshift;
% keep = ((r+rshift)<=size(moving_overlay,1)) & ((c+cshift)<=size(moving_overlay,2)) & ((r+rshift)>0) & ((c+cshift)>0);
% r_new = r_new(keep); c_new=c_new(keep);
% ind = sub2ind(size(moving_overlay),r_new,c_new);
% movingoverlay(ind) = handles.movingim(keep);
xOff = handles.Ypos*handles.scale-size(moving,2)/2;
yOff = handles.Xpos*handles.scale-size(moving,1)/2;
[xx, yy] = ndgrid(1:size(moving_overlay,1),1:size(moving_overlay,2));
mapX = round((xx-xOff)*cosd(handles.rotation) + (yy-yOff)*sind(handles.rotation) +size(moving,2)/2);
mapY = round((yy-yOff)*cosd(handles.rotation) - (xx-xOff)*sind(handles.rotation) +size(moving,1)/2);
zeromask = mapX<1 | mapX>size(moving,2) | mapY<1 | mapY>size(moving,1);
mapX(zeromask) = 0; mapY(zeromask)=0;
useind = find(~zeromask);
mapInd = sub2ind(size(moving),mapY(useind),mapX(useind));
moving_overlay(useind) = moving(mapInd);
moving_overlay = 255*moving_overlay/max(moving_overlay(:));
set(handles.text8,'String','Plotting...');
%imshowpair(handles.fixedim_resize,moving_overlay,'blend','Scaling','joint','Parent',handles.axes1)
meanIm = (handles.fixedim_resize + moving_overlay)/2;
imagesc(meanIm,'Parent',handles.axes1);
set(handles.axes1,'XTickLabel',{round(get(handles.axes1,'XTick')/handles.scale)})
set(handles.axes1,'YTickLabel',{round(get(handles.axes1,'YTick')/handles.scale)})
colormap gray
% hold on
% plot(yOff,xOff,'rx')
% hold off
set(handles.text8,'String','Done plotting');
guidata(hObject,handles)

% --- Outputs from this function are returned to the command line.
function varargout = align2PtoWF_GUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
handles.Xpos = round(get(handles.slider1,'Value'));
set(handles.edit1,'String',num2str(handles.Xpos));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
handles.Ypos = round(get(handles.slider2,'Value'));
set(handles.edit2,'String',num2str(handles.Ypos));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
handles.rotation = get(handles.slider3,'Value');
set(handles.edit3,'String',num2str(handles.rotation));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
handles.scale = get(handles.slider4,'Value');
set(handles.edit4,'String',num2str(handles.scale));
handles.fixedim_resize = imresize(handles.fixedim,handles.scale);
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
handles.Xpos = round(str2num(get(handles.edit1,'String')));
handles.Xpos = max(min(handles.Xpos,get(handles.slider1,'Max')),get(handles.slider1,'Min'));
set(handles.slider1,'Value',handles.scale);
set(handles.edit1,'String',num2str(handles.Xpos));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
handles.Ypos = round(str2num(get(handles.edit2,'String')));
handles.Ypos = max(min(handles.Ypos,get(handles.slider2,'Max')),get(handles.slider2,'Min'));
set(handles.slider2,'Value',handles.Ypos);
set(handles.edit2,'String',num2str(handles.Ypos));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
handles.rotation = str2num(get(handles.edit3,'String'));
handles.rotation = max(min(handles.rotation,get(handles.slider3,'Max')),get(handles.slider3,'Min'));
set(handles.slider3,'Value',handles.rotation);
set(handles.edit3,'String',num2str(handles.rotation));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
handles.scale = str2num(get(handles.edit4,'String'));
handles.scale = max(min(handles.scale,get(handles.slider4,'Max')),get(handles.slider4,'Min'));
set(handles.slider4,'Value',handles.scale);
set(handles.edit4,'String',num2str(handles.scale));
handles.fixedim_resize = imresize(handles.fixedim,handles.scale);
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushWiggleX.
function pushWiggleX_Callback(hObject, eventdata, handles)
set(handles.pushResetX,'Enable','on')
handles.Xpos_prev = handles.Xpos; % save previous
%wiggle
set(handles.text8,'String','Wiggling...');
moving = handles.movingim';
moving_overlay = 0*handles.fixedim_resize;
fixed = handles.fixedim_resize;
[xx, yy] = ndgrid(1:size(moving_overlay,1),1:size(moving_overlay,2));
xshift = -11:11; R=0*xshift;
for ix=1:length(xshift)
    moving_overlay = 0*moving_overlay;
    xOff = handles.Ypos*handles.scale-size(moving,2)/2;
    yOff = (handles.Xpos+xshift(ix))*handles.scale-size(moving,1)/2;
    mapX = round((xx-xOff)*cosd(handles.rotation) + (yy-yOff)*sind(handles.rotation) +size(moving,2)/2);
    mapY = round((yy-yOff)*cosd(handles.rotation) - (xx-xOff)*sind(handles.rotation) +size(moving,1)/2);
    zeromask = mapX<1 | mapX>size(moving,2) | mapY<1 | mapY>size(moving,1);
    mapX(zeromask) = 0; mapY(zeromask)=0;
    useind = find(~zeromask);
    mapInd = sub2ind(size(moving),mapY(useind),mapX(useind));
    moving_overlay(useind) = moving(mapInd);
    R(ix) = corr2(fixed,moving_overlay);
    [ix R]
end
ind = find(R==max(R),1);
handles.Xpos = handles.Xpos+xshift(ind);
set(handles.text8,'String','Done wiggling');
set(handles.slider1,'Value',handles.Xpos);
set(handles.edit1,'String',num2str(handles.Xpos));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)

% --- Executes on button press in pushWiggleY.
function pushWiggleY_Callback(hObject, eventdata, handles)
set(handles.pushResetY,'Enable','on')
handles.Ypos_prev = handles.Ypos; % save previous
%wiggle
set(handles.text8,'String','Wiggling...');
moving = handles.movingim';
moving_overlay = 0*handles.fixedim_resize;
fixed = handles.fixedim_resize;
[xx, yy] = ndgrid(1:size(moving_overlay,1),1:size(moving_overlay,2));
yshift = -11:11; R=0*yshift;
for iy=1:length(yshift)
    moving_overlay = 0*moving_overlay;
    xOff = (handles.Ypos+yshift(iy))*handles.scale-size(moving,2)/2;
    yOff = handles.Xpos*handles.scale-size(moving,1)/2;
    mapX = round((xx-xOff)*cosd(handles.rotation) + (yy-yOff)*sind(handles.rotation) +size(moving,2)/2);
    mapY = round((yy-yOff)*cosd(handles.rotation) - (xx-xOff)*sind(handles.rotation) +size(moving,1)/2);
    zeromask = mapX<1 | mapX>size(moving,2) | mapY<1 | mapY>size(moving,1);
    mapX(zeromask) = 0; mapY(zeromask)=0;
    useind = find(~zeromask);
    mapInd = sub2ind(size(moving),mapY(useind),mapX(useind));
    moving_overlay(useind) = moving(mapInd);
    R(iy) = corr2(fixed,moving_overlay);
    [iy R]
end
ind = find(R==max(R),1);
handles.Ypos = handles.Ypos+yshift(ind);
set(handles.text8,'String','Done wiggling');
set(handles.slider2,'Value',handles.Ypos);
set(handles.edit2,'String',num2str(handles.Ypos));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes on button press in pushWiggleR.
function pushWiggleR_Callback(hObject, eventdata, handles)
set(handles.pushResetR,'Enable','on')
handles.rotation_prev = handles.rotation; % save previous
%wiggle
set(handles.text8,'String','Wiggling...');
moving = handles.movingim';
moving_overlay = 0*handles.fixedim_resize;
fixed = handles.fixedim_resize;
[xx, yy] = ndgrid(1:size(moving_overlay,1),1:size(moving_overlay,2));
rshift = -5:0.5:5; R=0*rshift;
for ir=1:length(rshift)
    r = handles.rotation+rshift(ir)
    moving_overlay = 0*moving_overlay;
    xOff = handles.Ypos*handles.scale-size(moving,2)/2;
    yOff = handles.Xpos*handles.scale-size(moving,1)/2;
    mapX = round((xx-xOff)*cosd(rotation) + (yy-yOff)*sind(handles.rotation) +size(moving,2)/2);
    mapY = round((yy-yOff)*cosd(rotation) - (xx-xOff)*sind(handles.rotation) +size(moving,1)/2);
    zeromask = mapX<1 | mapX>size(moving,2) | mapY<1 | mapY>size(moving,1);
    mapX(zeromask) = 0; mapY(zeromask)=0;
    useind = find(~zeromask);
    mapInd = sub2ind(size(moving),mapY(useind),mapX(useind));
    moving_overlay(useind) = moving(mapInd);
    R(ir) = corr2(fixed,moving_overlay);
    [ir R]
end
ind = find(R==max(R),1);
handles.rotation = handles.rotation+rshift(ind);
set(handles.text8,'String','Done wiggling');
set(handles.slider3,'Value',handles.rotation);
set(handles.edit3,'String',num2str(handles.rotation));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes on button press in pushWiggleS.
function pushWiggleS_Callback(hObject, eventdata, handles)
set(handles.pushResetS,'Enable','on')
handles.scale_prev = handles.scale; % save previous
%wiggle
set(handles.text8,'String','Wiggling...');
moving = handles.movingim';
moving_overlay = 0*handles.fixedim_resize;
fixed = handles.fixedim_resize;
[xx, yy] = ndgrid(1:size(moving_overlay,1),1:size(moving_overlay,2));
yshift = -11:11; R=0*yshift;
for iy=1:length(yshift)
    moving_overlay = 0*moving_overlay;
    xOff = (handles.Ypos+yshift(iy))*handles.scale-size(moving,2)/2;
    yOff = handles.Xpos*handles.scale-size(moving,1)/2;
    mapX = round((xx-xOff)*cosd(handles.rotation) + (yy-yOff)*sind(handles.rotation) +size(moving,2)/2);
    mapY = round((yy-yOff)*cosd(handles.rotation) - (xx-xOff)*sind(handles.rotation) +size(moving,1)/2);
    zeromask = mapX<1 | mapX>size(moving,2) | mapY<1 | mapY>size(moving,1);
    mapX(zeromask) = 0; mapY(zeromask)=0;
    useind = find(~zeromask);
    mapInd = sub2ind(size(moving),mapY(useind),mapX(useind));
    moving_overlay(useind) = moving(mapInd);
    R(iy) = corr2(fixed,moving_overlay);
    [iy R]
end
ind = find(R==max(R),1);
handles.Ypos = handles.Ypos+yshift(ind);
set(handles.text8,'String','Done wiggling');
set(handles.slider2,'Value',handles.Ypos);
set(handles.edit2,'String',num2str(handles.Ypos));
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)


% --- Executes on button press in pushResetX.
function pushResetX_Callback(hObject, eventdata, handles)
handles.Xpos = handles.Xpos_prev;
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)

% --- Executes on button press in pushResetY.
function pushResetY_Callback(hObject, eventdata, handles)
handles.Ypos = handles.Ypos_prev;
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)

% --- Executes on button press in pushResetR.
function pushResetR_Callback(hObject, eventdata, handles)
handles.rotation = handles.rotation_prev;
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)

% --- Executes on button press in pushResetS.
function pushResetS_Callback(hObject, eventdata, handles)
handles.scale = handles.scale_prev;
guidata(hObject,handles)
update_axes(hObject, eventdata, handles)
