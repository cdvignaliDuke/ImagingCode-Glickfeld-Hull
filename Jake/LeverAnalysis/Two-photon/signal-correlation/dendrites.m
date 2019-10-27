function varargout = dendrites(varargin)
% DENDRITES M-file for dendrites.fig
% Copyright 2006 Ilker Ozden&Megan R Sullivan.
%      DENDRITES, by itself, creates a new DENDRITES or raises the existing
%      singleton*.
%
%      H = DENDRITES returns the handle to a new DENDRITES or the handle to
%      the existing singleton*.
%
%      DENDRITES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DENDRITES.M with the given input arguments.
%
%      DENDRITES('Property','Value',...) creates a new DENDRITES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dendrites02_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dendrites_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dendrites

% Last Modified by GUIDE v2.5 14-Jul-2009 18:04:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dendrites_OpeningFcn, ...
                   'gui_OutputFcn',  @dendrites_OutputFcn, ...
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


% --- Executes just before dendrites is made visible.
function dendrites_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dendrites (see VARARGIN)

% Choose default command line output for dendrites
handles.output = hObject;

handles.folder=['C:\data\'];
handles.mainmask=zeros(64,64,2);
handles.maskcnt=1;
handles.filtercounter=zeros(100,1);
handles.totaldendmask=zeros(64,64);
handles.dendfrnummat=0;
set(handles.brightness,'string',0.1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dendrites wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dendrites_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function fr_slider_Callback(hObject, eventdata, handles)
% hObject    handle to fr_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


b=(str2num(get(handles.brightness,'string')))/10;
frinit=str2num(get(handles.frinit,'string'));
mov=handles.mov;

frame_num = round(get(handles.fr_slider,'Value'));
set(handles.fr_slider,'value',frame_num);
set(handles.fr_num,'string',num2str(frame_num+frinit-1));
   
frame1=b*mov(:,:,:,frame_num);
axes(handles.axes1);
image(frame1)
colormap(gray)


% --- Executes during object creation, after setting all properties.
function fr_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function fr_num_Callback(hObject, eventdata, handles)
% hObject    handle to fr_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fr_num as text
%        str2double(get(hObject,'String')) returns contents of fr_num as a double


% --- Executes during object creation, after setting all properties.
function fr_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cfdf=get(handles.cfdfile,'value');
handles.filtercounter=zeros(100,1);
if cfdf==0;
    frinit=str2num(get(handles.frinit,'string'));
    frfinal=str2num(get(handles.frfinal,'string'));

    filterspec = {'*.tif','TIF files (.tif)';'*.mat','MAT files (.mat)'};
    title = 'load';
    folder=handles.folder; 
    [filename, path, index] = uigetfile (filterspec, title, folder);
    eval(['addpath' ' ' char(39) path char(39) ';']);

    folder=path;
    total_filename=[folder filename];
    
    
    movie_info=imfinfo(total_filename);
    frames_size=size(movie_info);
    frames=frames_size(2);
    
    try
        mi1=getfield(movie_info(1,1),'ImageDescription'); 
        mi2=findstr(mi1,'zoomones');

        zoomfac=mi1(1,mi2+9);    
        set(handles.zoom,'string',zoomfac);
    catch
    end

    if frfinal>frames;
        frfinal=frames;
    end

    frames=frfinal-frinit+1;

    set(handles.fr_num,'String',num2str(frfinal));
    pixels_y=movie_info(1,1).Width;
    pixels_x=movie_info(1,1).Height;



    if frames==1;
        set(handles.fr_slider,'max',frames);
        set(handles.fr_slider,'min',0);
        set(handles.fr_slider,'SliderStep',[1 1]);
    else
        set(handles.fr_slider,'max',frames);
        set(handles.fr_slider,'min',0);
        set(handles.fr_slider,'SliderStep',[1/frames 1/frames]);
    end

    size_filename=size(filename);
    fname=filename(1:size_filename(2)-4);
    handles.fn=filename((size_filename(2)-6):(size_filename(2)-4));

    mov=zeros(pixels_x,pixels_y,1,frames);
      
    if index == 1;
        set(handles.filename, 'String', fname);
        for frame=1:(frfinal-frinit+1);
            [mov(:,:,:,frame),map] = imread(fname, 'tiff', (frame+frinit-1));
        end
        frame1 = mov(:,:,:,1);
        set(handles.fr_slider,'value',1);
        cla(handles.axes1);
        axes(handles.axes1);
        imagesc(frame1)
        colormap(gray)
    end
end

if cfdf==1;
    filterspec = {'*.cfd','cfd files (.cfd)'};
    title = 'load';
    folder=handles.folder; 
    temp_file=folder;


    [filename, path, index] = uigetfile (filterspec, title, temp_file);
    eval(['addpath' ' ' char(39) path char(39) ';']);

    folder=path;

    total_filename=[folder filename];

    hd = headerreader(total_filename);
    frames=hd.frames;


    set(handles.fr_num,'String',num2str(frames));
    set(handles.fr_num2,'String',num2str(frames));

    set(handles.ismask,'value',0);


    pixels_x=hd.x;
    pixels_y=hd.y;

    if frames==1;
        set(handles.fr_slider,'max',frames);
        set(handles.fr_slider,'min',0);
        set(handles.fr_slider,'SliderStep',[1 1]);
    else
        set(handles.fr_slider,'max',frames);
        set(handles.fr_slider,'min',0);
        set(handles.fr_slider,'SliderStep',[1/frames 1/frames]);
    end

    set(handles.brightness,'string',0.1);

    size_filename=size(filename);
    fname=filename(1:size_filename(2)-4);
    
    handles.fn=filename((size_filename(2)-6):(size_filename(2)-4));    
    if index == 1;
        fid=fopen(total_filename);
        [A,count]=fread(fid,'uint8');
        set(handles.filename, 'String', fname);    
        if frames==1;
            data=A((hd.size+1):pixels_x*pixels_y+(hd.size));
        else        
            data=A((hd.size+1):(pixels_x*pixels_y*(frames)+(hd.size))); %removes file header
        end
    
        mov=reshape(data,pixels_x,pixels_y,1,frames); %reshapes for svd

        frame1 = mov(:,:,:,1);
        cla(handles.axes1);
        axes(handles.axes1);
        imagesc(frame1)
        colormap(gray) 

    end

end

cla(handles.axes2);
cla(handles.axes3);

handles.eventcounter=zeros(100,1);
handles.folder=folder;
handles.fname=fname;
handles.mov=mov;
handles.movor=mov;
handles.pixels_x=pixels_x;
handles.pixels_y=pixels_y;
handles.frames=frames;

guidata(hObject, handles);


function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in background.
function background_Callback(hObject, eventdata, handles)
% hObject    handle to background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mov=handles.mov;
frames=handles.frames;

mov=double(mov);
movie_mean1fr=mean(mov,4);

movie_mean=zeros(handles.pixels_x,handles.pixels_y,1,handles.frames);

for i=1:frames;
    movie_mean(:,:,1,i)=movie_mean1fr(:,:);
end

delF=get(handles.deltaF,'value');

if delF==1;
    movie_ac=100*((mov-movie_mean)./movie_mean);
end

if delF==0;
    movie_ac=(mov-movie_mean);
end

frm=get(handles.fr_slider,'value');
frame1 = movie_ac(:,:,:,frm);

cla(handles.axes1);
axes(handles.axes1);
imagesc(frame1)
colormap(gray)

handles.mov=movie_ac;

guidata(hObject, handles);


% --- Executes on button press in filtermovie.
function filtermovie_Callback(hObject, eventdata, handles)
% hObject    handle to filtermovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mov=handles.mov;
th1=str2num(get(handles.lfc,'string'));
th2=str2num(get(handles.hfc,'string'));

x=handles.pixels_x;
y=handles.pixels_y;
frames=handles.frames;

for i=1:x;
    for j=1:y;
        tr(1:frames)=mov(i,j,1,1:frames);
        trmax=min(tr(frames-10:frames));
        trmin=min(tr(1:10));
        tr(1:frames)=tr(1:frames)-(trmin+(0:frames-1)*((trmax-trmin)/frames));
        mov(i,j,1,1:frames)=tr(1:frames);
    end
end

mov=reshape(mov,x*y,frames);
mov=double(mov');
movfft=fft(mov);

if th1~=0;
    movfft(1:th1,:)=0;
    movfft((frames-th1+1):frames,:)=0;
end

if th2~=0;
    movfft((ceil(frames/2-th2)):(floor(frames/2+th2)),:)=0;
end

mov=real(ifft(movfft));

mov=mov';
mov=reshape(mov,x,y,1,frames);
handles.mov=mov;
guidata(hObject, handles);


function lfc_Callback(hObject, eventdata, handles)
% hObject    handle to lfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lfc as text
%        str2double(get(hObject,'String')) returns contents of lfc as a double


% --- Executes during object creation, after setting all properties.
function lfc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hfc_Callback(hObject, eventdata, handles)
% hObject    handle to hfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hfc as text
%        str2double(get(hObject,'String')) returns contents of hfc as a double


% --- Executes during object creation, after setting all properties.
function hfc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in rawmovie.
function rawmovie_Callback(hObject, eventdata, handles)
% hObject    handle to rawmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rawmovie

handles.mov=handles.movor;

guidata(hObject, handles);


% --- Executes on button press in chooseregion.
function chooseregion_Callback(hObject, eventdata, handles)
% hObject    handle to chooseregion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frames=handles.frames;
pixels_x=handles.pixels_x;
pixels_y=handles.pixels_y;

axes(handles.axes1);
reg=roipoly;
mov_fltrd=handles.mov;

realreg=zeros(pixels_x,pixels_y);
    

[reg_x,reg_y]=find(reg==1);
pixel_number=length(reg_x);
    
xmax=max(reg_x);
ymax=max(reg_y);
xmin=min(reg_x);
ymin=min(reg_y);

movie_cut=zeros(xmax-xmin+1,ymax-ymin+1,1,frames);
movie_reg=zeros(pixels_x,pixels_y,1,frames);

xsize=xmax-xmin+1;
ysize=ymax-ymin+1;

for i=1:pixel_number;
        movie_reg(reg_x(i),reg_y(i),1,:)=mov_fltrd(reg_x(i),reg_y(i),1,:);
end

for i=xmin:xmax;
    for j=ymin:ymax;
        movie_cut(i+1-xmin,j+1-ymin,1,:)=movie_reg(i,j,1,:);
    end
end

fr_num=str2num(get(handles.fr_num,'string'));
b=str2num(get(handles.brightness,'string'));
frame2=b*movie_reg(:,:,:,fr_num);

cla(handles.axes2);
axes(handles.axes2);
imagesc(frame2)
colormap(gray)

set(handles.fr_slider2,'value',0);
set(handles.fr_num2,'string',num2str(0));


for i=1:frames;
    t=sum(movie_cut(:,:,1,i));
    trace(i)=sum(t)/pixel_number;
end

time(1:frames)=2*pixels_x*(1:frames)/1000;
time=time';

axes(handles.axes3);
hold on
plot(time,trace)
zoom on

handles.axes3trace=trace;

handles.regpixelnumber=pixel_number;
handles.xsize=xsize;
handles.ysize=ysize;
handles.xmax=xmax;
handles.xmin=xmin;
handles.ymax=ymax;
handles.ymin=ymin;
handles.moviecut=movie_cut;
handles.moviereg=movie_reg;

guidata(hObject, handles);



% --- Executes on slider movement.
function fr_slider2_Callback(hObject, eventdata, handles)
% hObject    handle to fr_slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


pixx=handles.pixels_x;
pixy=handles.pixels_y;

mainmask=handles.mainmask;
dendfrnummat=handles.dendfrnummat;
siz=size(mainmask);
frame_num = round(get(handles.fr_slider2,'Value'));
set(handles.fr_slider2,'value',frame_num);
set(handles.fr_num2,'string',num2str(frame_num));
set(handles.dendnum,'string',num2str(frame_num));
if frame_num<siz(3);
    set(handles.vertpos,'string',num2str(handles.posx(frame_num)));
    set(handles.horpos,'string',num2str(handles.posy(frame_num)));
    set(handles.dendfrnum, 'string', num2str(dendfrnummat(frame_num)));
end

tempfr(:,:)=mainmask(:,:,siz(3));
tempfr2(:,:)=mainmask(:,:,frame_num);

if frame_num~=siz(3);

    for i=1:pixx;
        for j=1:pixy;
            if tempfr(i,j)>=2;
                tempfr(i,j)=1;
            end
            if tempfr2(i,j)>=2;
                tempfr2(i,j)=1;
            end
        end
    end
    frame2=tempfr(:,:)+tempfr2(:,:);
elseif frame_num==siz(3);
    frame2=mainmask(:,:,frame_num);
end


axes(handles.axes2);
imagesc(frame2)
colormap(gray)

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function fr_slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




function fr_num2_Callback(hObject, eventdata, handles)
% hObject    handle to fr_num2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fr_num2 as text
%        str2double(get(hObject,'String')) returns contents of fr_num2 as a double


% --- Executes during object creation, after setting all properties.
function fr_num2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_num2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in showregiontrace.
function showregiontrace_Callback(hObject, eventdata, handles)
% hObject    handle to showregiontrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

movie=handles.moviecut;
regpixelnumber=handles.regpixelnumber;

frames=size(movie,4);

for i=1:frames;
    t=sum(movie(:,:,1,i));
    trace(i)=sum(t)/regpixelnumber;
end

pixels_x=handles.pixels_x;
time(1:frames)=2*pixels_x*(1:frames)/1000;
time=time';

axes(handles.axes3);
hold on
plot(time,trace)
zoom on

handles.axes3trace=trace;

guidata(hObject, handles);




% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes3);



% --- Executes on button press in finddend.
function finddend_Callback(hObject, eventdata, handles)
% hObject    handle to finddend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mov=handles.moviecut;
th1=str2num(get(handles.corrth,'string'));
th2=str2num(get(handles.nth,'string'));
frames=handles.frames;
xsize=handles.xsize;
ysize=handles.ysize;
xmin=handles.xmin;
xmax=handles.xmax;
ymin=handles.ymin;
ymax=handles.ymax;
pixels_x=handles.pixels_x;
pixels_y=handles.pixels_y;

movie_size=size(mov);
sizex=movie_size(1);
sizey=movie_size(2);

mov=reshape(mov,sizex*sizey,frames);

mov=mov';
cormat=corrcoef(mov);
cormat_size=size(cormat,1);


for i=1:cormat_size;
    for j=1:cormat_size;
        if cormat(i,j)>=th1;
            cormat(i,j)=1;
        else
            cormat(i,j)=0;
        end
    end
end

cnt=1;

for i=1:cormat_size;
    t=sum(cormat(i,:));
    if t>=th2;
        cor_pix(cnt)=i;
        cnt=cnt+1;
    end
end

handles.dendpixelnumber=cnt-1;

movie_mask=zeros(frames,xsize*ysize);
movie_tempmask=zeros(xsize*ysize,1);

movie_mask(:,cor_pix(1:max(size(cor_pix))))=mov(:,cor_pix(1:max(size(cor_pix))));
movie_tempmask(cor_pix(1:max(size(cor_pix))))=1;

movie_dend=reshape(movie_mask',xsize,ysize,1,frames);    
movie_tempmask=reshape(movie_tempmask,xsize,ysize);  

dendrite_mask=zeros(pixels_x,pixels_y);

dendfrnum=get(handles.fr_num,'string');
set(handles.dendfrnum,'string',dendfrnum);

dendrite_mask(xmin:xmax,ymin:ymax)=movie_tempmask(:,:);
handles.moviedend=movie_dend;
handles.dendmask=dendrite_mask;

cla(handles.axes2);
axes(handles.axes2);
imagesc(dendrite_mask)
colormap(gray)


guidata(hObject, handles);


function corrth_Callback(hObject, eventdata, handles)
% hObject    handle to corrth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrth as text
%        str2double(get(hObject,'String')) returns contents of corrth as a double


% --- Executes during object creation, after setting all properties.
function corrth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nth_Callback(hObject, eventdata, handles)
% hObject    handle to nth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nth as text
%        str2double(get(hObject,'String')) returns contents of nth as a double


% --- Executes during object creation, after setting all properties.
function nth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in dendtrace.
function dendtrace_Callback(hObject, eventdata, handles)
% hObject    handle to dendtrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dendmask=handles.dendmask;
mov=handles.mov;

dendpixelnumber=handles.dendpixelnumber;

px=handles.pixels_x;
py=handles.pixels_y;
frames=handles.frames;
m=zeros(px,py,1,frames);

for i=1:px;
    for j=1:py;
        if dendmask(i,j)==1;
            m(i,j,1,:)=mov(i,j,1,:);
        end
    end
end

trc=zeros(frames,1);

for i=1:frames;
    t=sum(m(:,:,1,i));
    trc(i)=sum(t)/dendpixelnumber;
end

frames=length(trc);
pixels_x=handles.pixels_x;
time(1:frames)=2*pixels_x*(1:frames)/1000;
time=time';

axes(handles.axes3);
hold on
plot(time,trc,'r')
zoom on

handles.dendtrc=trc;
handles.axes3trace=trc;

guidata(hObject, handles);



% --- Executes on button press in savetrace.
function savetrace_Callback(hObject, eventdata, handles)
% hObject    handle to savetrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pixels_x=handles.pixels_x;
data=handles.axes3trace;
data=data';
frames=length(data);
time(1:frames)=2*pixels_x*(1:frames)/1000;
time=time';
filterspec = {'*.mat','mat files (.mat)'};
title = 'Save';
temp_file=[handles.fname 'tr'];
filename=handles.fname;

[filename, path, index] = uiputfile (filterspec, title, temp_file);
    
if index == 1
    save([path filename '.mat'],'data','time');
end

guidata(hObject, handles);



% --- Executes on button press in adddendrite.
function adddendrite_Callback(hObject, eventdata, handles)
% hObject    handle to adddendrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles.totaldendmask=handles.dendmask+handles.totaldendmask;
pixels_x=handles.pixels_x;
pixels_y=handles.pixels_y;
frames=handles.frames;
maskcnt=handles.maskcnt;

if maskcnt~=1;
    tempmask=handles.mainmask;
    tempdendfrnummat=handles.dendfrnummat;
end

if maskcnt==1;
    tempmask=zeros(pixels_x,pixels_y,maskcnt);
    tempdendfrnummat=0;
end

mainmask=zeros(pixels_x,pixels_y,maskcnt);
dendfrnummat=zeros(maskcnt,1);

if maskcnt~=1;
    for i=1:maskcnt-1;
        mainmask(:,:,i)=tempmask(:,:,i);
        dendfrnummat(i)=tempdendfrnummat(i);
    end
end

dendfrnummat(maskcnt)=str2num(get(handles.fr_num,'string'));
mainmask(:,:,maskcnt)=handles.dendmask(:,:);

handles.posx=zeros(maskcnt+1,1);
handles.posy=zeros(maskcnt+1,1);
handles.xLength=zeros(maskcnt+1,1);
handles.yLength=zeros(maskcnt+1,1);
handles.numpixel=zeros(maskcnt+1,1);

totaldendmask=zeros(pixels_x,pixels_y);

for i=1:maskcnt;
    totaldendmask(:,:)=totaldendmask(:,:)+mainmask(:,:,i);    
end


inds=find(totaldendmask>2);
totaldendmask(inds)=2;

totalmasktemp(:,:,1)=totaldendmask(:,:);
mainmask=cat(3,mainmask,totalmasktemp);

dendfrnummat(maskcnt)=str2num(get(handles.fr_num,'string'));

set(handles.fr_slider2,'max',(maskcnt+1));
set(handles.fr_slider2,'min',0);
set(handles.fr_slider2,'SliderStep',[1/(maskcnt+1) 1/(maskcnt+1)]);
set(handles.fr_num2,'String',num2str(maskcnt+1));

cla(handles.axes2);
axes(handles.axes2);
imagesc(totaldendmask)
colormap(gray)

handles.maskcnt=maskcnt+1;
set(handles.fr_slider2,'value',(maskcnt+1));  

masksize=size(mainmask);
numdendrites=masksize(3)-1;
posx=zeros(numdendrites,1);
posy=zeros(numdendrites,1);
xLength=zeros(numdendrites,1);
yLength=zeros(numdendrites,1);
numpixel=zeros(numdendrites,1);

for i=1:numdendrites;
    dend(:,:)=mainmask(:,:,i);
    [x,y]=find(dend);
    posx(i)=(max(x)+min(x))/2;
    posy(i)=round((max(y)+min(y))/2);
    xLength(i)=max(x)-min(x);
    yLength(i)=max(y)-min(y);    
end

tempmask=(zeros(masksize(1),masksize(2)));
tempdendfrnum=0;
qq=0;   

for i=1:numdendrites;
    [posval,posmin]=min(posx(i:numdendrites));
    if posmin==numdendrites-i+1&qq==0;
        lastdendrite=i;
        qq=qq+1;
    end
    aa=posy(posmin+i-1);
    posy(posmin+i-1)=posy(i);
    posy(i)=aa;
    
    bb=yLength(posmin+i-1);
    yLength(posmin+i-1)=yLength(i);
    yLength(i)=bb;
    
    cc=xLength(posmin+i-1);
    xLength(posmin+i-1)=xLength(i);
    xLength(i)=cc;
    
    posx(posmin+i-1)=posx(i);
    posx(i)=posval;
    tempmask(:,:)=mainmask(:,:,posmin+i-1);
    mainmask(:,:,posmin+i-1)=mainmask(:,:,i);
    mainmask(:,:,i)=tempmask(:,:);
    
    tempdendfrnum=dendfrnummat(posmin+i-1);
    dendfrnummat(posmin+i-1)=dendfrnummat(i);
    dendfrnummat(i)=tempdendfrnum;
end

handles.posy=posy;
handles.posx=posx;
handles.mainmask=mainmask;
handles.lastdendrite=lastdendrite;
handles.xLength=xLength;
handles.yLength=yLength;
handles.numpixel=numpixel;
handles.dendfrnummat=dendfrnummat;

guidata(hObject, handles);


% --- Executes on button press in clearalldendrites.
function clearalldendrites_Callback(hObject, eventdata, handles)
% hObject    handle to clearalldendrites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.totaldendmask=zeros(handles.pixels_x,handles.pixels_y);
handles.mainmask=zeros(handles.pixels_x,handles.pixels_y,1);
handles.dendfrnummat=0;
handles.maskcnt=1;
handles.eventcounter=zeros(100,1);

guidata(hObject, handles);

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
% hObject    handle to show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show
frinit=str2num(get(handles.frinit,'string'));
a=get(hObject,'Value');
s=size(handles.mainmask);

frnum=str2num(get(handles.fr_num2,'string'));

if frnum==0;
    dendmask=handles.dendmask;
else
    dendmask(:,:)=handles.mainmask(:,:,frnum);
end
mov=handles.mov;

m_size=size(mov);
px=m_size(1);
py=m_size(2);
frames=m_size(4);
m=zeros(m_size(1),m_size(2),1,m_size(4));

b=max(max(max(mov)));

if a==1;

        for i=1:px;
            for j=1:py;
                %for k=1:frames;
                    if dendmask(i,j)==0;
                        m(i,j,1,:)=mov(i,j,1,:);
                    end
                    if dendmask(i,j)~=0;
                        m(i,j,1,:)=b;
                    end
                %end
            end
        end
        handles.m=m;

end
b=(str2num(get(handles.brightness,'string')))/10;

if a==0;
    mov=handles.mov;
end
if a==1;
    mov=handles.m;
end


frame_num = round(get(handles.fr_slider,'Value'));
set(handles.fr_slider,'value',frame_num);
set(handles.fr_num,'string',num2str(frame_num+frinit-1));

frame1=b*mov(:,:,:,frame_num);
axes(handles.axes1);
image(frame1)
colormap(gray)


guidata(hObject, handles);


% --- Executes on button press in savealldendrites.
function savealldendrites_Callback(hObject, eventdata, handles)
% hObject    handle to savealldendrites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainmask=handles.mainmask;
dendfrnummat=handles.dendfrnummat;
filterspec = {'*.mat','mat files (.mat)'; '*.dat','dat files (.dat)'};
title = 'Save';
temp_file=[handles.fname 'dendriteall'];
filename=handles.fname;

positions=handles.posx;

xLen=handles.xLength;
yLen=handles.yLength;
numpixels=handles.numpixel;

aaa=[xLen yLen numpixels];

[filename, path, index] = uiputfile (filterspec, title, temp_file);

    
if index == 1
     save([path filename '.mat'],'mainmask');
     save([path filename 'dendfrnum.mat'],'dendfrnummat');
     eval(['save ' path filename 'positions.dat' ' positions -ascii']);
     eval(['save ' path filename 'dendsize.dat' ' aaa -ascii']);
end

guidata(hObject, handles);



function brightness_Callback(hObject, eventdata, handles)
% hObject    handle to brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of brightness as text
%        str2double(get(hObject,'String')) returns contents of brightness as a double


% --- Executes during object creation, after setting all properties.
function brightness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in filtertrace.
function filtertrace_Callback(hObject, eventdata, handles)
% hObject    handle to filtertrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filtertrace


tr=handles.axes3trace;

frames=length(tr);
th1=str2num(get(handles.filtertracelfc,'string'));
th2=str2num(get(handles.filtertracehfc,'string'));

trmax=min(tr(frames-10:frames));
trmin=min(tr(1:10));

for i=1:frames;
    tr(i)=tr(i)-(trmin+(i-1)*((trmax-trmin)/frames));
end

trfft=fft(tr);

if th1~=0;
    trfft(1:th1)=0;
    trfft((frames-th1+1):frames)=0;
end

if th2~=0;
    trfft((frames/2-th2):(frames/2+th2))=0;
end

tr_fltrd=real(ifft(trfft));

pixels_x=handles.pixels_x;
time(1:frames)=2*pixels_x*(1:frames)/1000;
time=time';

cla(handles.axes4);
axes(handles.axes4);
plot(time,tr_fltrd);
zoom on
handles.axes4trace=tr_fltrd;

guidata(hObject, handles);


function filtertracelfc_Callback(hObject, eventdata, handles)
% hObject    handle to filtertracelfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filtertracelfc as text
%        str2double(get(hObject,'String')) returns contents of filtertracelfc as a double


% --- Executes during object creation, after setting all properties.
function filtertracelfc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filtertracelfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filtertracehfc_Callback(hObject, eventdata, handles)
% hObject    handle to filtertracehfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filtertracehfc as text
%        str2double(get(hObject,'String')) returns contents of filtertracehfc as a double


% --- Executes during object creation, after setting all properties.
function filtertracehfc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filtertracehfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deltaF.
function deltaF_Callback(hObject, eventdata, handles)
% hObject    handle to deltaF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deltaF





function frinit_Callback(hObject, eventdata, handles)
% hObject    handle to frinit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frinit as text
%        str2double(get(hObject,'String')) returns contents of frinit as a double


% --- Executes during object creation, after setting all properties.
function frinit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frinit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frfinal_Callback(hObject, eventdata, handles)
% hObject    handle to frfinal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frfinal as text
%        str2double(get(hObject,'String')) returns contents of frfinal as a double


% --- Executes during object creation, after setting all properties.
function frfinal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frfinal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in saveproc.
function saveproc_Callback(hObject, eventdata, handles)
% hObject    handle to saveproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


mov=uint16(handles.mov);
filterspec = {'*.tif','tif files (.tif)'};
title = 'Save';
temp_file=[handles.fname 'processed'];
filename=handles.fname;
frames=size(mov,4);


[filename, path, index] = uiputfile (filterspec, title, temp_file);
if index == 1
    imwrite(mov(:,:,:,1),[path filename '.tif'],'TIF','compression', 'none')
    for i=2:frames;   
        imwrite(mov(:,:,:,i),[path filename '.tif'],'TIF','WriteMode','append','compression','none');
    end
end

guidata(hObject, handles);




% --- Executes on button press in loadmask.
function loadmask_Callback(hObject, eventdata, handles)
% hObject    handle to loadmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filterspec = {'*.mat','MAT files (.mat)'};
title = 'load';
folder='C:\matlabR2007a\work\'; 


[filename, path, index] = uigetfile (filterspec, title, folder);
eval(['addpath' ' ' char(39) path char(39) ';']);

folder=path;

total_filename=[folder filename];
fnlength=length(total_filename)
  
if index == 1;
    eval(['load ' total_filename]);
    try
        mainmask=dendmask;
    catch
    end
    dendfn=[total_filename(1:fnlength-4) 'dendfrnum.mat'];
    if exist(dendfn)==2;   
        eval(['load ' total_filename(1:fnlength-4) 'dendfrnum.mat']);
    else
        a=size(mainmask);
        dendfrnummat=zeros(a(3)-1,1);
    end
end



mainmask_size=size(mainmask);
numdend=mainmask_size(3)-1;
set(handles.dendnum,'string',num2str(numdend));

totaldendmask(:,:)=mainmask(:,:,(numdend+1));

set(handles.fr_slider2,'max',(numdend+1));
set(handles.fr_slider2,'min',0);
set(handles.fr_slider2,'SliderStep',[1/(numdend+1) 1/(numdend+1)]);

set(handles.fr_num2,'String',num2str(numdend+1));
set(handles.fr_slider2,'Value',(numdend+1));

cla(handles.axes2);
axes(handles.axes2);
imagesc(totaldendmask)
colormap(gray)

handles.dendfrnummat=dendfrnummat;
handles.eventcounter=zeros(100,1);
handles.inspkcounter=zeros(100,1);
handles.posx=zeros(numdend,1);
handles.posy=zeros(numdend,1);
handles.xLength=zeros(numdend,1);
handles.yLength=zeros(numdend,1);
handles.numpixel=zeros(numdend,1);
handles.maskcnt=numdend+1;
handles.numdend=numdend;
handles.mainmask=mainmask;

guidata(hObject, handles);




function dendnum_Callback(hObject, eventdata, handles)
% hObject    handle to dendnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dendnum as text
%        str2double(get(hObject,'String')) returns contents of dendnum as a double


% --- Executes during object creation, after setting all properties.
function dendnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dendnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in dFtraces.
function dFtraces_Callback(hObject, eventdata, handles)
% hObject    handle to dFtraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainmask=handles.mainmask;
mov=double(handles.movor);

handles.filtercounter=zeros(100,1);
handles.inspkcounter=zeros(100,1);
handles.rvrsinspkcounter=zeros(100,1);
a=size(mainmask);
numdend=a(3)-1;

m_size=size(mov);
px=m_size(1);
py=m_size(2);
frames=m_size(4);

tracemat=zeros(numdend,frames);
trc=zeros(1,frames);

for s=1:numdend;   
    [x,y]=find(mainmask(:,:,s));
    cntx=length(x);
    m=zeros(cntx,1,frames);   
    for i=1:cntx;        
            m(i,1,:)=mov(x(i),y(i),1,:);
    end
    for i=1:frames;
        trc(i)=sum(m(:,1,i));                        
    end
    movorsum=sum(trc);    
    trc(:)=trc(:)*frames/movorsum-1;    
    tracemat(s,:)=trc(:);

end

handles.numdend=numdend;
handles.tracemat=tracemat;
handles.tracematfltrd=zeros(numdend,frames);
handles.eventsmat=zeros(numdend,frames);
handles.eventsmatcorrected=zeros(numdend,frames);
guidata(hObject, handles);



% --- Executes on button press in showdendtrace.
function showdendtrace_Callback(hObject, eventdata, handles)
% hObject    handle to showdendtrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dendnumb=str2num(get(handles.dendnum,'String'));

set(handles.fr_slider2,'value',dendnumb);
set(handles.fr_num2,'string',num2str(dendnumb));

frame2=handles.mainmask(:,:,dendnumb);
axes(handles.axes2);
imagesc(frame2)
colormap(gray)

tracemat=handles.tracemat;
trc=tracemat(dendnumb,:);
 
cla(handles.axes3);
axes(handles.axes3);
plot(trc,'r')
zoom on

handles.axes3trc=trc;

guidata(hObject, handles);


% --- Executes on button press in savealldendtraces.
function savealldendtraces_Callback(hObject, eventdata, handles)
% hObject    handle to savealldendtraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filtercounter=handles.filtercounter;
pixels_x=handles.pixels_x;
frames=handles.frames;
a=size(handles.mainmask);
numdend=a(3)-1;
filterspec = {'*.mat','mat files (.mat)'};
title = 'Save';
temp_file=[handles.fname 'tracemat'];
filename=handles.fname;

[filename, path, index] = uiputfile (filterspec, title, temp_file);

fn=handles.fn;


if index == 1
    tracemat=handles.tracemat;
    tracematfltrd=handles.tracematfltrd;
    eval(['m' fn 'tracemat=handles.tracemat;']);
    eval(['m' fn 'tracematfltrd=handles.tracematfltrd;']);      
    save([path filename '.mat'],['m' fn 'tracemat']);
    save([path filename 'fltrd.mat'],['m' fn 'tracematfltrd']);

end
   
guidata(hObject, handles);    
        
  

% --- Executes on button press in filterdendtrace.
function filterdendtrace_Callback(hObject, eventdata, handles)
% hObject    handle to filterdendtrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dendnum=str2num(get(handles.dendnum,'string'));

tr=handles.tracemat(dendnum,:);

frames=length(tr);
th1=str2num(get(handles.dendfilterlfc,'string'));
th2=str2num(get(handles.dendfilterhfc,'string'));

trmax=min(tr(frames-10:frames));
trmin=min(tr(1:10));

tr=tr-(trmin+(0:frames-1)*((trmax-trmin)/frames));

trfft=fft(tr);

if th1~=0;
    trfft(1:th1)=0;
    trfft((frames-th1+1):frames)=0;
end

if th2~=0;
    trfft((ceil(frames/2-th2)):(floor(frames/2+th2)))=0;
end

tr_fltrd=real(ifft(trfft));

cla(handles.axes4);
axes(handles.axes4);
plot(tr_fltrd);
zoom on
handles.axes4trace=tr_fltrd;

handles.tracematfltrd(dendnum,:)=tr_fltrd;
handles.filtercounter(dendnum)=1;

guidata(hObject, handles);

function dendfilterlfc_Callback(hObject, eventdata, handles)
% hObject    handle to dendfilterlfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dendfilterlfc as text
%        str2double(get(hObject,'String')) returns contents of dendfilterlfc as a double





% --- Executes during object creation, after setting all properties.
function dendfilterlfc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dendfilterlfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dendfilterhfc_Callback(hObject, eventdata, handles)
% hObject    handle to dendfilterhfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dendfilterhfc as text
%        str2double(get(hObject,'String')) returns contents of dendfilterhfc as a double


% --- Executes during object creation, after setting all properties.
function dendfilterhfc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dendfilterhfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in removedendrite.
function removedendrite_Callback(hObject, eventdata, handles)
% hObject    handle to removedendrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainmask=handles.mainmask;
dendritenum=str2num(get(handles.dendnum,'string'));
dendfrnummat=handles.dendfrnummat;

a=size(mainmask);

tempmainmask=zeros(a(1),a(2),a(3)-1);
tempdendfrnummat=zeros(a(3)-2,1);

tempmainmask(:,:,1:(dendritenum-1))=mainmask(:,:,1:(dendritenum-1));
tempmainmask(:,:,((dendritenum):(a(3)-2)))=mainmask(:,:,((dendritenum+1):(a(3)-1)));

tempdendfrnummat(1:dendritenum-1)=dendfrnummat(1:dendritenum-1);
tempdendfrnummat(dendritenum:(a(3)-2))=dendfrnummat((dendritenum+1):(a(3)-1));

totaldendmask=zeros(a(1),a(2));

for i=1:(a(3)-2);
    totaldendmask(:,:)=totaldendmask(:,:)+tempmainmask(:,:,i);
end

tempmainmask(:,:,(a(3)-1))=totaldendmask(:,:);
handles.mainmask=tempmainmask;

handles.dendfrnummat=tempdendfrnummat;


set(handles.fr_slider2,'max',((a(3)-1)));
set(handles.fr_slider2,'min',0);
set(handles.fr_slider2,'SliderStep',[1/((a(3)-1)) 1/((a(3)-1))]);

set(handles.fr_num2,'String',num2str((a(3)-1)));
set(handles.fr_slider2,'Value',((a(3)-1)));


cla(handles.axes2);
axes(handles.axes2);
imagesc(totaldendmask)
colormap(gray)

handles.maskcnt=handles.maskcnt-1;

    
set(handles.fr_slider2,'value',(a(3)-1));  

guidata(hObject, handles);


% --- Executes on button press in findevents.
function findevents_Callback(hObject, eventdata, handles)
% hObject    handle to findevents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

whichmethod=get(handles.templatemethod,'value');
pixels_x=handles.pixels_x;
th=str2num(get(handles.eventth,'string'));
dendnum=str2num(get(handles.dendnum,'string'));

if whichmethod==0;
   
    dendnum=str2num(get(handles.dendnum,'string'));
    if handles.filtercounter(dendnum)==1;
        tr=handles.tracematfltrd(dendnum,:);
    else
        tr=handles.tracemat(dendnum,:);
    end
 
    t=handles.pixels_x;
    datasize=length(tr);
    
    i=1;
    cnt=1;

    while i<=datasize;
        if tr(i)>=th;
            while ((i+1<=datasize)&&(tr(i+1)>tr(i)));
                i=i+1;
            end
            peakpos(cnt)=i;
            cnt=cnt+1;
            while ((i+1<=datasize)&&(tr(i+1)<tr(i)) && (tr(i+1)>th));
                i=i+1;
            end
        end
        i=i+1;
    end

    events=zeros(datasize,1);
    events(peakpos(:))=1;

    handles.eventsmat(dendnum,:)=events(:);
    eval(['handles.trc' num2str(dendnum) 'eventindex=peakpos;']);

    handles.eventcounter(dendnum)=1;
    time=2*t*datasize/1000;
    freq=length(peakpos)/time;
    set(handles.eventf,'string',num2str(freq));
    set(handles.eventnum,'string',num2str(length(peakpos)));

    axes(handles.axes3)
    cla
    plot(tr)
    zoom on
    handles.axes3trace=tr;

    axes(handles.axes4)
    cla
    plot(events)
    zoom on
    handles.axes4trace=events;
end


if whichmethod==1;
      
    tr=handles.tr4events;
    datasize=length(tr);
    t=handles.pixels_x;

    i=1;
    cnt=1;

    while i<=datasize;
        if tr(i)>=th;
            while ((i+1<=datasize)&&(tr(i+1)>tr(i)));
                i=i+1;
            end
            peakpos(cnt)=i;
            cnt=cnt+1;
            while ((i+1<=datasize)&&(tr(i+1)<tr(i)) && (tr(i+1)>th));
                i=i+1;
            end
        end
        i=i+1;
    end

    events=zeros(datasize,1);
    events(peakpos(:))=1;


    handles.eventsmat(dendnum,:)=events(:);
    handles.eventcounter(dendnum)=1;

    time=2*t*datasize/1000;

    freq=length(peakpos)/time;
   

    set(handles.eventf,'string',num2str(freq));
    set(handles.eventnum,'string',num2str(length(peakpos)));

    axes(handles.axes3)
    cla
    plot(tr)
    zoom on
    handles.axes3trace=tr;

    axes(handles.axes4)
    cla
    plot(events)
    zoom on
    handles.axes4trace=events;
    
end

handles.events=events;
guidata(hObject, handles);


function eventth_Callback(hObject, eventdata, handles)
% hObject    handle to eventth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eventth as text
%        str2double(get(hObject,'String')) returns contents of eventth as a double


% --- Executes during object creation, after setting all properties.
function eventth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in saveevents.
function saveevents_Callback(hObject, eventdata, handles)
% hObject    handle to saveevents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

eventsmat=handles.eventsmat;

dendpos=handles.posx;

eventcounter=handles.eventcounter;
pixels_x=handles.pixels_x;
frames=handles.frames;
a=size(handles.mainmask);
numdend=a(3)-1;

eventsmatcorrected=zeros(numdend,frames);

for i=1:numdend;
    tr(1,:)=eventsmat(i,:);
    trc=zeros(1,frames+1);
    p=dendpos(i)/pixels_x;
    inds=find(tr>0);
    trc(inds+1)=p;
    trc(inds)=1-p;
    eventsmatcorrected(i,1:frames)=trc(1,2:frames+1);
end
    
    
filterspec = {'*.mat','mat files (.mat)'; '*.dat','dat files (.dat)'};
title = 'Save';
temp_file=handles.fname;
filename=handles.fname;
fn=handles.fn;

[filename, path, index] = uiputfile (filterspec, title, temp_file);

if index == 1;       
        eval(['m' fn 'eventsmat=eventsmat;']);
        save([path filename 'eventsmat.mat'],['m' fn 'eventsmat']);
        eval(['m' fn 'eventsmatcorrected=eventsmatcorrected;']);
        save([path filename 'eventsmatcorrected.mat'],['m' fn 'eventsmatcorrected']);
end
        
guidata(hObject, handles); 



function eventf_Callback(hObject, eventdata, handles)
% hObject    handle to eventf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eventf as text
%        str2double(get(hObject,'String')) returns contents of eventf as a double


% --- Executes during object creation, after setting all properties.
function eventf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function eventnum_Callback(hObject, eventdata, handles)
% hObject    handle to eventnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eventnum as text
%        str2double(get(hObject,'String')) returns contents of eventnum as a double


% --- Executes during object creation, after setting all properties.
function eventnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in finalizedendrites.
function finalizedendrites_Callback(hObject, eventdata, handles)
% hObject    handle to finalizedendrites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainmask=handles.mainmask;
sorthor=get(handles.sorthor,'value');
masksize=size(mainmask);
numdendrites=masksize(3)-1;

posx=zeros(numdendrites,1);
xLength=zeros(numdendrites,1);
yLength=zeros(numdendrites,1);
numpixel=zeros(numdendrites,1);

if sorthor==0;

    for i=1:numdendrites;
        dend(:,:)=mainmask(:,:,i);
        [x,y]=find(dend);
        posx(i)=(max(x)+min(x))/2;
        xLength(i)=max(x)-min(x);
        posy(i)=round((max(y)+min(y))/2);
        yLength(i)=max(y)-min(y);
    end

        tempmask=(zeros(masksize(1),masksize(2)));

    for i=1:numdendrites;
        [posval,posmin]=min(posx(i:numdendrites));
        aa=posy(posmin+i-1);
        posy(posmin+i-1)=posy(i);
        posy(i)=aa;

        bb=yLength(posmin+i-1);
        yLength(posmin+i-1)=yLength(i);
        yLength(i)=bb;

        cc=xLength(posmin+i-1);
        xLength(posmin+i-1)=xLength(i);
        xLength(i)=cc;

        posx(posmin+i-1)=posx(i);
        posx(i)=posval;
        tempmask(:,:)=mainmask(:,:,posmin+i-1);
        mainmask(:,:,posmin+i-1)=mainmask(:,:,i);
        mainmask(:,:,i)=tempmask(:,:);
    end

    for i=1:numdendrites;
        numpixel(i)=sum(sum(mainmask(:,:,i)));
    end

    handles.posy=posy;
    handles.posx=posx;
    handles.mainmask=mainmask;
    handles.xLength=xLength;
    handles.yLength=yLength;
    handles.numpixel=numpixel;
    
elseif sorthor==1;
    for i=1:numdendrites;
        dend(:,:)=mainmask(:,:,i);
        [y,x]=find(dend);
        posx(i)=(max(x)+min(x))/2;
        xLength(i)=max(x)-min(x);
        posy(i)=round((max(y)+min(y))/2);
        yLength(i)=max(y)-min(y);
    end

        tempmask=(zeros(masksize(1),masksize(2)));

    for i=1:numdendrites;
        [posval,posmin]=min(posx(i:numdendrites));
        aa=posy(posmin+i-1);
        posy(posmin+i-1)=posy(i);
        posy(i)=aa;

        bb=yLength(posmin+i-1);
        yLength(posmin+i-1)=yLength(i);
        yLength(i)=bb;

        cc=xLength(posmin+i-1);
        xLength(posmin+i-1)=xLength(i);
        xLength(i)=cc;

        posx(posmin+i-1)=posx(i);
        posx(i)=posval;
        tempmask(:,:)=mainmask(:,:,posmin+i-1);
        mainmask(:,:,posmin+i-1)=mainmask(:,:,i);
        mainmask(:,:,i)=tempmask(:,:);
    end

    for i=1:numdendrites;
        numpixel(i)=sum(sum(mainmask(:,:,i)));
    end

    handles.posy=posx;
    handles.posx=posy;
    handles.mainmask=mainmask;
    handles.xLength=yLength;
    handles.yLength=xLength;
    handles.numpixel=numpixel;
    
end

guidata(hObject, handles); 



function vertpos_Callback(hObject, eventdata, handles)
% hObject    handle to vertpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vertpos as text
%        str2double(get(hObject,'String')) returns contents of vertpos as a double


% --- Executes during object creation, after setting all properties.
function vertpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vertpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function horpos_Callback(hObject, eventdata, handles)
% hObject    handle to horpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of horpos as text
%        str2double(get(hObject,'String')) returns contents of horpos as a double


% --- Executes during object creation, after setting all properties.
function horpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to horpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


function hd = headerreader(filename);

%  This is the same function as headerreader, but if the file does not exist
%  it returns headersize=-1. Will be mainly used in the loops inside of other functions.
%  OPENS 'filename' file and reads cfd header
%  opens "File Open" dialog if 'filename' is not valid name of file
%  RETURNS structure with next fields:
%  .version     Version fo the CfNT program                 uint16
%  .name        Filename                                    char 14
%  .user        User                                        char 16
%  .time        Time                                        char 16
%  .date        Date                                        char 16
%  .zPos        Scan position z                             int32
%  .size        Size of the heqader (byte)                  int32
%  .xPos        Scan position x                             int32
%  .yPos        Scan position y                             int32
%  .x           x-size in pixels                            uint16
%  .y           y-size in pixels                            uint16
%  .frames      number of frames                            int16
%  .scadura     scan duration (msec)                        float32
%  .retrdura    retrace duration (msec)                     float32
%  .zStep       z axis increment in nm                      long (32)
%  .dwellTime	Time in ms to wait during Z/T series scan   int32
%  .avnum		Number of images averaged/step              int16
%  .imtoss  	Number of images to toss before average     int16
%  .AcqTime	    Time between acquisitions                   long (32)
%  .xStep		z axis increment in nm                      long (32)
%  .yStep   	z axis increment in nm                      long (32)
%  .zoom        zoom                                        float32



%   Thess parameters a calculated using other parameters
%  .scanline    scantime per line(msec)                     float32
%  .scanframe   scatime per frame(msec)                     float32

fid=fopen(filename,'r');
if fid==-1                                   % if filename is wrong opens the dailog box
   hd.size=-1;
else

hd.ver=fread(fid,1,'uint16');			
hd.name=char(fread(fid,14,'char')');
hd.user=char(fread(fid,16,'char')');    
hd.time=char(fread(fid,16,'char')');
hd.date=char(fread(fid,16,'char')');

fseek(fid,72,-1);
hd.zPos=fread(fid,1,'int32');
hd.size=fread(fid,1,'int32');                   % ussualy 768
hd.xPos=fread(fid,1,'int32');
hd.yPos=fread(fid,1,'int32');

fseek(fid,280,-1);
hd.x=fread(fid,1,'uint16');			
hd.y=fread(fid,1,'uint16');			

fseek(fid,8,0);							
hd.frames=fread(fid,1,'int16');		

fseek(fid,316,-1);
hd.scandura=fread(fid,1,'float32');
hd.retrdur=fread(fid,1,'float32');
hd.zStep=fread(fid,1,'int32');
hd.dwellTime=fread(fid,1,'uint32');                 % Somehow always 0  (Does not work in version 1.529)

fseek(fid,4,0);	
hd.avnum=fread(fid,1,'uint16');
hd.imtoss=fread(fid,1,'uint16');
hd.AcqTime=fread(fid,1,'int32');                    % Somehow always 0 (Does not work in version 1.529)

fseek(fid,4,0);
hd.xStep=fread(fid,1,'int32');
hd.yStep=fread(fid,1,'int32');

fseek(fid,528,-1);
hd.zoom=fread(fid,1,'float32');				% zoom
fclose(fid);

hd.scanline=hd.scandura+hd.retrdur;
hd.scanframe=hd.scanline*hd.x;
end


% --- Executes on button press in cfdfile.
function cfdfile_Callback(hObject, eventdata, handles)
% hObject    handle to cfdfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cfdfile


% --- Executes on button press in showfltrddendtrace.
function showfltrddendtrace_Callback(hObject, eventdata, handles)
% hObject    handle to showfltrddendtrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dendnumb=str2num(get(handles.dendnum,'String'));
set(handles.fr_slider2,'value',dendnumb);
set(handles.fr_num2,'string',num2str(dendnumb));

frame2=handles.mainmask(:,:,dendnumb);
axes(handles.axes2);
imagesc(frame2)
colormap(gray)

trc=handles.tracematfltrd(dendnumb,:);

cla(handles.axes3);
axes(handles.axes3);
plot(trc,'r')
zoom on

handles.axes3trace=trc;

guidata(hObject, handles);


% --- Executes on button press in removelastdendrite.
function removelastdendrite_Callback(hObject, eventdata, handles)
% hObject    handle to removelastdendrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainmask=handles.mainmask;
dendritenum=handles.lastdendrite;
dendfrnummat=handles.dendfrnummat;

a=size(mainmask);

tempdendfrnummat=zeros(a(3)-1,1);
tempmainmask=zeros(a(1),a(2),a(3)-2);

tempmainmask(:,:,1:(dendritenum-1))=mainmask(:,:,1:(dendritenum-1));
tempmainmask(:,:,((dendritenum):(a(3)-2)))=mainmask(:,:,((dendritenum+1):(a(3)-1)));

tempdendfrnummat(1:(dendritenum-1))=dendfrnummat(1:(dendritenum-1));
tempdendfrnummat(((dendritenum):(a(3)-2)))=dendfrnummat(((dendritenum+1):(a(3)-1)));

totaldendmask=zeros(a(1),a(2));

for i=1:(a(3)-2);
    totaldendmask(:,:)=totaldendmask(:,:)+tempmainmask(:,:,i);
end

tempmainmask(:,:,(a(3)-1))=totaldendmask(:,:);
handles.mainmask=tempmainmask;
handles.dendfrnummat=tempdendfrnummat;

set(handles.fr_slider2,'max',((a(3)-1)));
set(handles.fr_slider2,'min',0);
set(handles.fr_slider2,'SliderStep',[1/((a(3)-1)) 1/((a(3)-1))]);

set(handles.fr_num2,'String',num2str((a(3)-1)));
set(handles.fr_slider2,'Value',((a(3)-1)));

cla(handles.axes2);
axes(handles.axes2);
imagesc(totaldendmask)
colormap(gray)

handles.maskcnt=handles.maskcnt-1;    
set(handles.fr_slider2,'value',(a(3)-1));  
guidata(hObject, handles);


function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoom as text
%        str2double(get(hObject,'String')) returns contents of zoom as a double


% --- Executes during object creation, after setting all properties.
function zoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in usefilter.
function usefilter_Callback(hObject, eventdata, handles)
% hObject    handle to usefilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usefilter




% --- Executes on button press in overlay.
function overlay_Callback(hObject, eventdata, handles)
% hObject    handle to overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


cstrace=handles.axes3trace;
trace=handles.axes4trace;

cstrace=cstrace/(max(cstrace)-min(cstrace));
trace=trace/(max(trace)-min(trace));
trace=trace+1.05+abs(min(trace));

figure
plot(cstrace)
hold on
plot(trace,'r')
zoom on


% --- Executes on button press in applytemplate.
function applytemplate_Callback(hObject, eventdata, handles)
% hObject    handle to applytemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dendnum=str2num(get(handles.dendnum,'string'));

if handles.filtercounter(dendnum)==1;
    tr=handles.tracematfltrd(dendnum,:);
else
    tr=handles.tracemat(dendnum,:);
end

template=handles.finaltemplate;
tempsize=length(template)-1;
templ(1:tempsize)=template(2:(tempsize+1));

datasize=length(tr);
templsize=length(templ);
trmultiply=zeros(datasize,1);

for i=1:(datasize);
    s=0;
    trmax=tr(i);
    if i~=1;
        trampmin=tr(i-1);
        if i==2;
            trmin=min(tr((i-1):min(datasize,(i+templsize-1))));
        else
            trmin=min(tr((i-2):min(datasize,(i+templsize-1))));
        end
        tramp=trmax-trmin;
    else
        trampmin=tr(i);
        trmin=min(tr((i):min(datasize,(i+templsize))));
    end

    for j=1:templsize;
        if j+i-1<=datasize;
            s=s+(tr(i+j-1)-trmin)*templ(j);
        end
    end
    trmultiply(i)=s+trmin;
end


axes(handles.axes3)
cla
plot(tr)
zoom on
handles.axes3trace=tr;

axes(handles.axes4)
cla
plot(trmultiply)
zoom on

handles.axes4trace=trmultiply;
handles.orgtr=tr;
handles.tr4events=trmultiply;
guidata(hObject, handles);


% --- Executes on button press in createtemplate.
function createtemplate_Callback(hObject, eventdata, handles)
% hObject    handle to createtemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


dendnum=str2num(get(handles.dendnum,'string'));

if handles.filtercounter(dendnum)==1;
    tr=handles.tracematfltrd(dendnum,:);
else
    tr=handles.tracemat(dendnum,:);
end

datasize=length(tr);
peaknumber=str2num(get(handles.numberofpeaks,'string'));
tempsize=str2num(get(handles.templatesize,'string'));
set(handles.displayedtemp,'string',num2str(1));

for i=1:peaknumber;
    [peakvalue,peakpos]=max(tr(3:(datasize-tempsize)));
    peakpsnts(i)=peakpos+2;
    for j=1:tempsize;
        temp(j)=tr(peakpos+j);               
    end
    eval(['handles.template' num2str(i) '=zeros(tempsize,1);']);
    temp=(temp-min(temp))/max(temp-min(temp));            
    eval(['handles.template' num2str(i) '(1:tempsize)=temp;']);
    tr(peakpos+2)=0;
end

handles.numtemplates=peaknumber;
handles.tempcnt=1;    
numtemplates=handles.numtemplates;
finaltemplate=handles.template1;

for i=2:numtemplates;
    eval(['finaltemplate=finaltemplate+handles.template' num2str(i) ';']);
end


finaltemplate=finaltemplate/max(finaltemplate);
finaltemplate(1)=0;

axes(handles.axes4)
cla
plot(finaltemplate)
zoom on
axis auto

axes(handles.axes3)
cla
plot(handles.template1)
zoom on    
axis auto

handles.finaltemplate=finaltemplate;    

guidata(hObject, handles);



function numberofpeaks_Callback(hObject, eventdata, handles)
% hObject    handle to numberofpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numberofpeaks as text
%        str2double(get(hObject,'String')) returns contents of numberofpeaks as a double


% --- Executes during object creation, after setting all properties.
function numberofpeaks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numberofpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function templatesize_Callback(hObject, eventdata, handles)
% hObject    handle to templatesize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of templatesize as text
%        str2double(get(hObject,'String')) returns contents of templatesize as a double


% --- Executes during object creation, after setting all properties.
function templatesize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to templatesize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in shownextsubtemp.
function shownextsubtemp_Callback(hObject, eventdata, handles)
% hObject    handle to shownextsubtemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

templatecnt=handles.tempcnt;
numtemplates=handles.numtemplates;
templatecnt=mod(templatecnt,numtemplates)+1;
axes(handles.axes3)
cla
eval(['plot(handles.template' num2str(templatecnt) ')'])
set(handles.displayedtemp,'string',num2str(templatecnt));

handles.tempcnt=templatecnt;

guidata(hObject, handles);

% --- Executes on button press in excludesubtemplate.
function excludesubtemplate_Callback(hObject, eventdata, handles)
% hObject    handle to excludesubtemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

templatecnt=handles.tempcnt;
numtemplates=handles.numtemplates-1;

set(handles.displayedtemp,'string', num2str(templatecnt));

for i=templatecnt:numtemplates;
    eval(['handles.template' num2str(i) '=handles.template' num2str(i+1) ';']);
end

finaltemplate=handles.template1;

for i=2:numtemplates;
    eval(['finaltemplate=finaltemplate+handles.template' num2str(i) ';']);
end

finaltemplate=finaltemplate/max(finaltemplate);

finaltemplate(1)=0;
handles.finaltemplate=finaltemplate;

axes(handles.axes4)
cla
plot(finaltemplate)
zoom on

templatecnt=mod(templatecnt-1,numtemplates)+1;
eval(['currenttempl=handles.template' num2str(templatecnt) ';']);

axes(handles.axes3)
cla
plot(currenttempl)
zoom on

handles.templatecnt=templatecnt;
handles.numtemplates=numtemplates;
guidata(hObject, handles);

function displayedtemp_Callback(hObject, eventdata, handles)
% hObject    handle to displayedtemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displayedtemp as text
%        str2double(get(hObject,'String')) returns contents of displayedtemp as a double


% --- Executes during object creation, after setting all properties.
function displayedtemp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displayedtemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in templatemethod.
function templatemethod_Callback(hObject, eventdata, handles)
% hObject    handle to templatemethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of templatemethod

function [events, peakpos,cnt]=findpeaks(showregiontrace);

datasize=length(tr);
cnt=1;
i=1;

while i<=datasize;
        while ((i+1<=datasize)&&(tr(i+1)>=tr(i)));
            i=i+1;
        end
        peakpos(cnt)=i;
        cnt=cnt+1;
        while ((i+1<=datasize)&&(tr(i+1)<tr(i)));
            i=i+1;
        end

    i=i+1;
end

cnt=cnt-1;
events=zeros(datasize,1);
events(peakpos(:))=tr(peakpos(:));


% --- Executes on button press in trhist.
function trhist_Callback(hObject, eventdata, handles)
% hObject    handle to trhist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


numbins=str2num(get(handles.numbins,'string'));
usetemplate=get(handles.templatemethod,'value');

if usetemplate==0;
    dendnum=str2num(get(handles.dendnum,'string'));
    if handles.filtercounter(dendnum)==1;
        d=handles.tracematfltrd(dendnum,:);
    else
        d=handles.tracemat(dendnum,:);
    end
else
    d=handles.tr4events;
end

datasize=length(d);
maxdata=max(d);
mindata=min(d);

resolution=(maxdata-mindata)/(numbins-1);

trhist=zeros(numbins,2);
trhist(1:numbins,1)=resolution*(0:(numbins-1))+mindata;

for i=2:(datasize-1);
    
    if d(i)>=d(i-1) && d(i)>=d(i+1);
        pkval=d(i)-mindata;
        pkbin=floor(pkval/resolution)+1;
        trhist(pkbin,2)=trhist(pkbin,2)+1;
    end
end

figure
bar(trhist(:,1),trhist(:,2));

handles.tracehist=trhist;

guidata(hObject, handles);



function numbins_Callback(hObject, eventdata, handles)
% hObject    handle to numbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numbins as text
%        str2double(get(hObject,'String')) returns contents of numbins as a double


% --- Executes during object creation, after setting all properties.
function numbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in sorthor.
function sorthor_Callback(hObject, eventdata, handles)
% hObject    handle to sorthor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sorthor





function dendfrnum_Callback(hObject, eventdata, handles)
% hObject    handle to dendfrnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dendfrnum as text
%        str2double(get(hObject,'String')) returns contents of dendfrnum as a double


% --- Executes during object creation, after setting all properties.
function dendfrnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dendfrnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in showprevioussubtemplate.
function showprevioussubtemplate_Callback(hObject, eventdata, handles)
% hObject    handle to showprevioussubtemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

templatecnt=handles.tempcnt;
numtemplates=handles.numtemplates;
templatecnt=mod(templatecnt+numtemplates-2,numtemplates)+1;
axes(handles.axes3)
cla
eval(['plot(handles.template' num2str(templatecnt) ')'])
set(handles.displayedtemp,'string',num2str(templatecnt));

handles.tempcnt=templatecnt;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in cormat.
function cormat_Callback(hObject, eventdata, handles)
% hObject    handle to cormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


eventsmat=handles.eventsmat;

inds=find(eventsmat>0 & eventsmat <0);

dendpos=handles.posx;

if max(dendpos)>=64;
    pixels_x=128;
else
    pixels_x=64;
end

frames=size(eventsmat,2);
numdend=size(eventsmat,1);

if isempty(inds)==1;
    eventsmatcorrected=zeros(numdend,frames);
    for i=1:numdend;
        tr(1,:)=eventsmat(i,:);
        trc=zeros(1,frames+1);
        p=dendpos(i)/pixels_x;
        inds=find(tr>0);
        trc(inds+1)=p;
        trc(inds)=1-p;
        eventsmatcorrected(i,1:frames)=trc(1,2:frames+1);
    end
else
    eventsmatcorrected=eventsmat;
end
    
cormat=corrcoef(eventsmatcorrected');

corrmin=str2num(get(handles.corrmin,'string'));
corrmax=str2num(get(handles.corrmax,'string'));

handles.cormat=cormat;
axes(handles.axes6)
colormap(jet)
cla
imagesc(cormat, [corrmin corrmax]);

fovcoeff=str2num(get(handles.fov,'string'))/pixels_x;

corrvsdist=[];
for i=1:(numdend-1);
    for j=(i+1):numdend;
        dist=(dendpos(j)-dendpos(i))*fovcoeff;
        corval=cormat(i,j);
        corrvsdist=[corrvsdist; dist corval];
    end
end

handles.corrvsdist=corrvsdist;

guidata(hObject, handles); 




% --- Executes on button press in corvsdist.
function corvsdist_Callback(hObject, eventdata, handles)
% hObject    handle to corvsdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

corrvsdist=handles.corrvsdist;
axes(handles.axes3)
cla
scatter(corrvsdist(:,1), corrvsdist(:,2));
 


function corrmin_Callback(hObject, eventdata, handles)
% hObject    handle to corrmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmin as text
%        str2double(get(hObject,'String')) returns contents of corrmin as a double


% --- Executes during object creation, after setting all properties.
function corrmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmax_Callback(hObject, eventdata, handles)
% hObject    handle to corrmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmax as text
%        str2double(get(hObject,'String')) returns contents of corrmax as a double


% --- Executes during object creation, after setting all properties.
function corrmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savecorrs.
function savecorrs_Callback(hObject, eventdata, handles)
% hObject    handle to savecorrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cormat=handles.cormat;
corrvsdist=handles.corrvsdist;
   
filterspec = {'*.mat','mat files (.mat)'; '*.dat','dat files (.dat)'};
title = 'Save';

temp_file=handles.fname;
filename=handles.fname;
fn=handles.fn;

[filename, path, index] = uiputfile (filterspec, title, temp_file);


if index == 1;       
        eval(['m' fn 'cormat=cormat;']);
        save([path filename 'cormat.mat'],['m' fn 'cormat']);
        eval(['m' fn 'corrvsdist=corrvsdist;']);
        save([path filename 'corrvsdist.mat'],['m' fn 'corrvsdist']);
end


% --- Executes on button press in loadeventsmat.
function loadeventsmat_Callback(hObject, eventdata, handles)
% hObject    handle to loadeventsmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


filterspec = {'*.mat','MAT files (.mat)'};
title = 'load';
folder=handles.folder; 
[filename, path, index] = uigetfile (filterspec, title, folder);
eval(['addpath' ' ' char(39) path char(39) ';']);

folder=path;
total_filename=[folder filename];

if index == 1;
    fin=load(total_filename);
    names=fieldnames(fin);
    eventsmat=getfield(fin,names{1}) 
end

handles.eventsmat=eventsmat;
guidata(hObject, handles); 



% --- Executes on button press in loadpositions.
function loadpositions_Callback(hObject, eventdata, handles)
% hObject    handle to loadpositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filterspec = {'*.dat','dat files (.dat)';'*.xls','xls files (.xls)'};
title = 'load';
folder=handles.folder 
temp_file=folder;

[filename, path, index] = uigetfile (filterspec, title, temp_file);

eval(['addpath' ' ' char(39) path char(39) ';']);

folder=path;

total_filename=[folder filename];


   
if index == 1;
    eval(['a=open(' char(39) total_filename char(39) ');']);
    name=fieldnames(a);
    s=cell2struct(name(1,1),'trtemp',1);
    p=getfield(s,'trtemp');
    handles.posx=getfield(a,p);
end

guidata(hObject, handles);



function fov_Callback(hObject, eventdata, handles)
% hObject    handle to fov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fov as text
%        str2double(get(hObject,'String')) returns contents of fov as a double


% --- Executes during object creation, after setting all properties.
function fov_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


        



