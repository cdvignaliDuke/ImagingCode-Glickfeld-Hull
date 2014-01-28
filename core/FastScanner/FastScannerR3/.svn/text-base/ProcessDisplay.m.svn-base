function p = ProcessDisplay(object, event);

global fd;
persistent data; %TODO - check timing without it

p=0;

if fd.DAQ.framecounter==0 %allocate memory for variables used for test
%     data = zeros(fd.DAQ.InPointsToCollect*4,1,'uint16');
%     frameFinalData = zeros(fd.DAQ.InPointsToCollect,4,'uint16'); % vb: changed to uint16
%     frameFinalDataDouble = zeros(fd.DAQ.InPointsToCollect,4);
end

fd.DAQ.framecounter=fd.DAQ.framecounter+1;

% fd.tictoc(fd.DAQ.framecounter) = toc;
% tic;

fd.map = fopen(fd.DAQ.MemMapFile,'r');

if fd.map==-1
    set(fd.handles.lblStatus,'String',['Was not able to open  ' fd.DAQ.MemMapFile]);
    set(fd.handles.lblStatus,'ForegroundColor',[1 0 1]);
    return
end


datasize=fread(fd.map,1,'uint32');
if isempty(datasize) || datasize==0
    fclose(fd.map);
    return
end
shift=fread(fd.map,1,'uint32');
if isempty(shift)
    fclose(fd.map);
    return
end

dummydata = fread(fd.map,shift*4,'*uint16'); % now returns uint16

data = fread(fd.map,datasize*4,'*uint16'); % now returns uint16

fclose(fd.map);

set(fd.handles.lblPhase,'String',shift*1000/fd.DAQ.realinputRate);
set(fd.handles.lblPeriod,'String',datasize*1000/(120*fd.DAQ.realinputRate));
if length(data)~=datasize*4
    set(fd.handles.lblStatus,'String','Data size mismatch');
    set(fd.handles.lblStatus,'ForegroundColor',[1 0 1]);
    return
else
    set(fd.handles.lblStatus,'String','Running Display');
    set(fd.handles.lblStatus,'ForegroundColor',[1 0 0]);
end

frameFinalData= 4096 - reshape(data,4,datasize); % still uint16

step=2*pi*120/datasize;
ydd=(0+fd.timing.FastMirrorTriggerDelay*step+pi)*4000+1:step*4000:((datasize-1+fd.timing.FastMirrorTriggerDelay)*step+pi)*4000+1;
yda=0+fd.timing.FastMirrorTriggerDelay:datasize-1+fd.timing.FastMirrorTriggerDelay;

%%%cos_in=0:4000000;
%%%fd.cos_table=uint16((1-cos(cos_in/4000))*255/2+1);

xii=fd.cos_table(uint32(ydd));
yii=uint16((yda*239.999999/(datasize-1)+1)-0.5);


% fd.DAQ.ind=sub2ind_no_error_check([256,256],xii,yii);
fd.DAQ.ind=sub2ind_no_error_check([256,256],yii,xii); % indices are transposed
satind2=[];
if fd.img.Ch(2).use % process red channel first
    fd.img.Ch(2).CData1D(fd.DAQ.ind)=frameFinalData(2,:);
    fd.img.Ch(2).CData3D(:,:,fd.img.curInd+1) = reshape(fd.img.Ch(2).CData1D, 256,256);
    
    if fd.img.curInd == fd.img.AvgOnDisplay - 1
        fd.img.Ch(2).CData = mean(fd.img.Ch(2).CData3D,3);
        fd.dataShow(:,:,1)=min(fd.img.Ch(2).max,max(fd.img.Ch(2).min,fd.img.Ch(2).CData));    
        fd.dataShow(:,:,1)=(fd.dataShow(:,:,1)-fd.img.Ch(2).min)/(fd.img.Ch(2).max-fd.img.Ch(2).min);
        satind2=find(fd.img.Ch(2).CData>fd.img.Ch(2).max);
    end 
else
    fd.dataShow(:,:,1)=0;
end
satind1=[];
if fd.img.Ch(1).use % then process green channel
    fd.img.Ch(1).CData1D(fd.DAQ.ind)=frameFinalData(1,:);
    fd.img.Ch(1).CData3D(:,:,fd.img.curInd+1) = reshape(fd.img.Ch(1).CData1D, 256,256);
    
    if fd.img.curInd == fd.img.AvgOnDisplay - 1
        fd.img.Ch(1).CData = mean(fd.img.Ch(1).CData3D,3);
        fd.dataShow(:,:,2)=min(fd.img.Ch(1).max,max(fd.img.Ch(1).min,fd.img.Ch(1).CData));    
        fd.dataShow(:,:,2)=(fd.dataShow(:,:,2)-fd.img.Ch(1).min)/(fd.img.Ch(1).max-fd.img.Ch(1).min);
        satind1=find(fd.img.Ch(1).CData>fd.img.Ch(1).max);
%         set(fd.img.Ch(1).imH, 'CData', fd.dataShow);
    end
else
    fd.dataShow(:,:,2)=0;
end
fd.dataShow(:,:,3)=0;
if fd.img.Ch(1).use || fd.img.Ch(2).use
    if fd.img.curInd == fd.img.AvgOnDisplay - 1
        fd.dataShow(satind2)=1;
        fd.dataShow(satind2+256*256)=1;
        fd.dataShow(satind2+2*256*256)=1;
        fd.dataShow(satind1)=1;
        fd.dataShow(satind1+256*256)=1;
        fd.dataShow(satind1+2*256*256)=1;
        set(fd.img.Ch(1).imH, 'CData', fd.dataShow);
    end
end

drawnow;
yy=fd.stats.region.y1:fd.stats.region.y2;
xx=fd.stats.region.x1:fd.stats.region.x2;

if isfield(fd,'stats') && fd.stats.show
    if fd.img.Ch(1).use
        roi = fd.img.Ch(1).CData3D(yy,xx,:);
        fd.stats.Ch1 = [fd.stats.Ch1;prctile(roi(:),[5 50 95])];    
        fd.stats.Ch1(1,:)=[];
        set(fd.stats.plot1(1),'ydata',fd.stats.Ch1(:,1));
        set(fd.stats.plot1(2),'ydata',fd.stats.Ch1(:,2));
        set(fd.stats.plot1(3),'ydata',fd.stats.Ch1(:,3));
    end
    
    if fd.img.Ch(2).use
        roi = fd.img.Ch(2).CData3D(yy,xx,:);
        fd.stats.Ch2 = [fd.stats.Ch2;prctile(roi(:),[5 50 95])];    
        fd.stats.Ch2(1,:)=[];
        set(fd.stats.plot2(1),'ydata',fd.stats.Ch2(:,1));
        set(fd.stats.plot2(2),'ydata',fd.stats.Ch2(:,2));
        set(fd.stats.plot2(3),'ydata',fd.stats.Ch2(:,3));
    end
end

drawnow;

fd.img.curInd=mod(fd.img.curInd+1,fd.img.AvgOnDisplay);
