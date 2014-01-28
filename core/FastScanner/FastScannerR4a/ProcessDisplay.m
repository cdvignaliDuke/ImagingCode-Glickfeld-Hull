function p = ProcessDisplay(object, event);

global fd;
%persistent data; %TODO - check timing without it

p=0;

if fd.DAQ.framecounter==0 %allocate memory for variables used for test
    %     data = zeros(fd.DAQ.InPointsToCollect*4,1,'uint16');
    %     frameFinalData = zeros(fd.DAQ.InPointsToCollect,4,'uint16'); % vb: changed to uint16
    %     frameFinalDataDouble = zeros(fd.DAQ.InPointsToCollect,4);
end

fd.map = fopen(fd.DAQ.MemMapFile,'r');

if fd.map==-1
    set(fd.handles.lblStatus,'String',['Was not able to open  ' fd.DAQ.MemMapFile]);
    set(fd.handles.lblStatus,'ForegroundColor',[1 0 1]);
    %    pause(0.01)
    return
end

frameNumberInScanner=fread(fd.map,1,'uint32');

if isempty(frameNumberInScanner) || fd.DAQ.framecounterInScanner==frameNumberInScanner
    fclose(fd.map);
    return
end
fd.DAQ.framecounter=fd.DAQ.framecounter+1;

fd.tictoc(fd.DAQ.framecounter) = toc;
tic;


fd.DAQ.framecounterInScanner=frameNumberInScanner;

datasize=fread(fd.map,1,'uint32');
if isempty(datasize) || datasize==0
    fclose(fd.map);
    return
end
data_green = fread(fd.map,datasize,'*uint16'); % now returns uint16
data_red = fread(fd.map,datasize,'*uint16'); % now returns uint16

fclose(fd.map);

%set(fd.handles.lblPeriod,'String',datasize);



satind2=[];
if fd.img.Ch(2).use % process red channel first
    fd.img.Ch(2).CData2D(1:length(data_red),fd.img.curInd+1) = data_red;
    fd.img.Ch(2).CDataV2 = mean(fd.img.Ch(2).CData2D,2);
    tmpData=min(fd.img.Ch(2).max,max(fd.img.Ch(2).min,fd.img.Ch(2).CDataV2));
    fd.dataShow(1:length(tmpData))=(tmpData-fd.img.Ch(2).min)/(fd.img.Ch(2).max-fd.img.Ch(2).min);
    satind2=find(fd.img.Ch(2).CDataV2>fd.img.Ch(2).max | fd.img.Ch(2).CDataV2==4096);
else
    fd.dataShow(:,:,1)=0;
end
satind1=[];
if fd.img.Ch(1).use % then process green channel
    fd.img.Ch(1).CData2D(1:length(data_green),fd.img.curInd+1) = data_green;
    fd.img.Ch(1).CDataV2 = mean(fd.img.Ch(1).CData2D,2);
    tmpData=min(fd.img.Ch(1).max,max(fd.img.Ch(1).min,fd.img.Ch(1).CDataV2));
    fd.dataShow(256*256+1:256*256*2)=(tmpData-fd.img.Ch(1).min)/(fd.img.Ch(1).max-fd.img.Ch(1).min);
    satind1=find(fd.img.Ch(1).CDataV2>fd.img.Ch(1).max | fd.img.Ch(1).CDataV2==4096);
else
    fd.dataShow(:,:,2)=0;
end


fd.dataShow(:,:,3)=0;

if fd.img.Ch(1).use || fd.img.Ch(2).use
    %    if fd.img.curInd == fd.img.AvgOnDisplay - 1
    fd.dataShow(satind2)=1;
    fd.dataShow(satind2+256*256)=1;
    fd.dataShow(satind2+2*256*256)=1;
    fd.dataShow(satind1)=1;
    fd.dataShow(satind1+256*256)=1;
    fd.dataShow(satind1+2*256*256)=1;

    fd.dataShow(:,:,1)=rot90(fd.dataShow(:,:,1),3);
    fd.dataShow(:,:,2)=rot90(fd.dataShow(:,:,2),3);
    fd.dataShow(:,:,3)=rot90(fd.dataShow(:,:,3),3);

    set(fd.img.Ch(1).imH, 'CData', fd.dataShow);
    drawnow;
    %    end
end


yy=fd.stats.region.y1:fd.stats.region.y2;
xx=fd.stats.region.x1:fd.stats.region.x2;

if isfield(fd,'stats') && fd.stats.show
    if fd.img.Ch(1).use
        %        roi = fd.img.Ch(1).CData3D(yy,xx,:);
        fd.img.Ch(1).CData = reshape(fd.img.Ch(1).CDataV2, 256,256);
        roi = fd.img.Ch(1).CData(yy,xx);
        nel=numel(roi);
        if nel>0
            srtRoi=sort(roi(:));
            st05=srtRoi(round((nel-1)*0.05)+1);
            st50=srtRoi(round((nel-1)*0.5)+1);
            st95=srtRoi(round((nel-1)*0.95)+1);
            fd.stats.Ch1 = [fd.stats.Ch1;[st05 st50 st95]];
        end
        %        fd.stats.Ch1 = [fd.stats.Ch1;prctile(roi(:),[5 50 95])];
        fd.stats.Ch1(1,:)=[];
        set(fd.stats.plot1(1),'ydata',fd.stats.Ch1(:,1));
        set(fd.stats.plot1(2),'ydata',fd.stats.Ch1(:,2));
        set(fd.stats.plot1(3),'ydata',fd.stats.Ch1(:,3));
    end

    if fd.img.Ch(2).use
        fd.img.Ch(2).CData = reshape(fd.img.Ch(2).CDataV2, 256,256);
        roi = fd.img.Ch(2).CData(yy,xx);
        nel=numel(roi);
        if nel>0
            srtRoi=sort(roi(:));
            st05=srtRoi(round((nel-1)*0.05)+1);
            st50=srtRoi(round((nel-1)*0.5)+1);
            st95=srtRoi(round((nel-1)*0.95)+1);
            fd.stats.Ch2 = [fd.stats.Ch2;[st05 st50 st95]];
        end
        %       fd.stats.Ch2 = [fd.stats.Ch2;prctile(roi(:),[5 50 95])];
        fd.stats.Ch2(1,:)=[];
        set(fd.stats.plot2(1),'ydata',fd.stats.Ch2(:,1));
        set(fd.stats.plot2(2),'ydata',fd.stats.Ch2(:,2));
        set(fd.stats.plot2(3),'ydata',fd.stats.Ch2(:,3));
    end
end

drawnow;

fd.img.curInd=mod(fd.img.curInd+1,fd.img.AvgOnDisplay);
