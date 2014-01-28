function DoDisplay

%% ini parameters
global fd;


%%
fd.ParamMap = fopen(fd.DAQ.ParamMap,'r');
if fd.ParamMap==-1
    set(fd.handles.lblStatus,'String',['Was not able to open  ' fd.DAQ.ParamMap]);
    set(fd.handles.lblStatus,'ForegroundColor',[1 0 1]);
    set(fd.handles.lblSamplingRate,'ForegroundColor',[1 0 0]);
    fd.DAQ.realinputRate=fd.DAQ.inputRate;
%    h = warndlg(['Input rate set to ' num2str(fd.DAQ.realinputRate)],'Error reading map file','modal');
else
    fd.DAQ.realinputRate = fread(fd.ParamMap,1,'double');
    fclose(fd.ParamMap);
end
set(fd.handles.lblSamplingRate,'String',fd.DAQ.realinputRate);

if fd.img.Ch(1).use
    set(fd.img.Ch(1).ax,'XLim',[1 256]);
    set(fd.img.Ch(1).ax,'YLim',[1 256]);
    set(fd.img.Ch(1).imH, 'EraseMode', 'none', 'CDataMapping','scaled');    
end
if fd.img.Ch(2).use
    set(fd.img.Ch(2).ax,'XLim',[1 256]);
    set(fd.img.Ch(2).ax,'YLim',[1 256]);
    set(fd.img.Ch(2).imH, 'EraseMode', 'none', 'CDataMapping','scaled');
end




fd.img.curInd=0;
fd.img.Ch(1).CData =zeros(256,256);
fd.img.Ch(2).CData =zeros(256,256);

 fd.img.Ch(1).CDataV2 =zeros(256*256,1);
 fd.img.Ch(2).CDataV2 =zeros(256*256,1);

fd.dataShow=zeros(256,256,3);

fd.img.Ch(1).CData1D =zeros(256*256,1);
fd.img.Ch(2).CData1D =zeros(256*256,1);
fd.img.Ch(1).CData3D =zeros(256,256,fd.img.AvgOnDisplay);
fd.img.Ch(2).CData3D =zeros(256,256,fd.img.AvgOnDisplay);

fd.img.Ch(1).CData2D =zeros(256*256,fd.img.AvgOnDisplay);
fd.img.Ch(2).CData2D =zeros(256*256,fd.img.AvgOnDisplay);


% % delete any left-over timer
% timers = timerfind;
% for itimer = 1:length(timers);
%     if strcmp(func2str(timers(itimer).TimerFcn),'ProcessDisplay');
%         stop(timers(itimer));
%         delete(timers(itimer));
%     end
% end
    
%fd.DisplayTimer = timer('TimerFcn',@ProcessDisplay, 'Period', 0.062,'ExecutionMode','FixedRate');
fd.DisplayTimer = timer('TimerFcn',@ProcessDisplay, 'Period', 0.005,'ExecutionMode','FixedRate');

fd.DAQ.framecounter = 0;
fd.DAQ.datasize=0;
fd.DAQ.framecounterInScanner = 0;
%create cos LUT
cos_in=0:4000000;
fd.cos_table=uint16((1-cos(cos_in/4000))*255/2+1);

tic;

start(fd.DisplayTimer);

