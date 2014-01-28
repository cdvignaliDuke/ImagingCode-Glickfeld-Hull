function p=ProcessStream(object, event)
% called every 5 ms;

global fs;

%% next part was to restart pCell waveform, commented by SY
%{
a=strfind(fs.DAQ.aoPockels.Running,'On');
if isempty(a)
    if mod(fs.DAQ.framecounter,75)==10
        restartPCellWaveform
    end
end
%}

%% timing test for each timer event

%      if fs.DAQ.TestMode==1 % if in test mode
%       fs.DAQ.clockcounter=fs.DAQ.clockcounter+1;
%           t=toc;
%           tic;
% %         fs.sAcquired(fs.DAQ.clockcounter)=fs.DAQ.ai.SamplesAcquired;
%          fs.sAvailable(fs.DAQ.clockcounter)=fs.DAQ.ai.SamplesAvailable;
%          fs.sAvailableFrames(fs.DAQ.clockcounter)=fs.DAQ.ai.SamplesAvailable/fs.DAQ.InPointsToCollect;
% %  %       fs.tictoc(fs.DAQ.clockcounter)=t;
% % %         if fs.DAQ.clockcounter==1
% % %             fs.totaltime(1)=t;
% % %         else
% % %             fs.totaltime(fs.DAQ.clockcounter)=fs.totaltime(fs.DAQ.clockcounter-1)+t;
% % %         end
%      end


%% wait for complete frame
% beware: fs.DAQ.InPointsToCollect varies from one call to another.
if fs.DAQ.ai.SamplesAvailable < fs.DAQ.InPointsToCollect
    return;
end


%% get frame
fs.DAQ.framecounter=fs.DAQ.framecounter+1;
%size(fs.dataTail)
frameFinalDataTransposed=[fs.dataTail (getdata(fs.DAQ.ai, fs.DAQ.InPointsToCollect, 'native') + 2048)'];
%frameFinalDataTransposed=(getdata(fs.DAQ.ai, fs.DAQ.InPointsToCollect, 'native') + 2048)';
%frameFinalDataTransposed=[fs.dataTail frameFinalDataTransposed];

%% timing test for each new frame
if fs.DAQ.TestMode==1 % if in test mode
    t=toc;
    tic;
    fs.sAcquired(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAcquired;
    fs.sAvailable(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAvailable;
    fs.sAvailableFrames(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAvailable/fs.DAQ.InPointsToCollect;
    fs.tictoc(fs.DAQ.framecounter)=t;
    if fs.DAQ.framecounter==1
        fs.totaltime(1)=t;
    else
        fs.totaltime(fs.DAQ.framecounter)=fs.totaltime(fs.DAQ.framecounter-1)+t;
    end
end


%% find period, phase and create indexes
% find zero crossings
%    xpos = int16(frameFinalDataTransposed(3,:))-2048;
%    sg = xpos>=0;
sg=frameFinalDataTransposed(3,:)>=2048; %looks faster than two lines above

dsg=diff(sg);
tcrossup=find(dsg>0); % sample before upward zero crossing

% calculate period
period=(tcrossup(end)-tcrossup(1))/(length(tcrossup)-1);
%fs.DAQ.period=period;

%calculate first min of x position (line start)
%do it only for first frame, after that data starts at minimum
if fs.DAQ.framecounter==1
    if tcrossup(1)>period/4
        linestart=round(tcrossup(1)-period/4);
    else
        linestart=round(tcrossup(1)+period*3/4);
        tcrossup(1)=[];
    end
    frameFinalDataTransposed(:,1:linestart)=[];
    
%     fs.greenStack = ij.ImageStack(256,240);
%     fs.redStack = ij.ImageStack(256,240);
%             
    fs.greenStack = ij.ImageStack(256,256);
    fs.redStack = ij.ImageStack(256,256);
end


if length(tcrossup)<128
    fs.dataTail=frameFinalDataTransposed;
    %collect little bit more to make sure that we have data for one frame
    fs.DAQ.InPointsToCollect=round(period*130 -size(fs.dataTail,2)); 
    return
end

% information - one frame is from tcrossup(1)-0.25*period (assume it is start of data) to tcrossup(128)+0.75*period
if size(frameFinalDataTransposed,2)<round(tcrossup(128)+0.75*period)
    fs.dataTail=frameFinalDataTransposed;
    %collect little bit more to make sure that we have data for one frame
    fs.DAQ.InPointsToCollect=round(tcrossup(128)+2.75*period -size(fs.dataTail,2));
    return
elseif size(frameFinalDataTransposed,2)>round(tcrossup(128)+0.75*period)
    fs.dataTail=frameFinalDataTransposed(:,round(tcrossup(128)+0.75*period)+1:end);
    %        frameFinalDataTransposed=frameFinalDataTransposed(:,1:round(tcrossup(128)+0.75*period));
    fs.DAQ.InPointsToCollect=round(period*130 -size(fs.dataTail,2));
else
    fs.dataTail=[];
    fs.DAQ.InPointsToCollect=round(period*130);
end

if fs.DAQ.TestMode==1 % if in test mode
    %for hdw trigger test - remove later
    fs.phase(fs.DAQ.framecounter)=tcrossup(1);
end




datasize=round(tcrossup(120)+0.75*period); %using only 240 lines
if  abs(datasize-fs.datasize)>20 %if big change in data size
    frameFinalDataTransposed=frameFinalDataTransposed(:,1:datasize);    
    fs.datasize=datasize;

    step=2*pi/period; %sample size in rad units
    %4000 is step of cos_table per one rad
    ydd=(fs.timing.FastMirrorTriggerDelay*step+pi)*4000+1:step*4000:((datasize-1+fs.timing.FastMirrorTriggerDelay)*step+pi)*4000+1;
    %        yda=0+fs.timing.FastMirrorTriggerDelay:datasize-1+fs.timing.FastMirrorTriggerDelay;
    ydaa=0.5+fs.timing.FastMirrorTriggerDelay*239.999999/(datasize-1):239.99999/(datasize-1):239.99999+0.5+fs.timing.FastMirrorTriggerDelay*239.999999/(datasize-1);
    yii=uint16(ydaa);
    xii=fs.cos_table(uint32(ydd));
    %%       yii=uint16((yda*239.999999/(datasize-1)+1)-0.5);
    fs.DAQ.ind=sub2ind_no_error_check([256,256],yii,xii); % indices are transposed    
else
    frameFinalDataTransposed=frameFinalDataTransposed(:,1:fs.datasize);
    %data set has the same size as previous; no need to create new indexes
end


% %image reindexing
% datasize=size(frameFinalDataTransposed,2);
% 
% 
% 
% if datasize~=fs.datasize
%     fs.datasize=datasize;
%     fs.timing.FastMirrorTriggerDelay=0; %just for test
%     step=2*pi/period;
%     ydd=(0+fs.timing.FastMirrorTriggerDelay*step+pi)*4000+1:step*4000:((datasize-1+fs.timing.FastMirrorTriggerDelay)*step+pi)*4000+1;
%     %        yda=0+fs.timing.FastMirrorTriggerDelay:datasize-1+fs.timing.FastMirrorTriggerDelay;
%     ydaa=0.5+fs.timing.FastMirrorTriggerDelay:239.99999/(datasize-1):239.99999+0.5+fs.timing.FastMirrorTriggerDelay;
%     yii=uint16(ydaa);
%     xii=fs.cos_table(uint32(ydd));
%     %%       yii=uint16((yda*239.999999/(datasize-1)+1)-0.5);
%     fs.DAQ.ind=sub2ind_no_error_check([256,256],yii,xii); % indices are transposed
% end

%% create image matrixes
fs.DAQ.ch(1).CData1D(fs.DAQ.ind)=frameFinalDataTransposed(1,:)';
fs.DAQ.ch(2).CData1D(fs.DAQ.ind)=frameFinalDataTransposed(2,:)';

% fs.DAQ.ch(2).CData1D(fs.DAQ.ind)=frameFinalDataTransposed(3,:)';
% fs.DAQ.ch(1).CData1D(fs.DAQ.ind)=frameFinalDataTransposed(4,:)';

frame1=reshape(fs.DAQ.ch(1).CData1D,256,256);
frame2=reshape(fs.DAQ.ch(2).CData1D,256,256);

%frame1=reshape(fs.DAQ.ch(1).CData1D(1:240*256),240,256);
%frame2=reshape(fs.DAQ.ch(2).CData1D(1:240*256),240,256);


%frame1=frame1(1:240,:);
%frame2=frame2(1:240,:);

%% set size of next block of data
%     if fs.DAQ.framecounter==1
%         fs.DAQ.InPointsToCollect=round(linestart+period*128*2-100-fs.DAQ.InPointsToCollect);
%     else
%         fs.DAQ.InPointsToCollect=round(period*128+(tcrossup(1)-400)/4);
%   %      fs.DAQ.InPointsToCollect=round(period*128 +period*3/4-(length(xpos)-tcrossup(end)));
%     end
%

%% prepare and send data to display if have time only
if fs.DAQ.ai.SamplesAvailable < fs.DAQ.InPointsToCollect*100
    % write frame to ram disk
    frewind(fs.map);
    fwrite(fs.map,fs.DAQ.ch(1).CData1D,'uint16');
    fwrite(fs.map,fs.DAQ.ch(2).CData1D,'uint16');

end

%% write data to hard drive
if fs.DAQ.DoNotSave==0
%     fn1 = sprintf('%s%06d.tif',[fs.DAQ.BaseFileNameLong '_1_'],fs.DAQ.framecounter);
%     fn2 = sprintf('%s%06d.tif',[fs.DAQ.BaseFileNameLong '_2_'],fs.DAQ.framecounter);
%     overwrite=1;
%     wtifc(uint16(frame1), [], fn1, 'none', '', 72, overwrite, 'rgb');
%     wtifc(uint16(frame2), [], fn2, 'none', '', 72, overwrite, 'rgb');
%     set(fs.handles.lblCurrentFileName,'String',fn1);
%     
    
%    greenProcessor = ij.process.ByteProcessor(256,240);
    greenProcessor = ij.process.ByteProcessor(256,256);
    greenProcessor.setFloatArray(frame1'./16)
    fs.greenStack.addSlice('frame', greenProcessor);
    
%    redProcessor = ij.process.ByteProcessor(256,240);
    redProcessor = ij.process.ByteProcessor(256,256);
    
    redProcessor.setFloatArray(frame2'./16)
    fs.redStack.addSlice('frame', redProcessor);
    
    
    if mod(fs.DAQ.framecounter,fs.tifStackSize) == 0
        fs.DAQ.filecounter=fs.DAQ.filecounter+1;
        greenPlus = ij.ImagePlus('GreenStack',fs.greenStack);
        greenFileSaver = ij.io.FileSaver(greenPlus);
        fn1 = sprintf('%s%06d.tif',[fs.DAQ.BaseFileNameLong '_1_'],fs.DAQ.filecounter);   
        set(fs.handles.lblCurrentFileName,'String',fn1);
        greenFileSaver.saveAsTiffStack(fn1);
%        fs.greenStack = ij.ImageStack(256,240);
        fs.greenStack = ij.ImageStack(256,256);
        
    
        redPlus = ij.ImagePlus('GreenStack',fs.redStack);
        redFileSaver = ij.io.FileSaver(redPlus);
        fn2 = sprintf('%s%06d.tif',[fs.DAQ.BaseFileNameLong '_2_'],fs.DAQ.filecounter);        
        redFileSaver.saveAsTiffStack(fn2);
%        fs.redStack = ij.ImageStack(256,240);
        fs.redStack = ij.ImageStack(256,256);
        
     end
end

%% display information on GUI
[hour, minute, second] = sec2hms(etime(clock,fs.timing.startime));
str = sprintf('%s %02d:%02d:%02d %07d %07d',fs.DAQ.modeString,hour,minute,round(second),fs.DAQ.framecounter,fs.DAQ.pCellCounter-fs.DAQ.framecounter);
set(fs.handles.lblProgress,'string',str);
drawnow

%% restart if acquired n frames per step
if fs.cycle.Do
    if fs.cycle.NFramesPrStep>0
        if mod(fs.DAQ.framecounter,fs.cycle.NFramesPrStep)==0
            StopDAQ;
            if fs.DAQ.framecounter<fs.cycle.NFramesPrStep*fs.cycle.NSteps*fs.cycle.NReps
                %change parameters inside loop here - TODO
                loopN=ceil(fs.DAQ.framecounter/(fs.cycle.NFramesPrStep*fs.cycle.NSteps));
                stepN=ceil((fs.DAQ.framecounter-(loopN-1)*fs.cycle.NFramesPrStep*fs.cycle.NSteps)/fs.cycle.NFramesPrStep);
                if fs.stage.ready
                    if fs.cycle.ZStep~=0
                        pos=SetPositionZ(fs.cycle.ZStart+fs.cycle.ZStep*stepN);
                        fs.position.z=pos;
                        set(fs.handles.txtZPos,'String',num2str(fs.position.z));
                        set(fs.handles.lblCurrentZ,'String',num2str(fs.position.z));
                        if fs.position.z<get(fs.handles.sliderZPos,'Min')
                            set(fs.handles.sliderZPos,'Min',fs.position.z);
                        end
                        if fs.position.z>get(fs.handles.sliderZPos,'Max')
                            set(fs.handles.sliderZPos,'Max',fs.position.z);
                        end
                        set(fs.handles.sliderZPos,'Value',fs.position.z);
                    end
                end
                %restart if not all frames are collected
                %start flushes DAQ engine
                StartDAQ;
            else
                %stop aquisition
                fs.DAQ.Running=0;
            end
        end
    end
end

%% check whether stimulation done
if fs.DAQ.stopLineEnabled==1 && getvalue(fs.DAQ.stopLine)==1
    fs.DAQ.Running = 0;
    disp('Acquisition interrupted by visual stimulator');
end

%% return (to timer) or stop acquisition
if fs.DAQ.Running
    return;
end
if fs.DAQ.DoNotSave==0
%     fn1 = sprintf('%s%06d.tif',[fs.DAQ.BaseFileNameLong '_1_'],fs.DAQ.framecounter);
%     fn2 = sprintf('%s%06d.tif',[fs.DAQ.BaseFileNameLong '_2_'],fs.DAQ.framecounter);
%     overwrite=1;
%     wtifc(uint16(frame1), [], fn1, 'none', '', 72, overwrite, 'rgb');
%     wtifc(uint16(frame2), [], fn2, 'none', '', 72, overwrite, 'rgb');
%     set(fs.handles.lblCurrentFileName,'String',fn1);
%     
    
%    greenProcessor = ij.process.ByteProcessor(256,240);
    greenProcessor = ij.process.ByteProcessor(256,256);

    greenProcessor.setFloatArray(frame1')
    fs.greenStack.addSlice('frame', greenProcessor);
    
%    redProcessor = ij.process.ByteProcessor(256,240);
    redProcessor = ij.process.ByteProcessor(256,256);
    
    redProcessor.setFloatArray(frame2')
    fs.redStack.addSlice('frame', greenProcessor);
    
        fs.DAQ.filecounter=fs.DAQ.filecounter+1;
        greenPlus = ij.ImagePlus('GreenStack',fs.greenStack);
        greenFileSaver = ij.io.FileSaver(greenPlus);
        fn1 = sprintf('%s%06d.tif',[fs.DAQ.BaseFileNameLong '_1_'],fs.DAQ.framecounter);        
        greenFileSaver.saveAsTiffStack(fn1);

    
        redPlus = ij.ImagePlus('GreenStack',fs.redStack);
        redFileSaver = ij.io.FileSaver(redPlus);
        fn2 = sprintf('%s%06d.tif',[fs.DAQ.BaseFileNameLong '_2_'],fs.DAQ.framecounter);        
        redFileSaver.saveAsTiffStack(fn2);

    
end
%% stop acquisition
StopAcquisition;


%% plot debug information if in debug mode
if fs.DAQ.TestMode==1
    %    figure;
    %    plot(fs.tictoc(:));
    figure;
    plot(fs.sAvailable(:),'r');
    %    plot(fs.phase(:),'r');
    figure;
    plot(fs.sAvailableFrames(:),'g');

end

