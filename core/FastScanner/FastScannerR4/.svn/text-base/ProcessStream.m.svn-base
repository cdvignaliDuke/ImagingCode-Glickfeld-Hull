function p=ProcessStream(object, event)
% called every 5 ms;

dad_fclose_delay = 55;
dad_fclose_delay_ms = dad_fclose_delay*31;

global fs;

currT = toc;
%tic
%% check if waiting to start acq (as in a new z level)
if fs.DAQ.zIsWaiting
    %fprintf(1, 'waiting: diff %6.3f\n', fs.DAQ.zNextStartTime-currT);
    assert(fs.DAQ.started == false, 'bug: waiting requested but still running');
    if currT < fs.DAQ.zNextStartTime; %make sure that we had enouth time to finish previous file
        return
    else
        % close prev, start acq again, and start a new file
        fs.DAQ.zIsWaiting = false;

%         if ~isempty(fs.fid_prev)  % empty if the last file has < dad_fclose_delay frames; is closed in process frame loop below
%             fclose(fs.fid_prev);
%             fs.fid_prev = [];
%         end
        fs.DAQ.InPointsToCollect=round(1.5*(256)*(fs.DAQ.realinputRate*fs.timing.FastMirrorPeriod/2000));

        StartDAQ;  % initializes framcounter_thisacq to 0
        return % next iteration into this fcn will handle data
    end
end

%% wait for complete frame
sampAv=fs.DAQ.ai.SamplesAvailable;
% beware: fs.DAQ.InPointsToCollect varies from one call to another.
if sampAv < fs.DAQ.InPointsToCollect
%    toc
    return;
end

%% process frame
fs.DAQ.framecounter=fs.DAQ.framecounter+1;
fs.DAQ.framecounter_thisacq = fs.DAQ.framecounter_thisacq + 1;

%collect little bit more than one frame

frameFinalDataTransposed=(getdata(fs.DAQ.ai, fs.DAQ.InPointsToCollect, 'native') + 2048)';

if fs.DAQ.framecounter_thisacq==1 %check trigger position using y signal
    yPosMax=max(frameFinalDataTransposed(4,:));
    yPosMin=min(frameFinalDataTransposed(4,:));
    yPosTshld=yPosMin+(yPosMax-yPosMin)/18;
    ycross=find(frameFinalDataTransposed(4,:)<yPosTshld);
    if ycross(1)>1
        frameFinalDataTransposed=frameFinalDataTransposed(:,ycross(1)+100:end);
    end
    
else

  frameFinalDataTransposed=[fs.DAQ.data_leftover frameFinalDataTransposed];

end

if fs.DAQ.TestMode==1 % if in test mode
    fs.sAcquired(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAcquired;
    fs.sAvailable(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAvailable;
    fs.sAvailableFrames(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAvailable/fs.DAQ.InPointsToCollect;
end
%good=0;
if fs.DAQ.framecounter_thisacq==1

    tcrossup=find(frameFinalDataTransposed(3,1:end-6)<=2048 & frameFinalDataTransposed(3,2:end-5)>2048 & frameFinalDataTransposed(3,3:end-4)>2048 & frameFinalDataTransposed(3,4:end-3)>2048 & frameFinalDataTransposed(3,5:end-2)>2048 & frameFinalDataTransposed(3,6:end-1)>2048 & frameFinalDataTransposed(3,7:end)>2048);
    firstpoint=min(tcrossup(1)-1,20);
    bbb=regress(double(frameFinalDataTransposed(3,tcrossup(1)-firstpoint:tcrossup(1)+20)-2048)',[ones(firstpoint+21,1) (-firstpoint:20)']);
    tcrossup(1)=tcrossup(1)-bbb(1)/bbb(2);

    dcrossup=diff(tcrossup);
    smdc=length(find(dcrossup<800*fs.DAQ.realinputRate/3300000)); %count false crossings
    smdc2=length(find(dcrossup<20*fs.DAQ.realinputRate/3300000));
%    period=(tcrossup(end)-tcrossup(1))/(length(tcrossup)-1-smdc2-(smdc-smdc2)/2);
    period=845*fs.DAQ.realinputRate/3300000;
    
    fs.DAQ.lastperiod=period;
    firstpoint=20;
    X=[ones(firstpoint+21,1) (-firstpoint:20)'];
    [Q,R,perm]=qr(X,0);
    bbb=zeros(2,1);

for line=1:122
    approx_next_tcrossup=round(tcrossup(line)+fs.DAQ.lastperiod);
    y=double(frameFinalDataTransposed(3,approx_next_tcrossup-20:approx_next_tcrossup+20)-2048)';
    bbb(perm)=R\(Q'*y);
    
    tcrossup(line+1)=approx_next_tcrossup-bbb(1)/bbb(2);
    fs.DAQ.lastperiod=tcrossup(line+1)-tcrossup(line);
end  
    period=(tcrossup(101)-tcrossup(1))/100;    

%{    
    da=double(frameFinalDataTransposed(3,1:round(period)));
    fda=fft(da);
    n=10;
    fda(n+1:end-n)=0;
    nda=real(ifft(fda));
    tcrossup_nda=find(nda(1:end-6)<=2048 & nda(2:end-5)>2048 & nda(3:end-4)>2048 & nda(4:end-3)>2048 & nda(5:end-2)>2048 & nda(6:end-1)>2048 & nda(7:end)>2048);
    tcrossup(1)=tcrossup_nda(1);
%}    
%    pos101=round(tcrossup(1)+100*period);
%    bbb=regress(double(frameFinalDataTransposed(3,pos101-20:pos101+20)-2048)',[ones(41,1) (-20:20)']);
%    tcrossup(101)=pos101-bbb(1)/bbb(2);
%    period=(tcrossup(101)-tcrossup(1))/100;
else

    le=round(min(1000*fs.DAQ.realinputRate/3300000,size(frameFinalDataTransposed,2)));
    tcrossup=find(frameFinalDataTransposed(3,1:le-6)<=2048 & frameFinalDataTransposed(3,2:le-5)>2048 & frameFinalDataTransposed(3,3:le-4)>2048 & frameFinalDataTransposed(3,4:le-3)>2048 & frameFinalDataTransposed(3,5:le-2)>2048 & frameFinalDataTransposed(3,6:le-1)>2048 & frameFinalDataTransposed(3,7:le)>2048);
    firstpoint=min(tcrossup(1)-1,20);
    [bbb,bint,r]=regress(double(frameFinalDataTransposed(3,tcrossup(1)-firstpoint:tcrossup(1)+20)-2048)',[ones(firstpoint+21,1) (-firstpoint:20)']);

    
    %{    
%    good=sum(r.^2);
%     if good>400
%         fs.DAQ.calculating_delay=1;
%     end
    
    if good>300
    da=double(frameFinalDataTransposed(3,1:round(fs.DAQ.lastperiod)*8));
    fda=fft(da);
    n=8;
    fda(n+1:end-n)=0;
    nda=real(ifft(fda));
    tcrossup_nda=find(nda(1:end-6)<=2048 & nda(2:end-5)>2048 & nda(3:end-4)>2048 & nda(4:end-3)>2048 & nda(5:end-2)>2048 & nda(6:end-1)>2048 & nda(7:end)>2048);
    tcrossup(1)=tcrossup_nda(1);    
%    else
%}    
 %{   
     if sum(r.^2)>1000
         figure(200)
         hold off
         plot(frameFinalDataTransposed(3,tcrossup(1)-firstpoint:tcrossup(1)+20)-2048,'r');
         hold on
         plot(nda(tcrossup_nda(1)-firstpoint:tcrossup_nda(1)+20)-2048,'b');
         
     else
%         figure(200)
%         plot(frameFinalDataTransposed(3,tcrossup(1)-firstpoint:tcrossup(1)+20)-2048,'b');
%         
     end
    %}
     tcrossup(1)=tcrossup(1)-bbb(1)/bbb(2);
%     end

%{
    pos101=round(tcrossup(1)+100*fs.DAQ.lastperiod-fs.DAQ.lastperiod/2);    
    da=double(frameFinalDataTransposed(3,pos101:pos101+round(fs.DAQ.lastperiod)));
    fda=fft(da);
    n=6;
    fda(n+1:end-n)=0;
    nda=real(ifft(fda));
    tcrossup_nda=find(nda(1:end-6)<=2048 & nda(2:end-5)>2048 & nda(3:end-4)>2048 & nda(4:end-3)>2048 & nda(5:end-2)>2048 & nda(6:end-1)>2048 & nda(7:end)>2048);
    tcrossup(101)=tcrossup_nda(1)+pos101-1;    
%} 
%tic
    firstpoint=20;
    X=[ones(firstpoint+21,1) (-firstpoint:20)'];
    [Q,R,perm]=qr(X,0);
    bbb=zeros(2,1);

for line=1:122
    approx_next_tcrossup=round(tcrossup(line)+fs.DAQ.lastperiod);
%    firstpoint=min(tcrossup(line)-1,20);
    y=double(frameFinalDataTransposed(3,approx_next_tcrossup-20:approx_next_tcrossup+20)-2048)';

    bbb(perm)=R\(Q'*y);

    tcrossup(line+1)=approx_next_tcrossup-bbb(1)/bbb(2);
    fs.DAQ.lastperiod=tcrossup(line+1)-tcrossup(line);
end  
%toc
 %   bbb=regress(double(frameFinalDataTransposed(3,pos101-20:pos101+20)-2048)',[ones(41,1) (-20:20)']);
 %   tcrossup(101)=pos101-bbb(1)/bbb(2);

    period=(tcrossup(101)-tcrossup(1))/100;    
%    pos100=firsttcrossup(1)+round(100*fs.DAQ.lastperiod);
%    lasttcrossup=find(frameFinalDataTransposed(3,pos100-100:pos100+100-6)<=2048 & frameFinalDataTransposed(3,pos100-100+1:pos100+100-5)>2048 & frameFinalDataTransposed(3,pos100-100+6:pos100+100)>2048);
%    tcrossup(1)=firsttcrossup(1);
%    period=(lasttcrossup(1)+pos100-100-1-firsttcrossup(1))/100;
end

fs.DAQ.lastperiod=period;
%firstpoint=min(tcrossup(1)-1,20);
%bbb=regress(double(frameFinalDataTransposed(3,tcrossup(1)-firstpoint:tcrossup(1)+20)-2048)',[ones(firstpoint+21,1) (-firstpoint:20)']);
%tcrossup(1)=tcrossup(1)-bbb(1)/bbb(2);
%calculate first min of x position (line start)
maxdelay=50;
maxdelay=100;
if tcrossup(1)>period/4+maxdelay+1
    linestart=round(tcrossup(1)-period/4);
else
    linestart=round(tcrossup(1)+period*3/4);
end

if fs.DAQ.framecounter==1
    fs.DAQ.linestart=linestart;
    fs.DAQ.period=period;
end

if fs.DAQ.TestMode==1
    %for hdw trigger test - remove later
    fs.phase(fs.DAQ.framecounter)=tcrossup(1);
end

%% reconstruct frame

%calculate all line starts
cyclelen=uint32(round(period/2)*2);

all_linestarts=uint32(round(linestart+(0:128)*period-maxdelay));
%if fs.DAQ.framecounter_thisacq>1
    all_linestarts(1:123)=uint32(round(tcrossup(1:123)-tcrossup(1)+linestart-maxdelay));
%end

%framestart=all_linestarts(1);
framestop=all_linestarts(122);

if fs.DAQ.calculating_delay
    fs.DAQ.calculating_delay=0;
    if ~isfield(fs.DAQ,'delay')
        fs.DAQ.delay=NaN;
    end
    [delay,err] = calculatedelay(sum(frameFinalDataTransposed(1:2,1:framestop),'native'),all_linestarts(1:120),cyclelen,fs.DAQ.delay,10);
    
    fs.DAQ.delay=delay;    
    set(fs.handles.txtDelay,'String',fs.DAQ.delay);    
end

green_frame=formatframe(frameFinalDataTransposed(1,1:framestop),all_linestarts(1:120),cyclelen,uint32(fs.DAQ.delay));
red_frame=formatframe(frameFinalDataTransposed(2,1:framestop),all_linestarts(1:120),cyclelen,uint32(fs.DAQ.delay));


%{
%test mode
if fs.DAQ.calculating_delay
    fs.DAQ.calculating_delay=0;
    if isfield(fs.DAQ,'delay')
        [delay,err] = CalculateDelay(frameFinalDataTransposed(3,framestart:framestop),all_linestarts(1:120),cyclelen,fs.DAQ.delay,5);
        fs.DAQ.delay=delay;
    else
        fs.DAQ.delay=NaN;
        [delay,err] = CalculateDelay(frameFinalDataTransposed(3,framestart:framestop),all_linestarts(1:120),cyclelen,fs.DAQ.delay,5);
        fs.DAQ.delay=delay;
    end
    set(fs.handles.txtDelay,'String',fs.DAQ.delay);
end


green_frame=formatframe(frameFinalDataTransposed(3,framestart:framestop),all_linestarts(1:120),cyclelen,uint32(fs.DAQ.delay));
red_frame=formatframe(frameFinalDataTransposed(4,framestart:framestop),all_linestarts(1:120),cyclelen,uint32(fs.DAQ.delay));
%}    


%calculate index LUT

%format frame
    step=2*pi/period;
%    step=pi/period;
    
    yda = 0:round(period/2)-1;
%    yda = 0:round(period)-1;
    
    xii=uint32((1-cos(yda*step))*(256-1)/2+1);
%    xii=uint32((yda/(round(period/2)-1))*(256-1)+1);

 fs.DAQ.fgreen(:)=2048-binpixels(green_frame,uint16(xii),uint16(256))';
 fs.DAQ.fred(:)=2048-binpixels(red_frame,uint16(xii),uint16(256))';

%fs.DAQ.fgreen(:)=2048-binpixels(green_frame,uint16(xii),uint16(256));
%fs.DAQ.fred(:)=2048-binpixels(red_frame,uint16(xii),uint16(256));


%% write frame to ram disk
%if good<300
frewind(fs.map);
fwrite(fs.map,fs.DAQ.framecounter,'uint32'); %current frame number
fwrite(fs.map,length(fs.DAQ.fgreen(:)),'uint32'); %image size in DAQ points
fwrite(fs.map,fs.DAQ.fgreen,'uint16');
fwrite(fs.map,fs.DAQ.fred,'uint16');
%end
% set size of next block of data
fs.DAQ.data_leftover=frameFinalDataTransposed(:,all_linestarts(129)-round(period*0.1):end);
%fs.DAQ.data_leftover=frameFinalDataTransposed(:,all_linestarts(129):end);

lftle=size(fs.DAQ.data_leftover,2);

%if fs.DAQ.framecounter_thisacq==1
%    fs.DAQ.InPointsToCollect=round(linestart+period*128*2+period*0.1-fs.DAQ.InPointsToCollect);
% else
    fs.DAQ.InPointsToCollect=max(1,round(period*128+period*0.1-lftle));
%end

%end
%}

%% write data to hard drive
if fs.DAQ.DoNotSave==0

    if mod(fs.DAQ.framecounter_thisacq,fs.DAQ.dad_size_in_frames)==1  || (fs.cycle.Do && fs.cycle.NFramesPrStep>0 && mod(fs.DAQ.framecounter_thisacq,fs.cycle.NFramesPrStep)==1) %determine how many blocks are saved in one file
%        subStartNewDad;
fs.DAQ.filecounter=fs.DAQ.filecounter+1;
fs.DAQ.fnmg=[fs.DAQ.BaseFileNameLongGreen '_' sprintf('%06i',fs.DAQ.filecounter) '.tif'];
set(fs.handles.lblCurrentFileName,'String',fs.DAQ.fnmg);
fs.DAQ.fnmr=[fs.DAQ.BaseFileNameLongRed '_' sprintf('%06i',fs.DAQ.filecounter) '.tif'];
fs.DAQ.greenStack = ij.ImageStack(256,240);
fs.DAQ.redStack = ij.ImageStack(256,240);

    end
    %   looks like uint16 takes less time, we can do it if our signal is always positive
    %   DAQ board is 12 bit, so we can safely shift data up (see above)

%    fwrite(fs.fid_curr,frameFinalDataTransposed,'uint16');
greenProcessor=ij.process.ShortProcessor(256,240);

pixels=reshape(fs.DAQ.fgreen,256*240,1);

greenProcessor.setPixels(pixels);

fs.DAQ.greenStack.addSlice('frame',greenProcessor);

redProcessor=ij.process.ShortProcessor(256,240);
pixels=reshape(fs.DAQ.fred,256*240,1);
redProcessor.setPixels(pixels);
fs.DAQ.redStack.addSlice('frame',redProcessor);    

    if mod(fs.DAQ.framecounter_thisacq,fs.DAQ.dad_size_in_frames)==0 || (fs.cycle.Do && fs.cycle.NFramesPrStep>0 && mod(fs.DAQ.framecounter_thisacq,fs.cycle.NFramesPrStep)==0)
%msg='writing'
greenPlus=ij.ImagePlus('GreenStack',fs.DAQ.greenStack);
greenFileSaver=ij.io.FileSaver(greenPlus);
greenFileSaver.saveAsTiffStack(fs.DAQ.fnmg);
redPlus=ij.ImagePlus('RedStack',fs.DAQ.redStack);
redFileSaver=ij.io.FileSaver(redPlus);
redFileSaver.saveAsTiffStack(fs.DAQ.fnmr);

    end
end

[hour, minute, second] = sec2hms(etime(clock,fs.timing.runStartClock));

str = sprintf('%s %02d:%02d:%02d %07d %07d',fs.DAQ.modeString,hour,minute,round(second),fs.DAQ.framecounter,fs.DAQ.pCellCounter-fs.DAQ.framecounter);
set(fs.handles.lblProgress,'string',str);

% %% update log file
% if ~fs.DAQ.DoNotSave
%     if is_new_file
%         fprintf(fs.DAQ.fdLog, 'Opened file %s; %5d frames in %8.3f sec\n', ...
%             fs.DAQ.CurrDadNameBare, fs.DAQ.framecounter, second);
%     end
% end


%% if in z stack, see if it's time to move stage
if fs.cycle.Do
    if fs.cycle.NFramesPrStep>0
        if mod(fs.DAQ.framecounter_thisacq,fs.cycle.NFramesPrStep)==0
            %disp(fs.DAQ.framecounter_thisacq);
            StopDAQ;

            % set time to restart
            fs.DAQ.zNextStartTime = currT + dad_fclose_delay_ms/1000;
            fs.DAQ.zIsWaiting = true;

            % log position
            if ~fs.DAQ.DoNotSave
                fprintf(fs.DAQ.fdLog, 'Z step done. Z pos: %d, framecounter %d\n', fs.position.z, fs.DAQ.framecounter);
            end

            if fs.DAQ.framecounter<fs.cycle.NFramesPrStep*fs.cycle.NSteps*fs.cycle.NReps
                %change parameters inside loop here - TODO
                loopN=ceil(fs.DAQ.framecounter/(fs.cycle.NFramesPrStep*fs.cycle.NSteps));
                stepN=ceil((fs.DAQ.framecounter-(loopN-1)*fs.cycle.NFramesPrStep*fs.cycle.NSteps)/fs.cycle.NFramesPrStep);
                if fs.stage.ready && fs.cycle.ZStep~=0
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
                % when we exit this loop, we will have a delay to wait for
                % last file to close
            else
                %end of stack, stop aquisition
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
%toc
%{
wr_duration=toc; %SY1
    if mod(fs.DAQ.framecounter_thisacq,fs.DAQ.dad_size_in_frames)==0 
   wr_duration    
    end
%}
%% return (to timer) or stop acquisition
if fs.DAQ.Running
    return;
end

%% stop acquisition
frewind(fs.map); %SY commented for test only!
fwrite(fs.map,0,'uint32'); %image size in DAQ points

StopAcquisition;

if fs.DAQ.TestMode==1
    %    figure;
    %    plot(fs.tictoc(:));
    figure;
    %   plot(fs.sAvailable(:),'r');
    plot(fs.phase(:),'r');
    %    figure;
    %    plot(fs.sAvailableFrames(:),'g');

end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function subStartNewDad
% % begin a new dad file on disk; called for first and all subsequent dads
% global fs
% 
% fs.DAQ.filecounter=fs.DAQ.filecounter+1;
% %fnm=[fs.DAQ.BaseFileNameLong num2str(fs.DAQ.filecounter) '.dad'];
% fnm = makeDadName(fs.DAQ.BaseFileName, fs.DAQ.filecounter);
% set(fs.handles.lblCurrentFileName,'String',fnm);
% assert(isempty(fs.fid_prev), 'bug: Trying to start new file but prev. not yet closed');
% fs.fid_prev=fs.fid_curr;
% fs.fid_curr= fopen(fnm,'W');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

