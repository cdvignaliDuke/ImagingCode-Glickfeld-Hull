function p=ProcessStream(object, event)
% called every 5 ms;

dad_size_in_frames = 150;
dad_fclose_delay=55;

global fs;

if fs.cycle.NFramesPrStep>0 %each data files contains data for one z-slice
    dad_size_in_frames = min(dad_size_in_frames,fs.cycle.NFramesPrStep);
    dad_fclose_delay=min(dad_fclose_delay,fs.cycle.NFramesPrStep-1);
end


%% next part was to restart pCell waveform, commented by SY
%{
a=strfind(fs.DAQ.aoPockels.Running,'On');
if isempty(a)
    if mod(fs.DAQ.framecounter,75)==10
        restartPCellWaveform
    end
end
%} 

% fs.DAQ.clockcounter=fs.DAQ.clockcounter+1;
%     if fs.DAQ.TestMode==1 % if in test mode
%         fs.sAcquired(fs.DAQ.clockcounter)=fs.DAQ.ai.SamplesAcquired; 
%         fs.sAvailable(fs.DAQ.clockcounter)=fs.DAQ.ai.SamplesAvailable;
%         fs.sAvailableFrames(fs.DAQ.clockcounter)=fs.DAQ.ai.SamplesAvailable/fs.DAQ.InPointsToCollect;
%  %       fs.tictoc(fs.DAQ.clockcounter)=t;
% %         if fs.DAQ.clockcounter==1
% %             fs.totaltime(1)=t;
% %         else
% %             fs.totaltime(fs.DAQ.clockcounter)=fs.totaltime(fs.DAQ.clockcounter-1)+t;
% %         end
%     end


%% wait for complete frame
sampAv=fs.DAQ.ai.SamplesAvailable;
%sampAv
% beware: fs.DAQ.InPointsToCollect varies from one call to another.
if sampAv < fs.DAQ.InPointsToCollect
    return;
end

%% process frame
fs.DAQ.framecounter=fs.DAQ.framecounter+1; 

frameFinalDataTransposed=(getdata(fs.DAQ.ai, fs.DAQ.InPointsToCollect, 'native') + 2048)';

if sampAv < fs.DAQ.InPointsToCollect*100 %prepare and send data to display if have time only
    if fs.DAQ.TestMode==1 % if in test mode
        fs.sAcquired(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAcquired; 
        fs.sAvailable(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAvailable;
        fs.sAvailableFrames(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAvailable/fs.DAQ.InPointsToCollect;
%         fs.tictoc(fs.DAQ.framecounter)=t;
%         if fs.DAQ.framecounter==1
%             fs.totaltime(1)=t;
%         else
%             fs.totaltime(fs.DAQ.framecounter)=fs.totaltime(fs.DAQ.framecounter-1)+t;
%         end
    end
    % find zero crossings
    xpos = int16(frameFinalDataTransposed(3,:))-2048;
    sg = xpos>=0;
    dsg=diff(sg);
    tcrossup=find(dsg>0); % sample before upward zero crossing

    % calculate period
    period=(tcrossup(end)-tcrossup(1))/(length(tcrossup)-1);

    %calculate first min of x position (line start)
    if tcrossup(1)>period/4
        linestart=round(tcrossup(1)-period/4);
    else
        linestart=round(tcrossup(1)+period*3/4);
    end
    if fs.DAQ.framecounter==1
        fs.DAQ.linestart=linestart;
        fs.DAQ.period=period;
    end
    
    %for hdw trigger test - remove later
   fs.phase(fs.DAQ.framecounter)=tcrossup(1);
    % write frame to ram disk
	frewind(fs.map);
	fwrite(fs.map,round(period*120),'uint32'); %image size in DAQ points
	fwrite(fs.map,linestart,'uint32'); %phase shift
    fwrite(fs.map,frameFinalDataTransposed,'uint16');
    
    % set size of next block of data
    if fs.DAQ.framecounter==1
        fs.DAQ.InPointsToCollect=round(linestart+period*128*2-100-fs.DAQ.InPointsToCollect);
    else
        fs.DAQ.InPointsToCollect=round(period*128+(tcrossup(1)-400)/4);
    end

end

%% write data to hard drive
if fs.DAQ.DoNotSave==0
    if mod(fs.DAQ.framecounter,dad_size_in_frames)==1 %determine how many blocks are saved in one file
        fs.DAQ.filecounter=fs.DAQ.filecounter+1;
        fnm=[fs.DAQ.BaseFileNameLong num2str(fs.DAQ.filecounter) '.dad'];
        set(fs.handles.lblCurrentFileName,'String',fnm);
        fs.fid_prev=fs.fid_curr;
        fs.fid_curr= fopen(fnm,'W');
    end
    %   looks like uint16 takes less time, we can do it if our signal is always positive
    %   DAQ board is 12 bit, so we can safely shift data up (see above)
    
    fwrite(fs.fid_curr,frameFinalDataTransposed,'uint16');
    if mod(fs.DAQ.framecounter,dad_size_in_frames)==dad_fclose_delay
        tic;
        if ~isempty(fs.fid_prev)
            fclose(fs.fid_prev);
        end
        t=toc;
        fs.tictoc(fs.DAQ.framecounter)=t;
    end
end

[hour, minute, second] = sec2hms(etime(clock,fs.timing.startime));

str = sprintf('%s %02d:%02d:%02d %07d %07d',fs.DAQ.modeString,hour,minute,round(second),fs.DAQ.framecounter,fs.DAQ.pCellCounter-fs.DAQ.framecounter);
set(fs.handles.lblProgress,'string',str);

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

%% stop acquisition
frewind(fs.map);
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
