function p=ProcessStream(object, event)
% called every 5 ms;

dad_fclose_delay = 55;
dad_fclose_delay_ms = dad_fclose_delay*31;

global fs;

currT = toc; 

%% check if waiting to start acq (as in a new z level)
if fs.DAQ.zIsWaiting
    %fprintf(1, 'waiting: diff %6.3f\n', fs.DAQ.zNextStartTime-currT);
    assert(fs.DAQ.started == false, 'bug: waiting requested but still running');
    if currT < fs.DAQ.zNextStartTime; %make sure that we had enouth time to finish previous file
        return
    else
        % close prev, start acq again, and start a new file
        fs.DAQ.zIsWaiting = false;

        if ~isempty(fs.fid_prev)  % empty if the last file has < dad_fclose_delay frames; is closed in process frame loop below
            fclose(fs.fid_prev);
            fs.fid_prev = [];
        end
        StartDAQ;  % initializes framcounter_thisacq to 0
        return % next iteration into this fcn will handle data
    end
end

%% wait for complete frame
sampAv=fs.DAQ.ai.SamplesAvailable;
% beware: fs.DAQ.InPointsToCollect varies from one call to another.
if sampAv < fs.DAQ.InPointsToCollect
    return;
end
%tic % tmp SY
%% process frame
fs.DAQ.framecounter=fs.DAQ.framecounter+1;
fs.DAQ.framecounter_thisacq = fs.DAQ.framecounter_thisacq + 1;

frameFinalDataTransposed=(getdata(fs.DAQ.ai, fs.DAQ.InPointsToCollect, 'native') + 2048)';

if fs.DAQ.framecounter_thisacq==1 %check trigger position using y signal
    yPosMax=max(frameFinalDataTransposed(4,:));
    yPosMin=min(frameFinalDataTransposed(4,:));
    yPosTshld=yPosMin+(yPosMax-yPosMin)/18;
    ycross=find(frameFinalDataTransposed(4,:)<yPosTshld);
    if ycross(1)>1
        frameFinalDataTransposed=frameFinalDataTransposed(:,ycross(1)+1:end);
        fs.DAQ.InPointsToCollect=fs.DAQ.InPointsToCollect-ycross(1);
    end
end

%if sampAv < fs.DAQ.InPointsToCollect*100 %prepare and send data to display if have time only
if fs.DAQ.TestMode==1 % if in test mode
    fs.sAcquired(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAcquired;
    fs.sAvailable(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAvailable;
    fs.sAvailableFrames(fs.DAQ.framecounter)=fs.DAQ.ai.SamplesAvailable/fs.DAQ.InPointsToCollect;
end
% find zero crossings
%    xpos = int16(frameFinalDataTransposed(3,:))-2048;
% calculate period

%    tcrossup=find(xpos(1:end-6)<=0 & xpos(2:end-5)>0 & xpos(7:end)>0);

if fs.DAQ.framecounter_thisacq==1
    tcrossup=find(frameFinalDataTransposed(3,1:end-6)<=2048 & frameFinalDataTransposed(3,2:end-5)>2048 & frameFinalDataTransposed(3,7:end)>2048);
    dcrossup=diff(tcrossup);
    smdc=length(find(dcrossup<800*fs.DAQ.realinputRate/3300000)); %count false crossings
    smdc2=length(find(dcrossup<20*fs.DAQ.realinputRate/3300000));
    period=(tcrossup(end)-tcrossup(1))/(length(tcrossup)-1-smdc2-(smdc-smdc2)/2);
else
    le=round(min(1000*fs.DAQ.realinputRate/3300000,size(frameFinalDataTransposed,2)));
    firsttcrossup=find(frameFinalDataTransposed(3,1:le-6)<=2048 & frameFinalDataTransposed(3,2:le-5)>2048 & frameFinalDataTransposed(3,7:le)>2048);
    pos100=firsttcrossup(1)+round(100*fs.DAQ.lastperiod);
    lasttcrossup=find(frameFinalDataTransposed(3,pos100-100:pos100+100-6)<=2048 & frameFinalDataTransposed(3,pos100-100+1:pos100+100-5)>2048 & frameFinalDataTransposed(3,pos100-100+6:pos100+100)>2048);
    tcrossup(1)=firsttcrossup(1);
    period=(lasttcrossup(1)+pos100-100-1-firsttcrossup(1))/100;
end

fs.DAQ.lastperiod=period;
firstpoint=min(tcrossup(1)-1,10);
bbb=regress(double(frameFinalDataTransposed(3,tcrossup(1)-firstpoint:tcrossup(1)+10)-2048)',[ones(firstpoint+11,1) (-firstpoint:10)']);
tcrossup(1)=tcrossup(1)-bbb(1)/bbb(2);

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

if fs.DAQ.TestMode==1
    %for hdw trigger test - remove later
    fs.phase(fs.DAQ.framecounter)=tcrossup(1);
end

% write frame to ram disk
frewind(fs.map);
fwrite(fs.map,fs.DAQ.framecounter,'uint32'); %current frame number
fwrite(fs.map,round(period*120),'uint32'); %image size in DAQ points
fwrite(fs.map,linestart,'uint32'); %phase shift
fwrite(fs.map,frameFinalDataTransposed,'uint16');

% set size of next block of data

if fs.DAQ.framecounter_thisacq==1
    fs.DAQ.InPointsToCollect=round(linestart+period*128*2-period*0.1-fs.DAQ.InPointsToCollect);
else
    fs.DAQ.InPointsToCollect=round(period*128+(tcrossup(1)-period*0.4)/4);
end
%end
%}
%% write data to hard drive
if fs.DAQ.DoNotSave==0
    if mod(fs.DAQ.framecounter_thisacq,fs.DAQ.dad_size_in_frames)==1 %determine how many blocks are saved in one file
        subStartNewDad;
    end
    %   looks like uint16 takes less time, we can do it if our signal is always positive
    %   DAQ board is 12 bit, so we can safely shift data up (see above)

    fwrite(fs.fid_curr,frameFinalDataTransposed,'uint16');
    if mod(fs.DAQ.framecounter_thisacq,fs.DAQ.dad_size_in_frames)==dad_fclose_delay
        if ~isempty(fs.fid_prev)
            fclose(fs.fid_prev);
            fs.fid_prev = [];
        end
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
%toc %SY
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

function subStartNewDad
% begin a new dad file on disk; called for first and all subsequent dads
global fs

fs.DAQ.filecounter=fs.DAQ.filecounter+1;
%fnm=[fs.DAQ.BaseFileNameLong num2str(fs.DAQ.filecounter) '.dad'];
fnm = makeDadName(fs.DAQ.BaseFileName, fs.DAQ.filecounter);
set(fs.handles.lblCurrentFileName,'String',fnm);
assert(isempty(fs.fid_prev), 'bug: Trying to start new file but prev. not yet closed');
fs.fid_prev=fs.fid_curr;
fs.fid_curr= fopen(fnm,'W');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

