function [out]=parsePresentationLog(s);
% [out]=parsePresentationLog(s);
% out.PictureTimes,out.PictureCodes,out.FrameTimes,out.PictureSeq

% extract experiment info
ManIndices = find(strcmp('Manual',[s.EventType]));
ManCodes = [s(ManIndices).Code];

[ParameterNames,r]=strtok(ManCodes,'=');
[discard,r]=strtok(r);
[ParameterValues,r]=strtok(r);

ParameterNames = deblank(ParameterNames);
ParameterValues = deblank(ParameterValues);

for ind = 1:length(ParameterValues)
    arr = str2num(ParameterValues{ind});
    if ~isempty(arr)
       ParameterValues{ind}=arr; 
    end
end

info =  cell2struct(ParameterValues,ParameterNames,2);

out.info = info;

% when stimulus presentation started
PulseIndices = find(strcmp('Pulse',[s.EventType]));
StartPulseIndex = PulseIndices(1);
out.StartTime = double(s(StartPulseIndex).Time);
out.StartTimeUncertainty = double(s(StartPulseIndex).TimeUncertainty);

% two-photon frames
FramePulseIndex = PulseIndices(2:end);
FrameTime = double([s(FramePulseIndex).Time]);

% stimulus frames
PicturePulseIndex= find(strcmp('Picture',[s.EventType]));
PictureTime = double([s(PicturePulseIndex).Time]);
PictureTimeUncertainty = double([s(PicturePulseIndex).TimeUncertainty]);
out.PictureCodes = double([s(PicturePulseIndex).iFrame]);
out.PictureTrials = double([s(PicturePulseIndex).iTrial]);
out.PictureStims = double([s(PicturePulseIndex).iStim]);
out.PictureBlanks = double([s(PicturePulseIndex).bBlank]);

deltas = diff(FrameTime);
missed = find(deltas>625);

if length(missed)
    fprintf(1,'WARNING: Some frame pulses have been missed\n!');
end
    
FrameTimeCorrected = FrameTime;

% % correct for missed frames
% FrameTimeCorrected = [];
% for index = 1:length(FrameTime)-1
%     FrameTimeCorrected(end+1)=FrameTime(index);
%     if deltas(index) > 625
%         FrameTimeCorrected(end+1) = round(FrameTime(index) + deltas(index) / 2);
%     end
% end

% % correct for hickups
% deltas= (diff(FrameTimeCorrected));
% ddeltas = [diff(deltas)];
% hickups = find(ddeltas>1);
% % hickups(find(diff(hickups)<=6)+1)=[];
% offsets = ddeltas(hickups)/2;

FrameTimeFiltered = FrameTimeCorrected;

% for index = 2:length(hickups)
%     sel = hickups(index)+1:hickups(index)+3;
%     FrameTimeFiltered(sel) = FrameTimeFiltered(sel) - ddeltas(sel);    
% end
% 
% figure;
% subplot(2,1,1);
% plot(diff(FrameTime)/10,'r')
% hold on;
% plot(diff(FrameTimeCorrected)/10);
% xlabel('Frame #');
% ylabel('Frame Interval (ms)');
% title('Correction for missed frame pulses');
% legend('Orig','Corrected');
% subplot(2,1,2);
% plot(diff(FrameTimeCorrected)/10);
% hold on;
% plot(diff(FrameTimeFiltered)/10,'k.-');
% ylabel('Frame Interval (ms)');
% xlabel('Frame #');
% title('Pulse time filtering');
% legend('Corrected','Filtered');

% estimate of frame period
FramePeriod=(FrameTimeFiltered(end)-FrameTimeFiltered(1))/...
            (length(FrameTimeFiltered)-1);

%%
nPicsTotal= max(out.PictureCodes);
%nPicsPerStim = nPicsTotal/nStimIncBlank; % assumes all stim inc blank same number of pics
% nRefreshes = length(out.PictureCodes);
% nStimExBlank = nStimIncBlank - (mode>0) ;
% nTrials = ceil(nRefreshes/nPicsTotal/(1+(mode>0)));

% fprintf('nPicsPerStim=%i nPicsTotal=%i nRefreshes=%i nStimExBlank=%i nTrials=%i\n',...
%          nPicsPerStim,nPicsTotal,nRefreshes,nStimExBlank,nTrials);

% time of each frame from trigger
out.FrameTimes = FrameTimeFiltered(find(FrameTimeFiltered>out.StartTime));
out.PictureTimes = PictureTime(find(PictureTime>out.StartTime));
out.PictureTimesUncertainty = PictureTimeUncertainty(find(PictureTime>out.StartTime));

out.FramePicCodes = zeros(size(out.FrameTimes));
out.FrameStims = zeros(size(out.FrameTimes));
out.FrameTrials = zeros(size(out.FrameTimes));
out.FrameBlanks = zeros(size(out.FrameTimes));

nFrames = length(out.FrameTimes);

for iFrame = 1:nFrames
    
    % find all pics presented up to current scan frame    
    sel =find(out.PictureTimes<out.FrameTimes(iFrame));
    % picture that was on during scan frame
    this = sel(end);
    
    out.FramePicCodes(iFrame) = out.PictureCodes(this);
    out.FrameStims(iFrame) = out.PictureStims(this);
    out.FrameTrials(iFrame) = out.PictureTrials(this);
    out.FrameBlanks(iFrame) = out.PictureBlanks(this);
    
end %

return;
