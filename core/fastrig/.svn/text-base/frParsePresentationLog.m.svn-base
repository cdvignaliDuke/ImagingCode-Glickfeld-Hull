function [stims,blanks,epochs,codes]=frParsePresentationLog(s,nStimIncBlank,nFrames,mode);
% [stims,blanks,epochs]=frParsePresentationLog(s,nStimIncBlank,nFrames,mode);
% where input argument s is log as structure array
%       nStimIncBlank is number of distinct stimuli including blank (see prot struct)
%       nFrames is number of imaged frames
%       mode is presentation mode (see presentation scenario)
%            mode = 0 , no blanks interleaved
%            mode = 1 , 1:1 blanks interleaved, blank first, stim second
%            mode = 2 , 1:1 blanks interleaved, stim first, blank second
%  where output argument stims is array of stimulation frames
%                        blanks is array of blank frames
%                        epochs is (compound) array of stimulus + blanks frames
%

%%
if nargin < 4;mode = 2;end;

%%
% when stimulus presentation started
StartPulseIndex = find(strcmp('Pulse',[s.EventType]));
StartTime = s(StartPulseIndex).Time;
% measured two-photon images
FramePulseIndex = find(strcmp('Port Input',[s.EventType]));
FrameTime = single([s(FramePulseIndex).Time]);

% stimulus frames
PicturePulseIndex= find(strcmp('Picture',[s.EventType]));
PicturePulseIndex(1)=[]; % discard first picture, not part of stimulus
PictureTime = single([s(PicturePulseIndex).Time]);
PictureTimeUncertainty = single([s(PicturePulseIndex).TimeUncertainty]);
PictureCode = single([s(PicturePulseIndex).Code]);

% correct for missed frames
deltas = diff(FrameTime);
missed = find(deltas>625);
FrameTimeCorrected = [];
for index = 1:length(FrameTime)-1
    FrameTimeCorrected(end+1)=FrameTime(index);
    if deltas(index) > 625
        FrameTimeCorrected(end+1) = round(FrameTime(index) + deltas(index) / 2);
    end
end

% correct for hickups
deltas= (diff(FrameTimeCorrected));
ddeltas = [diff(deltas)];
hickups = find(ddeltas>1);
% hickups(find(diff(hickups)<=6)+1)=[];
offsets = ddeltas(hickups)/2;

FrameTimeFiltered = FrameTimeCorrected;

for index = 2:length(hickups)
    sel = hickups(index)+1:hickups(index)+3;
    FrameTimeFiltered(sel) = FrameTimeFiltered(sel) - ddeltas(sel);    
end

figure;
subplot(2,1,1);
plot(diff(FrameTime)/10,'r')
hold on;
plot(diff(FrameTimeCorrected)/10);
xlabel('Frame #');
ylabel('Frame Interval (ms)');
title('Correction for missed frame pulses');
legend('Orig','Corrected');
subplot(2,1,2);
plot(diff(FrameTimeCorrected)/10);
hold on;
plot(diff(FrameTimeFiltered)/10,'k.-');
ylabel('Frame Interval (ms)');
xlabel('Frame #');
title('Pulse time filtering');
legend('Corrected','Filtered');

% estimate of frame period
FramePeriod=(FrameTimeFiltered(end)-FrameTimeFiltered(1))/...
            (length(FrameTimeFiltered)-1);

%%
nPicsTotal= max(PictureCode);
nPicsPerStim = nPicsTotal/nStimIncBlank; % assumes all stim inc blank same number of pics
nRefreshes = length(PictureCode);
nStimExBlank = nStimIncBlank - (mode>0) ;
nTrials = ceil(nRefreshes/nPicsTotal/(1+(mode>0)));

fprintf('nPicsPerStim=%i nPicsTotal=%i nRefreshes=%i nStimExBlank=%i nTrials=%i\n',...
         nPicsPerStim,nPicsTotal,nRefreshes,nStimExBlank,nTrials);


StimulusIDs = zeros(size(PictureCode));

%%
if nPicsPerStim > 1
    
    %% figure out which stimulus was on from the picture codes, this is robust
    % Assumes stimulus IDs and picture codes are related by the following
    %  Picture codes 1 to nPicsPerStim corresponds to Stimulus 1
    %  Picture codes nPicsPerStim+1 to 2(nPicsPerStim corresponds to Stimulus 2
    %  etc
    
    for istim = 1:nStimIncBlank
        sel = PictureCode>=(istim-1)*nPicsPerStim + 1 & PictureCode<=istim*nPicsPerStim;
        StimulusIDs(sel)=istim;  
    end
    
    %% find out when each stimulus was presented
    
    % infer stimulus transitions based on picture codes
    % (missing transitions with contiguous pictures)
    onsetIndA = [1 (find(abs(diff(PictureCode))>=nPicsPerStim*.75))+1];
    % infer transitions based on picture ids 
    % (missing transitions between identical stimuli)
    onsetIndB = [1 find(diff(StimulusIDs))+1];
    % combine them, should be complete
    onsetInd = union(onsetIndA,onsetIndB);  

%%
else
    onsetIndA = find(diff(abs(diff(PictureCode)))>0);
    onsetIndB = find(diff(abs(diff(PictureCode)))<0) + 1;
    onsetIndEst = union(onsetIndA,onsetIndB);   
    
    dOnsetInd  = diff(onsetIndEst);
    
    % a monstruous hack
    % fix missing blank onsets assumming onset interval periodic    
    holes = find(dOnsetInd >= 2*median(dOnsetInd))   ;
    onsetInd = [];
    
    for ihole = 1 :length(holes)
        if ihole > 1
            sel = holes(ihole-1)+1:holes(ihole);
        else
            sel = 1: holes(ihole);
        end
        
        onsetInd = [onsetInd onsetIndEst(sel)];
        
        ds = dOnsetInd(fliplr(sel(1:end-2)));
        
        sel2 = find(cumsum(ds)<dOnsetInd(holes(ihole)));
        
        ds(sel2)
        onsetInd = [onsetInd onsetInd(end)+cumsum(ds(sel2))];
 
    end
    
    onsetInd = [1 onsetInd onsetIndEst(holes(ihole)+1:end)];
    
    
    for ind = 2:length(onsetInd)
        sel = onsetInd(ind-1):onsetInd(ind)-1;
        ucodes = unique(PictureCode(sel));
        
        switch length(ucodes)
            case 1
                StimulusIDs(sel) = ucodes;
            case 2
                StimulusIDs(sel)  = setdiff(ucodes,nPicsTotal);
            otherwise
                error('something wrong happen');
        end
                
    end
    
    onsetInd(end)=[];
    
    figure;   subplot(2,1,1); plot(diff(onsetInd),'.-k');     
    subplot(2,1,2);plot(diff(onsetIndEst),'.-b');
    
  
end

%%
stimuliSeq = StimulusIDs(onsetInd);
onsetTimes = PictureTime(onsetInd);

%% corrects for (some) gaps in presentation log
% first picture displayed at stimulus onset
loggedOnsetPics = PictureCode(onsetInd);

% first picture that should have been presented
expectedOnsetPics = (stimuliSeq-1)*nPicsPerStim+1;

% assumes stimulus was displayed but not logged
dev = loggedOnsetPics-expectedOnsetPics;
gaps = find(dev);

dPictureTime = diff(PictureTime);

% measure refresh period, exclude outliers
pPictureTime = prctile(dPictureTime,[1 99]);
inlier = find(dPictureTime>=pPictureTime(1) & dPictureTime<=pPictureTime(2));
interval = mean(dPictureTime(inlier));

% onset estimation delay
delay = round(dev*interval);

onsetTimesCorrect = onsetTimes-delay;

%% fig summary of stimulus onset calculation
figure;
subplot(3,1,1);plot(PictureTime(onsetIndA(2:end))/10e3,diff(onsetIndA),'r.-');
hold on;
plot(PictureTime(onsetIndB(2:end))/10e3,diff(onsetIndB),'b.-');
plot(PictureTime(onsetInd(2:end))/10e3,diff(onsetInd),'k.-');
legend('picture codes','stimulus ids','output');
ylabel('Stimulus onset interval (indices)');
title('Inferred onset indices');

subplot(3,1,2);
plot(PictureTime(onsetInd(2:end))/10e3,diff(onsetTimes)/10e3,'k');
ylabel('Stimulus onset interval (s)');
title('Measured onset times');

subplot(3,1,3);
plot(PictureTime(onsetInd(2:end))/10e3,diff(onsetTimesCorrect)/10e3,'k');
ylabel('Stimulus onset interval (s)');
title('Corrected onset times');
xlabel('Time (s)');

%%
% time of each frame from trigger
FrameTimeAcq = FrameTimeFiltered(find(FrameTimeFiltered>StartTime));

% what picture on display at each scan frame
FramePicCode = zeros(size(FrameTimeAcq));
for iFrame = 1:nFrames
    sel =find(PictureTime<FrameTimeAcq(iFrame));
    if length(sel)
        FramePicCode(iFrame) = PictureCode(max(sel));
    else
        FramePicCode(iFrame) = max(PictureCode);
    end
end

if abs(length(FrameTimeAcq)-nFrames)>1
    fprintf(1,'WARNING: %i frames acquired vs. %i pulses measured\n',...
            length(FrameTimeAcq),nFrames);
end

% first frame following stimulus onset
for ind = 1: length(onsetTimesCorrect);
    sel= find(onsetTimesCorrect(ind)>=FrameTimeAcq);
    onsetFrameCount(ind)=length(sel)+1;
end

onsetFrameInd = onsetFrameCount+1;onsetFrameInd(end+1)=nFrames;

% takes care of dropped frame at beginning, a hell of a fudge
DummyFrames = 0;
if PictureCode(1)~=1
    fprintf(1,'Warning: First picture presented has code %i. First epoch may be corrupt.\n',PictureCode(1));
    fprintf(1,'Trying to rescue it anyways\n',PictureCode(1));

    onsetFrameInd(1) = 1;
    fprintf(1,'(1)Set first stimulus onset to acq trigger \n');
    
    % number of imaging frames that should have elapsed since acq trigger
    ExpectedFirstPicCode = floor((PictureTime(1)-double(StartTime))/median(diff(PictureTime)))+1;
    DroppedPictures = PictureCode(1)-ExpectedFirstPicCode;
    DummyFrames = round(DroppedPictures*median(diff(PictureTime))/median(diff(FrameTimeAcq)))-1;
    fprintf(1,'(2)%i pictures dropped. Added % i dummy frames\n',DroppedPictures, DummyFrames);
end

onsetFrameInd(1) = onsetFrameInd(1) - DummyFrames;

% checks for dropped frames later on
ucode = (unique(diff(PictureCode)));
bugs = find(abs(ucode)<nPicsPerStim-1&abs(ucode)>1);

for ibug = 1:length(bugs)
    bugPictureInd = find(diff(PictureCode)==ucode(bugs(ibug)));
    fprintf(1,'Timing off by ');
    fprintf(1,'%2.1f ',(PictureTime(bugPictureInd+1)-PictureTime(bugPictureInd))/median(diff(PictureTime))-1)
    fprintf(1,' frames. ');
    fprintf(1,'Presentation dropped %i frames\n',    ucode(bugs(ibug))-1);
end

stims = cell(nTrials,nStimIncBlank);
blanks = cell(nTrials,nStimIncBlank);
epochs = cell(nTrials,nStimIncBlank);
codes = cell(nTrials,nStimIncBlank);

for itrial = 1:nTrials
    for istim = 1:nStimIncBlank        
        switch mode
            case 0 % no blanks
                ind = (itrial-1)*nStimIncBlank+istim;
                stims{itrial,istim}=max(onsetFrameInd(ind):onsetFrameInd(ind+1)-1,1);  
                epochs{itrial,istim}=stims{itrial,istim};                
            case 1 % interleaved blank first
                ind = (itrial-1)*nStimIncBlank*2+(istim-1)*2+1;
                blanks{itrial,istim}=max(onsetFrameInd(ind):onsetFrameInd(ind+1)-1,1);        
                ind = (itrial-1)*nStimIncBlank*2+(istim-1)*2+2;
                stims{itrial,istim}=max(onsetFrameInd(ind):onsetFrameInd(ind+1)-1,1);
                epochs{itrial,istim}=[stims{itrial,istim} blanks{itrial,istim}];
            case 2 % interleaved stim first
                ind = (itrial-1)*nStimIncBlank*2+(istim-1)*2+1;
                stims{itrial,istim}=max(onsetFrameInd(ind):onsetFrameInd(ind+1)-1,1);
                ind = (itrial-1)*nStimIncBlank*2+(istim-1)*2+2;
                blanks{itrial,istim}=max(onsetFrameInd(ind):onsetFrameInd(ind+1)-1,1);        
                epochs{itrial,istim}=[stims{itrial,istim} blanks{itrial,istim}];
        end
        
        codes{itrial,istim}=FramePicCode(epochs{itrial,istim});
    end
end

fprintf('Acquisition pulse interval = %2.1f ms\n',median(diff(FrameTimeAcq))/10);
fprintf('Video refresh interval = %2.1f ms\n',interval/10);

trialdur = inf;
for iStim = 1:nStimIncBlank
    for iTrial = 1:nTrials
        trialdur=min(trialdur,length(epochs{iTrial,iStim}));
    end
end

return;
