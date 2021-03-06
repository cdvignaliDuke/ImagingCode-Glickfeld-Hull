% processing for cerebellarStim(_2P) data that processes errors in trial
% length that point to dropped or added frames and builds an index

dataStruct = struct;
%% formatting assumes "input" is a cerebellarStim data set
ds = input;
% and also enforces it
assert(isfield(ds, 'spCounter1TimesUs'), 'Input File is not a cerebellarStim data set');
% takes tr by tr
cTimes = input.counterTimesUs;
cVals = input.counterValues;
%  compresses all times/vals
totalTimes = cell2mat(cTimes);
totalVals =  cell2mat(cVals);
totalTimes = totalTimes(3:end); % removes start-stop-start artifact and 0 value initial state value
totalVals= totalVals(3:end); % removes start-stop-start artifact and 0 value initial state value
lastTrCounter = cell2mat(input.counter);

% extract all of our stimulus crap
stimOnFrame = cell2mat(input.ttCounter);
blank = zeros(1, length(totalVals))
blank(stimOnFrame)=1;
stimOnIx = blank; 


imagingRate = input.frameImagingFrequencyHz; % in Hz
frameIntUs = 1000000/imagingRate;
%% Check that counter was not reset mid-trial
assert(ismonotonic(totalTimes), 'counterTimesUs is not monotonically increasing: was experiment reset?');
assert(ismonotonic(totalVals), 'counterValues is not monotonically increasing: was experiment reset?');
%%
nTrs = length(input.counter);
for i=2:nTrs,
    lastVal = cVals{i-1}(end);
    thisTrVal = cVals{i}(1);
    framesBetweenTrs(i-1) = thisTrVal - lastVal;
end
pulseDiff = diff(totalTimes);
pulsePauseIx = pulseDiff>(1.5*frameIntUs);
problematic = pulseDiff(pulsePauseIx);
problemFrameNumbers = find(pulsePauseIx);
normal = pulseDiff(~pulsePauseIx);
avgFrame = mean(normal);

howManyMissed = problematic./(2*avgFrame);
% alert if data has a hole greater than 1 second
assert(sum(howManyMissed>imagingRate)<1, '*** Do not use data: there is a skip of %3d second skip between frames. ****', problematic(howManyMissed>15)/1000000);
hold on
%% Replicate graphing to allow frame skipping feedback
if length(problemFrameNumbers)>0,
    [a b c] = plotyy(1:length(pulseDiff), pulseDiff, problemFrameNumbers, howManyMissed)
    title('Estimates of Missed Frames')
    set(a(1), 'Ylim', [0 max(pulseDiff)])
    set(a(2), 'Ylim', [0 max(howManyMissed)])
    set(c, 'Marker', 'x')
    set(c, 'LineStyle', 'none')
    ylabel('Inter-Frame Latency (us)')
    xlabel('Frames')
    set(get(a(2),'Ylabel'),'String','Frames Missed')
else
    plot(1:length(pulseDiff), pulseDiff)
    ylim([0 max(pulseDiff)])
    xlabel('Frames')
    ylabel('Inter-Frame Latency (us)')
end

%% Generate the good frames index
shift = 0;
blankIx = ones(1,length(totalVals)+sum(howManyMissed));
for i = 1:length(problemFrameNumbers),
    insertFrames = howManyMissed(i);
    startFrame = problemFrameNumbers(i)+1; %to correct for 1 frame "diff" offset
    blankIx((startFrame+shift):(startFrame+insertFrames+shift)) = 0;
    stimOnIx((startFrame+shift):(startFrame+insertFrames+shift)) = 0;
    shift = shift+insertFrames;
end
dataStruct.goodFramesIx = boolean(blankIx);
dataStruct.stimOnIx = boolean(stimOnIx)
%% WRITE CODE HERE TO CROP OUT TRIALS WITH 1 SECOND HOLES

%% WRITE CODE HERE TO IMPORT IMAGE
%dataStruct.img = imread('');
%% WRITE CODE HERE TO INDEX IMAGE

%% WRITE CODE HERE TO CHECK IF IMAGE HAS THE SAME NUMBER OF FRAMES AS BEHAVIOR FILE HAS IN COUNTERS

%% IF doLED ACTIVE, SEARCH FOR HIGH LUMINANCE SPIKES

% if isfield(input, 'doLED')
%     if input.doLED==1,
    % do LED searchy type things: probably something like "all the luminous frames please stand up"
    % and then correlate that with trial starts/ends to see if you're all
    % matched up
%     end
% end
