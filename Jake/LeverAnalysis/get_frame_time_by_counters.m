function [ frame_times frame_count sub_input ] = get_frame_time_by_counters(input, info)
% get the frame timing information from counter event times
%Function editted by JH 3/18/16 to make flexible for old and new 2P
%lever data sets. Older datsets ran when the labjack updated every 15ms
%have sever counter value/times alignment issues when imaging at 30Hz
 nframes = info.config.frames;
    % find start of collection
    ntrials = size(input.trialOutcomeCell,2);   
    for i = 1:ntrials
        if length(input.counterValues{i})>2;
            trial_start = i;
            break
        end
    end
    
    % find end of collection
    trial_end = 1;
    for i = ntrials:-1:1
        if find(input.counterValues{i}==nframes,1,'last');
            trial_end = i;
            break
        end
    end
    % in case behavior ends before imaging
    if trial_end == 1;
        trial_end = ntrials;
    end
    
    % create structure with only imaged trials
    sub_input = trialChopper(input, [trial_start trial_end]);
    ntrials = size(sub_input.trialOutcomeCell,2);
    
    % find all frame times and counts
    counterTimesMs = [];
    counterValues = [];
    for i = 1:ntrials
        counterTimesMs = [counterTimesMs round(sub_input.counterTimesUs{i}./1000)];
        counterValues = [counterValues sub_input.counterValues{i}];
    end
    
    %the mean ifi for datasets using the old labjack refresh rate(15ms) is not
    %the most mode ifi. This occurs because of the oscillation between 30ms and
    %45ms when using an intended FR of 30Hz (33.3ms)
    modeDiff = mode(diff(counterTimesMs(1:400)));
    meanDiff = mode(diff(counterTimesMs(1:400)));
    assert(modeDiff<1.2*meanDiff & modeDiff>0.8*meanDiff);
    %All diff(counterTimesMs) values that are 2.5 times as large as the should be
    %have an associated missed counter value. Inserting counterValues and
    %estimated counterTimes to eliminate these outliers.
    missedFrames = find(diff(counterTimesMs)>(2.1*modeDiff));
    missingTimes = ceil((counterTimesMs(missedFrames+1)-counterTimesMs(missedFrames))./2)+counterTimesMs(missedFrames);
    %insert the missing counter value/time into counterValues/TimesMs
    insertionShift=0;
    for i=1:length(missedFrames);
        counterValues = [counterValues(1:[missedFrames(i)+insertionShift]), counterValues(missedFrames(i)+insertionShift)+1, counterValues([missedFrames(i)+insertionShift+1]:end)];
        counterTimesMs = [counterTimesMs(1:[missedFrames(i)+insertionShift]), missingTimes(i), counterTimesMs([missedFrames(i)+insertionShift+1]:end)];
        insertionShift = insertionShift+1;
    end
    
    %at this point every single missed frameCount or repeated frameCount has an
    %associated long or short diff(counterTimesMs) value. Using the
    %diff(counterValues) outliers to ID the diff(counterTimesMs) 
    missedCounters = find(diff(counterValues)>1.5);
    missingTimes = ceil((counterTimesMs(missedCounters+1)-counterTimesMs(missedCounters))./2)+counterTimesMs(missedCounters); %overwritting missingTimes here.
    
    % There are rare occurances where the diff(counterTimesMs) is twice
    %what it should be but does not have an associated missed counter. These
    %are always immmediately followed by a diff(counterTimesMs) value of 0.5
    %its typical value 
    insertionShift=0;
    for i=1:length(missedCounters);
        counterValues = [counterValues(1:[missedCounters(i)+insertionShift]), counterValues(missedCounters(i)+insertionShift)+1, counterValues([missedCounters(i)+insertionShift+1]:end)];
        counterTimesMs = [counterTimesMs(1:[missedCounters(i)+insertionShift]), missingTimes(i), counterTimesMs([missedCounters(i)+insertionShift+1]:end)];
        insertionShift = insertionShift+1;
    end
    
    % Now every single repeated counterValue has an associated short counter
    % Time. Using the location of these repeated counter values to remove
    % both that value and the associated shortened counter time.
    repeatedCounters = find(diff(counterValues)<0.5);
    counterValues(repeatedCounters)=[];
    counterTimesMs(repeatedCounters)=[];
    
    %Now the diff(counterValues)=1 for all frames. However there are still some
    %isolated cases where a diff(counterTimesMs) value is twice what is should
    %be and the following value is 0.5 times what it should be. Fixing this so that 
    %the first value is the standard 1.5*ifi error and the second value is equal 
    %to the ifi. This will more accurately track the real frame times.
    diffTimes = diff(counterTimesMs);
    oscillations = find(diff(counterTimesMs)>(modeDiff*1.75));
    for i = 1:length(oscillations);
        assert(diffTimes(oscillations(i)+1)<(modeDiff*0.66));
    end
    counterTimesMs(oscillations+1) = counterTimesMs(oscillations+1)-modeDiff/2;
    
    %there are also smaller oscillations that occure with a 1.5, 0.5, 1.5 times
    %ifi pattern. Fixing these by reasssigning the first 1.5*ifi value to a
    %value=ifi. Thus restoring a 1, 1, 1.5 pattern.
    diffTimes = diff(counterTimesMs);
    minorOsc = find(diffTimes(1:end-2)<(modeDiff*0.66));
    for i =1:length(minorOsc);
        assert(diffTimes(minorOsc(i)-1)>[modeDiff*1.2]);
    end
    counterTimesMs(minorOsc) = counterTimesMs(minorOsc)-modeDiff/2;
    
    %adjust for extra two counters at end of session
    if counterValues(end) == info.config.frames+2
        counterValues = counterValues(1:end-2);
        counterTimesMs = counterTimesMs(1:end-2);
    end
    
    %zero time stamp to start of first frame
    frame_times = counterTimesMs-counterTimesMs(1);
    frame_count = counterValues;
end

