function [counterValuesFixed, counterTimesMsFixed] = counter_fixer(counterValues, counterTimesMs);
%This function takes the frame times and values and corrects for any
%misalignment issues. It overwrites the old vectors. 

%check to see if the timecourses are properly aligned. If so then exit the function
if ismonotonic(counterValues, 1) & ismonotonic(counterTimesMs);
    if max(diff(counterValues))==1 & min(diff(counterValues))==1
        if max((counterTimesMs))<110 & min(diff(counterTimesMs))>85;
            disp('counterValues and counterTimesMs already well aligned. Check plots to make sure.')
            figure; plot(diff(counterTimesMs)); title('difference between counter TIMES'); xlabel('frame #'); ylabel('ms');
            figure; plot(diff(counterValues)); title('difference between counter VALUES'); xlabel('frame #'); ylabel('frame # diff');
            counterValuesFixed = counterValues;
            counterTimesMsFixed = counterTimesMs; 
            return
        end
    end
end
disp('error detected in counterTimes or counterValues: diagnosing');
%plot the data to show the confounds in there original form. 
figure; plot(diff(counterTimesMs)); title('difference between counter TIMES before correction'); xlabel('frame #'); ylabel('ms');
figure; plot(diff(counterValues)); title('difference between counter VALUES before correction'); xlabel('frame #'); ylabel('frame # diff');

%check for a particular type of error which occured with img52_160725.
if counterValues(1:4) == [0,1,0,1];  %Imaging started on first trial. 
    counterValues = counterValues(4:end);
    counterTimesMs = counterTimesMs(4:end);
end

%remove zeroes from counter and counter times
IFI = mode(diff(counterTimesMs));  % finds the InterFrame Interval
if max(diff(counterValues))>1 | min(diff(counterValues))<1 | min(diff(counterTimesMs))<(IFI*0.9);
    zeroInd = find(counterValues==0); %find and remove all counter values = 0
    if ~isempty(zeroInd)
        counterValues(zeroInd) = [];
        counterTimesMs(zeroInd)= [];
        disp('Zeroes Removed: counter values of 0 located and removed along with their associated counter times');
    end
    repeat_ind = [];
    for i = 2:length(counterValues)
        if find(counterValues(i)==counterValues(i-1)); %finds repeated counter values and associated counterTimes
            repeat_ind = [repeat_ind, i];
        end
    end
    if ~isempty(repeat_ind);
        counterValues(repeat_ind) = [];
        counterTimesMs(repeat_ind)= [];
        disp('Repeats Removed: counter values which repeat were located and removed along with their associated counter times')
    end
end

%check for a specific type of error occuring in 161212_img73 cue reward pairing
diff_counter_times = diff(counterTimesMs);
shortInts = find(diff_counter_times < IFI*0.6 & diff_counter_times > IFI*0.4);
for ii = length(shortInts)-1:-1:1
    if shortInts(ii)+1 == shortInts(ii+1) %if there are two consecutive short intervals
        diff_counter_values = diff(counterValues);
        if diff_counter_values(shortInts(ii)) == 2  %if there is also a skipped counter value at this time...
            counterTimesMs(shortInts(ii+1)) = [];  
            counterValues((shortInts(ii)+1):end)  = counterValues((shortInts(ii)+1):end)-1;     
            counterValues(end) = []; 
        end
    end
end

%check for skipped frame times and counters. Fill them in.
if max(diff(counterValues)) == 2; 
    skipped_counters = find(diff(counterValues)==2);
    skipped_times = find(diff(counterTimesMs)>180 & diff(counterTimesMs)<220);
    if length(skipped_counters) == length(skipped_times)
        for ii = length(skipped_counters):-1:1;
            counterValues = [counterValues(1:skipped_counters(ii)), counterValues(skipped_counters(ii))+1, counterValues(skipped_counters(ii)+1:end)];
            counterTimesMs = [counterTimesMs(1:skipped_counters(ii)), counterTimesMs(skipped_counters(ii))+IFI ,counterTimesMs(skipped_counters(ii)+1:end)];
        end
        disp('skipped frames located and filled in')
    else
        error('number of skipped frame times does not equal number of skipped frame values')
    end
end

%error found in 160129_img36. 
diff_counter_times = diff(counterTimesMs);
frame_osc = find(diff_counter_times==110) ;
if ~isempty(frame_osc)
    for ii = 1:length(frame_osc)
        if diff_counter_times(frame_osc(ii)+1)<95;
            counterTimesMs(frame_osc(ii)+1) = counterTimesMs(frame_osc(ii)+1) - (diff_counter_times(frame_osc(ii))-IFI);
        end
    end
end

%corrections specific to the datasets collected when labjack only updated
%every 15ms. 
diff_counter_times = diff(counterTimesMs);
if length(counterTimesMs)*.2 < length(find(diff_counter_times<(IFI-10))); %Identify old datssets: if more than 20% of IFIs are 10ms less than the most frequent IFI then this is an old dataset. 
diff_counter_times = diff(counterTimesMs);
shifted_frame_inx = find(diff(counterTimesMs)>115 & diff(counterTimesMs)<125);
if ~isempty(shifted_frame_inx);
    for ii = 1:length(shifted_frame_inx)
        if diff_counter_times(shifted_frame_inx(ii)+1) < 85;
            counterTimesMs(shifted_frame_inx(ii)+1) = counterTimesMs(shifted_frame_inx(ii)+1) - (diff_counter_times(shifted_frame_inx(ii))-100); %subtract off the difference from the IFI
        end
        if diff_counter_times(shifted_frame_inx(ii)+1) < 95 & diff_counter_times(shifted_frame_inx(ii)+2) < 95;
            counterTimesMs(shifted_frame_inx(ii)+1) = counterTimesMs(shifted_frame_inx(ii)+1) - (diff_counter_times(shifted_frame_inx(ii))-IFI);
        end
        if diff_counter_times(shifted_frame_inx(ii)+1) < 95 & diff_counter_times(shifted_frame_inx(ii)-1) < 95 & diff_counter_times(shifted_frame_inx(ii)-2) > IFI-2; % this error found in 151022_img29
            counterTimesMs(shifted_frame_inx(ii)+1) = counterTimesMs(shifted_frame_inx(ii)+1) - (diff_counter_times(shifted_frame_inx(ii))-IFI);
        end
    end
end
end

%if statement for error seen in 161216_img73
if length(counterValues) == 40003 & length(counterTimesMs) == 40003;
    diff_counter_times = diff(counterTimesMs);
    shortVals = find(diff_counter_times< 30);
    if ~isempty(shortVals)
        for ii = length(shortVals):-1:1
            if diff_counter_times(shortVals(ii)+1)<110 & diff_counter_times(shortVals(ii)-1)<110;
                counterTimesMs(shortVals(ii)) = [];
                counterValues(shortVals(ii):end) = counterValues(shortVals(ii):end)-1
                counterValues(shortVals(ii)) = [];
            end
        end
        %store values and check to see if alignment worked
        counterValuesFixed = counterValues;
        counterTimesMsFixed = counterTimesMs;
        assert(length(counterTimesMs) == length(counterValues));
        return
    end
end

%store values and check to see if alignment worked
counterValuesFixed = counterValues;
counterTimesMsFixed = counterTimesMs;
assert(length(counterTimesMs) == length(counterValues));

%if the vectors pass all these if statements they should be properly aligned.
if ismonotonic(counterValues, 1) & ismonotonic(counterTimesMs);
    if max(diff(counterValues))==1 & min(diff(counterValues))==1
        if max(diff(counterTimesMs))<110 & min(diff(counterTimesMs))>85;
            figure; plot(diff(counterTimesMs)); title('difference between counter TIMES'); xlabel('frame #'); ylabel('ms');
            figure; plot(diff(counterValues)); title('difference between counter VALUES'); xlabel('frame #'); ylabel('frame # diff');
            return
        else
            error('counterValues or counterTimes not well aligned');
        end
    else
        error('counterValues or counterTimes not well aligned');
    end
else
    error('counterValues or counterTimes not well aligned');
end
end

