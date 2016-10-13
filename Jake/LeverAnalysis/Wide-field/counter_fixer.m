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

