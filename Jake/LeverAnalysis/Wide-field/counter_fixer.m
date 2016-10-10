function [counterValuesFixed, counterTimesMsFixed] = counter_fixer(counterValues, counterTimesMs);
%generate a suite of functions. Need one to insert missed frames in HAD_2P_frames. For HAD_2P_frames all the missing frames happen inbetween trials. 
%need one to correct for articfacts found in img53's datasets
%may need more to correct for artifacts found before labjack was updated
%------------------------------
%check for a particular type of error which occured with img52_160725.
%Imaging started on first trial. 
if counterValues(1:4) == [0,1,0,1];
    counterValues = counterValues(4:end);
    counterTimesMs = counterTimesMs(4:end)
end
%------------------------------
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
%----------------------------------
counterValuesFixed = counterValues;
counterTimesMsFixed = counterTimesMs;
assert(ismonotonic(counterValuesFixed, 'Warning: counter_fixer was unable to fix counterValues'));
assert(ismonotonic(counterTimesMsFixed, 'Warning: counter_fixer was unable to fix counterTimesMs'));
end
