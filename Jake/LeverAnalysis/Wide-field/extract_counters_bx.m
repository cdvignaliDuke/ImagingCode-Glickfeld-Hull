function [counterValues, counterTimesMs, f_frame_trial_num, l_frame_trial_num] = extract_counters_bx(b_data);
%extracts frame counter times and values from the behavior file
%checks to make sure the trial actually has frame counters
%Excludes spurious counters from all non-imaged trials
%Does not address issues with repeating counters or repeating values of 0
counterValues = [];
counterTimesUs = []; 
for i = 1:length(b_data.counterValues)
    counters_this_trial = b_data.counterValues{i};
    if sum(counters_this_trial) > 2;   %sums counter values. Must be greater than 2 to ensure counters are actually logging frame # and not spurious counts that reset (eg 0,1,0)
        if length(counters_this_trial)>1; %after imaging has stopped some trials have a single counter value even though there were no frames
            if isempty(counterValues) 
                f_frame_trial_num = i;  %stores the trial num of the first trial with valid counters. 
            end
            l_frame_trial_num = i; %stores the trial num of the last trial with valid counter values
            counterValues  = [counterValues, counters_this_trial];
            counterTimesUs = [counterTimesUs, b_data.counterTimesUs{i}];
        end
    end
end
counterTimesMs = round(counterTimesUs/1000);
end