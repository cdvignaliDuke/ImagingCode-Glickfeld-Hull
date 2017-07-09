function [counterValues, counterTimesMs] = counter_fixer_2P(counterValues, counterTimesMs)
%This function takes the frame times and values and corrects for any
%misalignment issues. It overwrites the old vectors. 

% remove duplicate counters
[~,index] = unique(counterValues,'first');
counterValues = counterValues(sort(index));
counterTimesMs = counterTimesMs(sort(index));
counterTimesMs(counterValues == 0) = [];
counterValues(counterValues == 0) = [];

counterValue_temp = min(counterValues):max(counterValues);

counterValue_miss = setdiff(counterValue_temp, counterValues) - min(counterValues) + 1;

counterTime_miss = counterTimesMs(counterValue_miss - int64(1:length(counterValue_miss))) + mode(diff(counterTimesMs));

counter_logic = false(1,length(counterTimesMs)+length(counterTime_miss));

counter_logic(counterValue_miss)=true;

counterTimes = nan(size(counter_logic));

counterTimes(~counter_logic) = counterTimesMs;

counterTimes(counter_logic) = counterTime_miss;

counterValues = counterValue_temp;
counterTimesMs = counterTimes;
end

