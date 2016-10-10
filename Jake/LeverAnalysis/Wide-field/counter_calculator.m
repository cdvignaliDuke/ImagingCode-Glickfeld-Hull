function counter = counter_calculator(counterValues, counterTimesMs);
%calculating frame.counter where each slot represents 1ms of time for the
%duration of imaging and the value in that slot represents the frame being
%collected at that time. e.g 1,1,1,1,1,...2,2,2,2,2,2,...3,3,3,3,3,3,...
counter_times_diff = diff(counterTimesMs);
counter=[];
IFI = mode(counter_times_diff); %duration of first frame must be estimated.
counter = repmat(counterValues(1),1,IFI); %beginning to create counter. using the estimated duration of first frame. Then concatenating on the duration of each frame after that in the forloop
for ii=1:length(counter_times_diff);
    counter = cat(2, counter, repmat(counterValues(ii+1),1,counter_times_diff(ii)));
end
end