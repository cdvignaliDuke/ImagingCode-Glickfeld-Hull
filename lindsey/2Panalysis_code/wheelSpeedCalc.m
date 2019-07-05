function [wheel_speed] = wheelSpeedCalc(input,ticks,wheel);

%outputs wheel_speed in cm/s for each frame
%input is (1) input structure, (2) number of ticks on encoder (32, 64 or
%1000), (3) and type of wheel ('r' for red).

if wheel == 'r'
    diameter = 13.6;
    circ = pi.*diameter;
end
cm_per_tick = circ./(ticks.*4);
ntrials = length(input.tGratingDirectionDeg);
counterTimes = [];
counterValues = [];
wheelSpeedTimes = [];
wheelSpeedValues = [];
for itrial = 1:ntrials
    counterTimes = [counterTimes input.counterTimesUs{itrial}];
    counterValues = [counterValues input.counterValues{itrial}];
    wheelSpeedTimes = [wheelSpeedTimes input.wheelSpeedTimesUs{itrial}];
    wheelSpeedValues = [wheelSpeedValues input.wheelSpeedValues{itrial}];
end
nframes = input.counterValues{end}(end);
wheel_speed = nan(1,nframes);
for iframe = 1:nframes-1
    fr_time_start = counterTimes(find(counterValues==iframe,1,'last'));
    fr_time_end = counterTimes(find(counterValues==iframe+1));
    ind = intersect(find(wheelSpeedTimes>=fr_time_start),find(wheelSpeedTimes<=fr_time_end));
    if length(ind)>0
        wheel_speed(:,iframe) = mean(wheelSpeedValues(ind),2);
    else
        wheel_speed(:,iframe) = 0;
    end
end
wheel_speed(:,nframes) = wheel_speed(:,nframes-1);
frame_dur = mean(diff(counterTimes),2)./1000000;
wheel_speed = (wheel_speed.*cm_per_tick)./frame_dur;