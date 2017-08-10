function [speedMS, speedTrStart,speedTarOn] = getWheelSpeedFromInput(expt_input)
    wv = expt_input.wheelSpeedValues; % pulses per second
    rotSpeedS = cellfun(@(x) x/32, wv, 'unif',0); % there are 32 pulses per rotation   
    speedMS = cell2mat(cellfun(@(x) x*0.36, rotSpeedS, 'unif',0)); % inside track of red wheel = 36cm.
    
    wt = expt_input.wheelSpeedTimesUs;
    trStart = expt_input.cFirstStim;
    target = expt_input.cTargetOn;
    cVals = expt_input.counterValues;
    cTimes = expt_input.counterTimesUs;
    
    start_counter_ind = cellfun(@(x,y) find(x == y),cVals,trStart,'unif',0);
    target_counter_ind = cellfun(@(x,y) find(x == y),cVals,target,'unif',0);
    
    start_time = cellfun(@(x,y) x(y), cTimes, start_counter_ind, 'unif',0);
    target_time = cellfun(@(x,y) x(y), cTimes, target_counter_ind, 'unif',0);
    
    [~, start_ind] = cellfun(@(x,y) min(abs(x - y)), wt, start_time, 'unif', 0);
    [~, target_ind] = cellfun(@(x,y) min(abs(x - y)), wt, target_time, 'unif', 0);
    
    trLengthCmlvSum = cumsum(cellfun(@length,wt));
    
    speedTrStart = cell2mat(start_ind) + [0 trLengthCmlvSum(2:end)];
    speedTarOn = cell2mat(target_ind) + [0 trLengthCmlvSum(2:end)];
end