function baseline_timesMs = find_baseline_times(b_data);
%Determine baseline times: Time windows in which to take an F for df/f
%assumes lever solenoid is up during the iti

%define some useful variables
nTrials = length(b_data.trialOutcomeCell);
holdStartsMs = round(cell2mat(b_data.holdStartsMs));
totalReqHoldMs = cell2mat(b_data.tTotalReqHoldTimeMs);
baseline_times = zeros(2,nTrials);   %baseline_times will be matrix containing the times (beg and end) of each window during the iti

%if statements to determine the type of experiment conducted 
if b_data.doLever == 0;  % REWARD ONLY CONDITION or CUE REWARD PAIRING CONDITION. both use the same rule for establishing a baseline interval. 
    %if isfield(b_data,'doFakeMouseSuccessOnly') & b_data.doFakeMouseSuccessOnly==1 & b_data.postRewardMs >= 4000; %CUE REWARD PAIRING CONDITION
    for iT = 1:nTrials-1;
         baseline_times(1,iT) = holdStartsMs(iT)+totalReqHoldMs(iT)-900; %identifies a 700ms interval just before the appearance of the cue
         baseline_times(2,iT) = holdStartsMs(iT)+totalReqHoldMs(iT)-200;
    end
else
    %for lever trials I want to use a period just before leverHoldStart, not trialStart. 
    for iT = 1:nTrials-1
        baseline_times(1,iT) = holdStartsMs(iT)-900; %identifies a 700ms interval just before the appearance of the cue
        baseline_times(2,iT) = holdStartsMs(iT)-200;
    end
end 
baseline_timesMs = baseline_times; 
end 

