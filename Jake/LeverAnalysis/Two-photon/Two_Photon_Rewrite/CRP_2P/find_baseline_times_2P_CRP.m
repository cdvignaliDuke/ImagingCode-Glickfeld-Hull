function [baseline_timesMs] = find_baseline_times_2P_CRP(b_data)
%Determine baseline times: Time windows in which to take an F for df/f

nTrials = length(b_data.trialOutcomeCell);
baseline_times = zeros(2,nTrials);   %baseline_times will be matrix containing the times (beg and end) of each window during the iti
%if session is a no-lever control...
if  b_data.doLever == 0 || (isfield(b_data, 'doAnalogLever') && b_data.doAnalogLever == 0)  %if no lever event occurs during that trial   i.e. REWARD ONLY CONDITION
    for iT = 1:nTrials;
        baseline_times(1,iT) = round(b_data.holdStartsMs{iT}+b_data.tTotalReqHoldTimeMs{iT}-1000)*1000;  %in reward only conditions the window for F is the last 700ms of the iti
        baseline_times(2,iT) = round(b_data.holdStartsMs{iT}+b_data.tTotalReqHoldTimeMs{iT})*1000;
    end
else
    disp('error: not a classical conditioning experiment');
    pause
end
ex_trials = length(find(isnan(baseline_times(1,:))));
baseline_timesMs = round(baseline_times/1000);
end