function [baseline_timesMs, trial_outcome] = find_baseline_times_2P(b_data, trial_outcome, holdT_min)
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
    %identify windows of at least 400ms between a lever release and subsequent press. Only takes windows during the iti.
    %Sometimes, even when the solenoid is on, the first lever press can occur roughly 1ms before the iti ends and the trial begins.
    for iT = 1:nTrials-1
        ind_press = find(cell2mat(b_data.leverValues(iT))==1);      %finds locations of each lever press in the cell. Only really necessary for sessions without a solenoid
        ind_release = find(cell2mat(b_data.leverValues(iT))==0);
        leverTimes = cell2mat(b_data.leverTimesUs(iT));   %time of each lever event
        trialStart = cell2mat(b_data.holdStartsMs(iT)).*1000;
        if isempty(cell2mat(b_data.leverTimesUs(40:50)));
            leverTimes = trialStart;
        end
        ind_prerelease = leverTimes<=trialStart;   %misnomer. Lever events that occur during the iti. Again only relevant in sessions w/o solenoid
        if isempty(ind_prerelease);
            baseline_times(:,iT) = [NaN; NaN];
            trial_outcome.ind_press_prerelease(iT) = NaN;
        elseif sum(ind_prerelease,2) == 1 || sum(ind_prerelease,2) == 0
            trial_outcome.ind_press_prerelease(iT) = 1;
            if leverTimes(1)- (cell2mat(b_data.tThisTrialStartTimeMs(iT))*1000) > holdT_min + 100000  %Mainly relevant for non-solenoid trials. First lever event must occur 1sec or more after the start of the trials in order to collect a baseline time from that iti
                baseline_times(1,iT) = leverTimes(1) - holdT_min - 100000;
                baseline_times(2,iT) = leverTimes(1) - 100000;
                continue
            elseif baseline_times(:,iT) == 0
                baseline_times(:,iT) = [NaN; NaN];
            end
        else
            tag_prerelease = fliplr(find(ind_prerelease)); %again I think this is most relevant for non-solenoid trials as on solenoid trials find(ind_prerelease)=1
            trial_outcome.ind_press_prerelease(iT) = tag_prerelease(1);
            for i = 1:length(tag_prerelease)
                ip = tag_prerelease(i);
                if sum(ismember(ind_press,ip),2)
                    if ip == 1
                        if leverTimes(1)- (cell2mat(b_data.tThisTrialStartTimeMs(iT))*1000) > holdT_min + 100000;
                            baseline_times(1,iT) = leverTimes(1) - holdT_min - 100000;
                            baseline_times(2,iT) = leverTimes(1) - 100000;
                            break
                        elseif baseline_times(:,iT) == 0
                            baseline_times(:,iT) = [NaN; NaN];
                        end
                    elseif ip>1
                        if leverTimes(ip)- leverTimes(ip-1) > holdT_min + 100000
                            baseline_times(1,iT) = leverTimes(ip) - holdT_min - 100000;
                            baseline_times(2,iT) = leverTimes(ip) - 100000;
                            break
                        elseif leverTimes(1)- (cell2mat(b_data.tThisTrialStartTimeMs(iT))*1000) > holdT_min + 100000;
                            baseline_times(1,iT) = leverTimes(1) - holdT_min - 100000;
                            baseline_times(2,iT) = leverTimes(1) - 100000;
                            break
                        end
                    end
                end
            end
        end
    end
end
ex_trials = length(find(isnan(baseline_times(1,:))));

baseline_timesMs = round(baseline_times/1000);
end

