%trial over trial analysis
load([dest_sub '_evoked_events.mat'])
load([dest '_frame_times.mat'])
load([dest '_parse_behavior.mat'])
load([dest_sub '_pvals.mat'])

ifi = (frame_times(end)-frame_times(1))/length(frame_times);
pre_buffer = ceil(1000/double(ifi));
short_win = ceil(100/double(ifi));

%index trials
successIx = strcmp(input.trialOutcomeCell, 'success');
failIx = strcmp(input.trialOutcomeCell, 'failure');
pressIx = zeros(size(successIx));
pressIx(intersect(find(cell2mat(input.holdTimesMs)>=500),find(~isnan(trial_outcome.ind_press_prerelease)))) = 1;
reactTime = double(cell2mat(input.reactTimesMs));
holdTime = double(cell2mat(input.holdTimesMs));

%relate index back to success_tc framework
success_ind = cumsum(successIx);
fail_ind = cumsum(failIx);
press_ind = cumsum(pressIx);

%determine probability of success after each behavioral event
nextTrOutcome_success = [];
nextTrOutcome_fail = [];
nextTrOutcome_press = [];

for is = 1:size(success_tc,3)
    tr = find(success_ind == is, 1);
    if tr < size(successIx,2)
        nextTrOutcome_success = [nextTrOutcome_success successIx(1,tr+1)];
    end
end
for ie = 1:size(fail_tc,3)
    tr = find(fail_ind == ie, 1);
    if tr < size(successIx,2)
        nextTrOutcome_fail = [nextTrOutcome_fail successIx(1,tr+1)];
    end
end
for ip = 1:size(press_tc,3)
    tr = find(press_ind == ip, 1);
    if tr < size(successIx,2)
        nextTrOutcome_press = [nextTrOutcome_press successIx(1,tr+1)];
    end
end

secondTrOutcome_success = [];
secondTrOutcome_fail = [];
secondTrOutcome_press = [];

for is = 1:size(success_tc,3)
    tr = find(success_ind == is, 1);
    if tr < size(successIx,2)-1
        secondTrOutcome_success = [secondTrOutcome_success successIx(1,tr+2)];
    end
end
for ie = 1:size(fail_tc,3)
    tr = find(fail_ind == ie, 1);
    if tr < size(successIx,2)-1
        secondTrOutcome_fail = [secondTrOutcome_fail successIx(1,tr+2)];
    end
end
for ip = 1:size(press_tc,3)
    tr = find(press_ind == ip, 1);
    if tr < size(successIx,2)-1
        secondTrOutcome_press = [secondTrOutcome_press successIx(1,tr+2)];
    end
end

prevTrOutcome_success = [];
prevTrOutcome_fail = [];
prevTrOutcome_press = [];

for is = 1:size(success_tc,3)
    tr = find(success_ind == is, 1);
    if tr < size(successIx,2)-1
        prevTrOutcome_success = [prevTrOutcome_success successIx(1,tr+2)];
    end
end
for ie = 1:size(fail_tc,3)
    tr = find(fail_ind == ie, 1);
    if tr < size(successIx,2)-1
        prevTrOutcome_fail = [prevTrOutcome_fail successIx(1,tr+2)];
    end
end
for ip = 1:size(press_tc,3)
    tr = find(press_ind == ip, 1);
    if tr < size(successIx,2)-1
        prevTrOutcome_press = [prevTrOutcome_press successIx(1,tr+2)];
    end
end

[pctCorrect_all ci95_all] = binofit(sum(successIx,2),length(successIx));
[pctCorrect_success ci95_success] = binofit(sum(nextTrOutcome_success,2),length(nextTrOutcome_success));
[pctCorrect_fail ci95_fail] = binofit(sum(nextTrOutcome_fail,2),length(nextTrOutcome_fail));
[pctCorrect_press ci95_press] = binofit(sum(nextTrOutcome_press,2),length(nextTrOutcome_press));

save([dest_sub '_nextTrOutcome.mat'], 'nextTrOutcome_success', 'nextTrOutcome_fail', 'nextTrOutcome_press', 'successIx', 'failIx', 'pressIx', 'secondTrOutcome_success', 'secondTrOutcome_fail', 'secondTrOutcome_press', 'prevTrOutcome_success', 'prevTrOutcome_fail', 'prevTrOutcome_press');

%extract trial by trial outcomes according to CS
for ic = 1:size(success_tc,2)
    for iCS = 1:4
        if or(iCS ==1, iCS==3)
            event = 1;
        else
            event = 0;
        end
        if or(iCS ==1, iCS==2)
            win = round(200./double(ifi));
        else
            win = ceil(100./double(ifi));
        end
        success(ic).CS(iCS).event = event;
        fail(ic).CS(iCS).event = event;
        press(ic).CS(iCS).event = event;
        success(ic).CS(iCS).win = win;
        fail(ic).CS(iCS).win = win;
        press(ic).CS(iCS).win = win;
        success(ic).CS(iCS).NextTrOutcome = [];
        fail(ic).CS(iCS).NextTrOutcome = [];
        press(ic).CS(iCS).NextTrOutcome = [];
        press(ic).CS(iCS).ThisTrOutcome = [];
        success(ic).CS(iCS).SecondTrOutcome = [];
        fail(ic).CS(iCS).SecondTrOutcome = [];
        press(ic).CS(iCS).SecondTrOutcome = [];
        success(ic).CS(iCS).PrevTrOutcome = [];
        fail(ic).CS(iCS).PrevTrOutcome = [];
        press(ic).CS(iCS).PrevTrOutcome = [];
        success(ic).CS(iCS).RandTrOutcome = [];
        fail(ic).CS(iCS).RandTrOutcome = [];
        press(ic).CS(iCS).RandTrOutcome = [];
        success(ic).CS(iCS).NextTrReactTime = [];
        fail(ic).CS(iCS).NextTrReactTime = [];
        press(ic).CS(iCS).NextTrReactTime = [];
        success(ic).CS(iCS).NextTrHoldTime = [];
        fail(ic).CS(iCS).NextTrHoldTime = [];
        press(ic).CS(iCS).NextTrHoldTime = [];
    end
    for i = 1:2
        for is = 1:size(success_tc,3)
            tr = find(success_ind == is, 1);
            if tr < size(successIx,2)
                if i == 1
                    if success(ic).event(1,is)
                        iCS = 1;
                    else
                        iCS = 2;
                    end
                else
                    win = success(ic).CS(iCS).win;
                    if isnan(success(ic).event_ind(1,is))
                        iCS = 4;
                    elseif success(ic).event_ind(1,is)-pre_buffer < win
                        iCS = 3;
                    else
                        iCS = 4;
                    end
                end
                success(ic).CS(iCS).NextTrOutcome = [success(ic).CS(iCS).NextTrOutcome successIx(1,tr+1)];
                success(ic).CS(iCS).NextTrReactTime = [success(ic).CS(iCS).NextTrReactTime double(reactTime(1,tr+1))];
                success(ic).CS(iCS).NextTrHoldTime = [success(ic).CS(iCS).NextTrHoldTime double(holdTime(1,tr+1))];
                success(ic).CS(iCS).RandTrOutcome = [success(ic).CS(iCS).RandTrOutcome successIx(1,randsample(size(successIx,2),1))];
                if tr<size(successIx,2)-1
                    success(ic).CS(iCS).SecondTrOutcome = [success(ic).CS(iCS).SecondTrOutcome successIx(1,tr+2)];
                end
                if tr>1
                    success(ic).CS(iCS).PrevTrOutcome = [success(ic).CS(iCS).PrevTrOutcome successIx(1,tr-1)];
                end
            end
        end
        for ie = 1:size(fail_tc,3)
            tr = find(fail_ind == ie, 1);
            if tr < size(successIx,2)
                if i == 1
                    if fail(ic).event(1,ie)
                        iCS = 1;
                    else
                        iCS = 2;
                    end
                else
                    win = fail(ic).CS(iCS).win;
                    if isnan(fail(ic).event_ind(1,ie))
                        iCS = 4;
                    elseif fail(ic).event_ind(1,ie)-pre_buffer < win
                        iCS = 3;
                    else
                        iCS = 4;
                    end
                end
                fail(ic).CS(iCS).NextTrOutcome = [fail(ic).CS(iCS).NextTrOutcome successIx(1,tr+1)];
                fail(ic).CS(iCS).NextTrReactTime = [fail(ic).CS(iCS).NextTrReactTime reactTime(1,tr+1)];
                fail(ic).CS(iCS).NextTrHoldTime = [fail(ic).CS(iCS).NextTrHoldTime holdTime(1,tr+1)];
                fail(ic).CS(iCS).RandTrOutcome = [fail(ic).CS(iCS).RandTrOutcome successIx(1,randsample(size(successIx,2),1))];
                if tr<size(successIx,2)-1
                    fail(ic).CS(iCS).SecondTrOutcome = [fail(ic).CS(iCS).SecondTrOutcome successIx(1,tr+2)];
                end
                if tr>1
                    fail(ic).CS(iCS).PrevTrOutcome = [fail(ic).CS(iCS).PrevTrOutcome successIx(1,tr-1)];
                end
            end
        end
        for ip = 1:size(press_tc,3)
            tr = find(press_ind == ip, 1);
            if tr < size(successIx,2)
                if i == 1
                    if press(ic).event(1,ip)
                        iCS = 1;
                    else
                        iCS = 2;
                    end
                else
                    win = press(ic).CS(iCS).win;
                    if isnan(press(ic).event_ind(1,ip))
                        iCS = 4;
                    elseif (pre_buffer-press(ic).event_ind(1,ip) < win) & (pre_buffer-press(ic).event_ind(1,ip) >= 0)
                        iCS = 3;
                    else
                        iCS = 4;
                    end
                end
                press(ic).CS(iCS).ThisTrOutcome = [press(ic).CS(iCS).ThisTrOutcome successIx(1,tr)];
                press(ic).CS(iCS).NextTrOutcome = [press(ic).CS(iCS).NextTrOutcome successIx(1,tr+1)];
                press(ic).CS(iCS).NextTrReactTime = [press(ic).CS(iCS).NextTrReactTime reactTime(1,tr+1)];
                press(ic).CS(iCS).NextTrHoldTime = [press(ic).CS(iCS).NextTrHoldTime holdTime(1,tr+1)];
                press(ic).CS(iCS).RandTrOutcome = [press(ic).CS(iCS).RandTrOutcome successIx(1,randsample(size(successIx,2),1))];
                if tr<size(successIx,2)-1
                    press(ic).CS(iCS).SecondTrOutcome = [press(ic).CS(iCS).SecondTrOutcome successIx(1,tr+2)];
                end
                if tr>1
                    press(ic).CS(iCS).PrevTrOutcome = [press(ic).CS(iCS).PrevTrOutcome successIx(1,tr-1)];
                end
            end
        end
    end
    for iCS = 1:4
        if length(success(ic).CS(iCS).NextTrOutcome) > 20
            randTrials = randsample(size(success(ic).CS(iCS).NextTrOutcome,2),20);
            success(ic).CS(iCS).pctCorrectNext = binofit(sum(success(ic).CS(iCS).NextTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(success(ic).CS(iCS).SecondTrOutcome,2),20);
            success(ic).CS(iCS).pctCorrectSecond = binofit(sum(success(ic).CS(iCS).SecondTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(success(ic).CS(iCS).PrevTrOutcome,2),20);
            success(ic).CS(iCS).pctCorrectPrev = binofit(sum(success(ic).CS(iCS).PrevTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(success(ic).CS(iCS).RandTrOutcome,2),20);
            success(ic).CS(iCS).pctCorrectRand = binofit(sum(success(ic).CS(iCS).RandTrOutcome(1,randTrials),2),20);
            success(ic).CS(iCS).avgReact = mean(success(ic).CS(iCS).NextTrReactTime,2);
            randTrials = randsample(size(success(ic).CS(iCS).NextTrReactTime,2),20);
            success(ic).CS(iCS).semReact = std(success(ic).CS(iCS).NextTrReactTime(1,randTrials),[],2)./sqrt(20); 
            randTrials = randsample(size(success(ic).CS(iCS).NextTrHoldTime,2),20);
            success(ic).CS(iCS).avgHold = mean(success(ic).CS(iCS).NextTrHoldTime,2);
            success(ic).CS(iCS).semHold = std(success(ic).CS(iCS).NextTrHoldTime(1,randTrials),[],2)./sqrt(20);
        else
            success(ic).CS(iCS).pctCorrectNext = NaN;
            success(ic).CS(iCS).pctCorrectSecond = NaN;
            success(ic).CS(iCS).pctCorrectPrev = NaN;
            success(ic).CS(iCS).pctCorrectRand = NaN;
            success(ic).CS(iCS).avgReact = NaN;
            success(ic).CS(iCS).semReact = NaN;
            success(ic).CS(iCS).avgHold = NaN;
            success(ic).CS(iCS).semHold = NaN;
        end
        if length(fail(ic).CS(iCS).NextTrOutcome) > 20
            randTrials = randsample(size(fail(ic).CS(iCS).NextTrOutcome,2),20);
            fail(ic).CS(iCS).pctCorrectNext = binofit(sum(fail(ic).CS(iCS).NextTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(fail(ic).CS(iCS).SecondTrOutcome,2),20);
            fail(ic).CS(iCS).pctCorrectSecond = binofit(sum(fail(ic).CS(iCS).SecondTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(fail(ic).CS(iCS).PrevTrOutcome,2),20);
            fail(ic).CS(iCS).pctCorrectPrev = binofit(sum(fail(ic).CS(iCS).PrevTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(fail(ic).CS(iCS).RandTrOutcome,2),20);
            fail(ic).CS(iCS).pctCorrectRand = binofit(sum(fail(ic).CS(iCS).RandTrOutcome(1,randTrials),2),20);
            fail(ic).CS(iCS).avgReact = mean(fail(ic).CS(iCS).NextTrReactTime,2);
            randTrials = randsample(size(fail(ic).CS(iCS).NextTrReactTime,2),20);
            fail(ic).CS(iCS).semReact = std(fail(ic).CS(iCS).NextTrReactTime(1,randTrials),[],2)./sqrt(20); 
            randTrials = randsample(size(fail(ic).CS(iCS).NextTrHoldTime,2),20);
            fail(ic).CS(iCS).avgHold = mean(fail(ic).CS(iCS).NextTrHoldTime,2);
            fail(ic).CS(iCS).semHold = std(fail(ic).CS(iCS).NextTrHoldTime(1,randTrials),[],2)./sqrt(20);
        else
            fail(ic).CS(iCS).pctCorrectNext = NaN;
            fail(ic).CS(iCS).pctCorrectSecond = NaN;
            fail(ic).CS(iCS).pctCorrectPrev = NaN;
            fail(ic).CS(iCS).pctCorrectRand = NaN;
            fail(ic).CS(iCS).avgReact = NaN;
            fail(ic).CS(iCS).semReact = NaN;
            fail(ic).CS(iCS).avgHold = NaN;
            fail(ic).CS(iCS).semHold = NaN;
        end
        release(ic).CS(iCS).NextTrOutcome = [fail(ic).CS(iCS).NextTrOutcome success(ic).CS(iCS).NextTrOutcome];
        release(ic).CS(iCS).SecondTrOutcome = [fail(ic).CS(iCS).SecondTrOutcome success(ic).CS(iCS).SecondTrOutcome];
        release(ic).CS(iCS).PrevTrOutcome = [fail(ic).CS(iCS).PrevTrOutcome success(ic).CS(iCS).PrevTrOutcome];
        release(ic).CS(iCS).RandTrOutcome = [fail(ic).CS(iCS).RandTrOutcome success(ic).CS(iCS).RandTrOutcome];
        release(ic).CS(iCS).NextTrReactTime = [fail(ic).CS(iCS).NextTrReactTime success(ic).CS(iCS).NextTrReactTime];
        release(ic).CS(iCS).NextTrHoldTime = [fail(ic).CS(iCS).NextTrHoldTime success(ic).CS(iCS).NextTrHoldTime];
        if length(release(ic).CS(iCS).NextTrOutcome) > 20
            randTrials = randsample(size(release(ic).CS(iCS).NextTrOutcome,2),20);
            release(ic).CS(iCS).pctCorrectNext = binofit(sum(release(ic).CS(iCS).NextTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(release(ic).CS(iCS).SecondTrOutcome,2),20);
            release(ic).CS(iCS).pctCorrectSecond = binofit(sum(release(ic).CS(iCS).SecondTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(release(ic).CS(iCS).PrevTrOutcome,2),20);
            release(ic).CS(iCS).pctCorrectPrev = binofit(sum(release(ic).CS(iCS).PrevTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(release(ic).CS(iCS).RandTrOutcome,2),20);
            release(ic).CS(iCS).pctCorrectRand = binofit(sum(release(ic).CS(iCS).RandTrOutcome(1,randTrials),2),20);
            release(ic).CS(iCS).avgReact = mean(release(ic).CS(iCS).NextTrReactTime,2);
            randTrials = randsample(size(release(ic).CS(iCS).NextTrReactTime,2),20);
            release(ic).CS(iCS).semReact = std(release(ic).CS(iCS).NextTrReactTime(1,randTrials),[],2)./sqrt(20); 
            randTrials = randsample(size(release(ic).CS(iCS).NextTrHoldTime,2),20);
            release(ic).CS(iCS).avgHold = mean(release(ic).CS(iCS).NextTrHoldTime,2);
            release(ic).CS(iCS).semHold = std(release(ic).CS(iCS).NextTrHoldTime(1,randTrials),[],2)./sqrt(20);
        else
            release(ic).CS(iCS).pctCorrectNext = NaN;
            release(ic).CS(iCS).pctCorrectSecond = NaN;
            release(ic).CS(iCS).pctCorrectPrev = NaN;
            release(ic).CS(iCS).pctCorrectRand = NaN;
            release(ic).CS(iCS).avgReact = NaN;
            release(ic).CS(iCS).semReact = NaN;
            release(ic).CS(iCS).avgHold = NaN;
            release(ic).CS(iCS).semHold = NaN;
        end
        if length(press(ic).CS(iCS).NextTrOutcome) > 20
            randTrials = randsample(size(press(ic).CS(iCS).NextTrOutcome,2),20);
            press(ic).CS(iCS).pctCorrectNext = binofit(sum(press(ic).CS(iCS).NextTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(press(ic).CS(iCS).SecondTrOutcome,2),20);
            press(ic).CS(iCS).pctCorrectSecond = binofit(sum(press(ic).CS(iCS).SecondTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(press(ic).CS(iCS).PrevTrOutcome,2),20);
            press(ic).CS(iCS).pctCorrectPrev = binofit(sum(press(ic).CS(iCS).PrevTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(press(ic).CS(iCS).RandTrOutcome,2),20);
            press(ic).CS(iCS).pctCorrectRand = binofit(sum(press(ic).CS(iCS).RandTrOutcome(1,randTrials),2),20);
            randTrials = randsample(size(press(ic).CS(iCS).ThisTrOutcome,2),20);
            press(ic).CS(iCS).pctCorrectThis = binofit(sum(press(ic).CS(iCS).ThisTrOutcome(1,randTrials),2),20);
            press(ic).CS(iCS).avgReact = mean(press(ic).CS(iCS).NextTrReactTime,2);
            randTrials = randsample(size(press(ic).CS(iCS).NextTrReactTime,2),20);
            press(ic).CS(iCS).semReact = std(press(ic).CS(iCS).NextTrReactTime(1,randTrials),[],2)./sqrt(20); 
            randTrials = randsample(size(press(ic).CS(iCS).NextTrHoldTime,2),20);
            press(ic).CS(iCS).avgHold = mean(press(ic).CS(iCS).NextTrHoldTime,2);
            press(ic).CS(iCS).semHold = std(press(ic).CS(iCS).NextTrHoldTime(1,randTrials),[],2)./sqrt(20);
        else
            press(ic).CS(iCS).pctCorrectThis = NaN;
            press(ic).CS(iCS).pctCorrectNext = NaN;
            press(ic).CS(iCS).pctCorrectSecond = NaN;
            press(ic).CS(iCS).pctCorrectPrev = NaN;
            press(ic).CS(iCS).pctCorrectRand = NaN;
            press(ic).CS(iCS).avgReact = NaN;
            press(ic).CS(iCS).semReact = NaN;
            press(ic).CS(iCS).avgHold = NaN;
            press(ic).CS(iCS).semHold = NaN;
        end
    end
end

save([dest_sub '_event_outcome_summary.mat'], 'press', 'release', 'success', 'fail');
