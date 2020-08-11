pre_event_time = 1000;
post_event_time = 4000;
preevent_frames = ceil(pre_event_time*(frame_rate/1000));
postevent_frames = ceil(post_event_time*(frame_rate/1000));

data_stim = nan(preevent_frames+postevent_frames,nCells,nTrials);
data_choice = nan(preevent_frames+postevent_frames,nCells,nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial)+29 < nframes
        data_stim(:,:,itrial) = npSub_tc(1+cStimOn(itrial)-preevent_frames:cStimOn(itrial)+postevent_frames,:);
    end
    if cDecision(itrial)
        if cDecision(itrial)+29 < nframes
            data_choice(:,:,itrial) = npSub_tc(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
        end
    end
end
dataf = mean(data_stim(1:preevent_frames,:,:),1);
data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
data_choice_dfof = bsxfun(@rdivide, bsxfun(@minus, data_choice, dataf), dataf);

data_stim_long = reshape(permute(data_stim_dfof,[1 3 2]), [(preevent_frames+postevent_frames)*nTrials nCells]);
data_choice_long = reshape(permute(data_choice_dfof,[1 3 2]), [(preevent_frames+postevent_frames)*nTrials nCells]);

rad_stim_long = reshape(rad_mat_start, [(preevent_frames+postevent_frames)*nTrials 1]);
rad_choice_long = reshape(rad_mat_decide, [(preevent_frames+postevent_frames)*nTrials 1]);

cen_stim_long = reshape(permute(centroid_mat_start,[1 3 2]), [(preevent_frames+postevent_frames)*nTrials 2]);
cen_choice_long = reshape(permute(centroid_mat_decide,[1 3 2]), [(preevent_frames+postevent_frames)*nTrials 2]);


ind_stim = find(isnan(rad_stim_long));
ind_choice = find(isnan(rad_choice_long));
rad_stim_temp = rad_stim_long;
rad_stim_temp(ind_stim,:) = [];
data_stim_temp = data_stim_long;
data_stim_temp(ind_stim,:) = [];
rad_choice_temp = rad_choice_long;
rad_choice_temp(ind_choice,:) = [];
data_choice_temp = data_choice_long;
data_choice_temp(ind_choice,:) = [];
r_stim_rad = zeros(1,nCells);
r_choice_rad = zeros(1,nCells);
for iCell = 1:nCells
    stim_temp = corrcoef(data_stim_temp(:,iCell), rad_stim_temp);
    r_stim_rad(1,iCell) = stim_temp(1,2);
    choice_temp = corrcoef(data_choice_temp(:,iCell), rad_choice_temp);
    r_choice_rad(1,iCell) = choice_temp(1,2);
end

cen_stim_temp = cen_stim_long;
cen_stim_temp(ind_stim,:) = [];
data_stim_temp = data_stim_long;
data_stim_temp(ind_stim,:) = [];
cen_choice_temp = cen_choice_long;
cen_choice_temp(ind_choice,:) = [];
data_choice_temp = data_choice_long;
data_choice_temp(ind_choice,:) = [];
r_stim_cen = zeros(2,nCells);
r_choice_cen = zeros(2,nCells);
for iCell = 1:nCells
    for i = 1:2
        stim_temp = corrcoef(data_stim_temp(:,iCell), cen_stim_temp(:,i));
        r_stim_cen(i,iCell) = stim_temp(1,2);
        choice_temp = corrcoef(data_choice_temp(:,iCell), cen_choice_temp(:,i));
        r_choice_cen(i,iCell) = choice_temp(1,2);
    end
end


n = floor(size(data_stim_temp,1)./2);
for iCell = 1:nCells
    stim_temp = corrcoef(data_stim_temp(1:n,iCell), rad_stim_temp(1:n,1));
    r_stim_rad_1(1,iCell) = stim_temp(1,2);
    stim_temp = corrcoef(data_stim_temp(n+1:end,iCell), rad_stim_temp(n+1:end,1));
    r_stim_rad_2(1,iCell) = stim_temp(1,2);
end


%% wheel
pre_frames = 100;
post_frames = 200;
tt = (1-pre_frames:post_frames).*frame_rate;
Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
qVals_frames_stim = nan(1+pre_frames+post_frames, uint16(length(input.trialOutcomeCell)));
for trN = 1:length(input.trialOutcomeCell)-1
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        qVals = double([input.quadratureValues{trN} input.quadratureValues{trN+1}]);
        cTimes = double([input.counterTimesUs{trN} input.counterTimesUs{trN+1}]./1000);
        cVals = double([input.counterValues{trN} input.counterValues{trN+1}]);
        stimVal = double(input.qStimOn{trN});
        ind_on = find(cVals == cStimOn(trN));
        cTime_min = cTimes(1+ind_on-pre_frames);
        cTime_max = cTimes(ind_on+post_frames);
        time_ind = find(qTimes>= cTime_min  & qTimes<=cTime_max);
        if length(time_ind)>2
            qTimes_sub = qTimes(time_ind);
            qVals_sub = qVals(time_ind)-stimVal;
            cTimes_sub = cTimes(ind_on-pre_frames:ind_on+post_frames);
            rep_ind = find(diff(qTimes_sub)==0);
            qTimes_sub(rep_ind) = [];
            qVals_sub(rep_ind) = [];
            qVals_frames_stim(:,trN) = interp1(qTimes_sub, qVals_sub, cTimes_sub)';
        else
            return
        end
    end
end

pre_frames = 100;
post_frames = 200;
tt = (1-pre_frames:post_frames).*frame_rate;
Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
qVals_frames_choice = nan(1+pre_frames+post_frames, uint16(length(input.trialOutcomeCell)));
for trN = 1:length(input.trialOutcomeCell)-1
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        qVals = double([input.quadratureValues{trN} input.quadratureValues{trN+1}]);
        cTimes = double([input.counterTimesUs{trN} input.counterTimesUs{trN+1}]./1000);
        cVals = double([input.counterValues{trN} input.counterValues{trN+1}]);
        stimVal = double(input.qStimOn{trN});
        ind_on = find(cVals == cDecision(trN));
        cTime_min = cTimes(1+ind_on-pre_frames);
        cTime_max = cTimes(ind_on+post_frames);
        time_ind = find(qTimes>= cTime_min  & qTimes<=cTime_max);
        if length(time_ind)>2
            qTimes_sub = qTimes(time_ind);
            qVals_sub = qVals(time_ind)-stimVal;
            cTimes_sub = cTimes(ind_on-pre_frames:ind_on+post_frames);
            rep_ind = find(diff(qTimes_sub)==0);
            qTimes_sub(rep_ind) = [];
            qVals_sub(rep_ind) = [];
            qVals_frames_choice(:,trN) = interp1(qTimes_sub, qVals_sub, cTimes_sub)';
        else
            return
        end
    end
end

qVals_stim = qVals_frames_stim(1+pre_frames-preevent_frames:pre_frames+postevent_frames,:);
qVals_choice = qVals_frames_choice(1+pre_frames-preevent_frames:pre_frames+postevent_frames,:);

qVals_stim_long = reshape(qVals_stim, [(preevent_frames+postevent_frames)*nTrials 1]);
qVals_choice_long = reshape(qVals_choice, [(preevent_frames+postevent_frames)*nTrials 1]);


ind_stim = find(isnan(qVals_stim_long));
ind_choice = find(isnan(qVals_choice_long));
q_stim_temp = qVals_stim_long;
q_stim_temp(ind_stim,:) = [];
data_stim_temp = data_stim_long;
data_stim_temp(ind_stim,:) = [];
q_choice_temp = qVals_choice_long;
q_choice_temp(ind_choice,:) = [];
data_choice_temp = data_choice_long;
data_choice_temp(ind_choice,:) = [];
r_stim_q = zeros(1,nCells);
r_choice_q = zeros(1,nCells);
for iCell = 1:nCells
    stim_temp = corrcoef(data_stim_temp(2:end,iCell), diff(q_stim_temp));
    r_stim_q(1,iCell) = stim_temp(1,2);
    choice_temp = corrcoef(data_choice_temp(2:end,iCell), diff(q_choice_temp));
    r_choice_q(1,iCell) = choice_temp(1,2);
end
