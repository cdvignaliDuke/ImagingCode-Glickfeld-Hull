wheelSpeedVec=[];
wheelTimeVec = [];
counterVec = [];
counterTimeVec = [];
ind = find(input.counterValues{1}==1);
if length(ind)>1
    ind = ind(end);
end
for i = 1:ntrials
    wheelSpeedVec = [wheelSpeedVec input.wheelSpeedValues{i}];
    counterVec = [counterVec input.counterValues{i}];
    wheelTimeVec = [wheelTimeVec input.wheelSpeedTimesUs{i}-input.counterTimesUs{1}(ind)];
    counterTimeVec = [counterTimeVec input.counterTimesUs{i}-input.counterTimesUs{1}(ind)];
end
cStimOn = celleqel2mat_padded(input.cStimOn);
tWheelSpeed = cell(1,ntrials);
for i = 1:ntrials
    if cStimOn(i)+round(4000./frameRateHz) < counterVec(end)
        [val_s, ind_s] = min(abs(wheelTimeVec-counterTimeVec(cStimOn(i)-round(1000/frameRateHz))),[],2);
        [val_e, ind_e] = min(abs(wheelTimeVec-counterTimeVec(cStimOn(i)+round(4000/frameRateHz))),[],2);
        tWheelSpeed{:,i} = wheelSpeedVec(:,ind_s:ind_e);
        avgSpeed(1,i) = mean(tWheelSpeed{i},2);
        cvSpeed(1,i) = var(tWheelSpeed{i},[],2);
    end
end

figure;
ind_run = find(avgSpeed>10);
ind_norun = find(avgSpeed<=10);
subplot(2,2,1)
plot(tt, squeeze(nanmean(nanmean(temp_resp_tc(:,good_ind,1,ind_run),2),4)))
ylabel('Normalized dF/F')
hold on
subplot(2,2,2)
plot(tt, squeeze(nanmean(nanmean(temp_resp_tc(:,good_ind,1,ind_norun),2),4)))
ylabel('Normalized dF/F')
subplot(2,1,2)
errorbar(8000, nanmean(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,1,ind_run),1),2),4), nanstd(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,1,ind_run),1),4),[],2)./sqrt(length(good_ind)), 'or');
ylabel('Normalized dF/F')
hold on
errorbar(8000, nanmean(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,1,ind_norun),1),2),4), nanstd(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,1,ind_norun),1),4),[],2)./sqrt(length(good_ind)), 'ok');
for ioff = 1:noff
    ind_off_run = intersect(ind_run, find(tFramesOff == offs(ioff)));
    ind_off_norun = intersect(ind_norun, find(tFramesOff == offs(ioff)));
    subplot(2,2,1)
    plot(tt, squeeze(nanmean(nanmean(temp_resp_tc(:,good_ind,2,ind_off_run),2),4)))
    hold on
    title('Running')
    subplot(2,2,2)
    plot(tt, squeeze(nanmean(nanmean(temp_resp_tc(:,good_ind,2,ind_off_norun),2),4)))
    hold on
    title('Stationary')
    subplot(2,1,2)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,2,ind_off_run),1),2),4), nanstd(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,2,ind_off_run),1),4),[],2)./sqrt(length(good_ind)), 'or');
    hold on;
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,2,ind_off_norun),1),2),4), nanstd(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,2,ind_off_norun),1),4),[],2)./sqrt(length(good_ind)), 'ok');
end

norm_resp_tc_run = bsxfun(@rdivide, temp_resp_tc, nanmean(mean(temp_resp_tc(resp_win,:,1,ind_run),1),4));
norm_resp_tc_norun = bsxfun(@rdivide, temp_resp_tc, nanmean(mean(temp_resp_tc(resp_win,:,1,ind_norun),1),4));
figure;
subplot(2,2,1)
plot(tt, squeeze(nanmean(nanmean(norm_resp_tc_run(:,good_ind,1,ind_run),2),4)))
ylabel('Normalized dF/F')
hold on
subplot(2,2,2)
plot(tt, squeeze(nanmean(nanmean(norm_resp_tc_norun(:,good_ind,1,ind_norun),2),4)))
ylabel('Normalized dF/F')
subplot(2,1,2)
errorbar(8000, nanmean(nanmean(nanmean(norm_resp_tc_run(resp_win,good_ind,1,ind_run),1),2),4), nanstd(nanmean(nanmean(norm_resp_tc_run(resp_win,good_ind,1,ind_run),1),4),[],2)./sqrt(length(good_ind)), 'or');
ylabel('Normalized dF/F')
hold on
errorbar(8000, nanmean(nanmean(nanmean(norm_resp_tc_norun(resp_win,good_ind,1,ind_norun),1),2),4), nanstd(nanmean(nanmean(norm_resp_tc_norun(resp_win,good_ind,1,ind_norun),1),4),[],2)./sqrt(length(good_ind)), 'ok');
for ioff = 1:noff
    ind_off_run = intersect(ind_run, find(tFramesOff == offs(ioff)));
    ind_off_norun = intersect(ind_norun, find(tFramesOff == offs(ioff)));
    subplot(2,3,ioff)
    plot(tt, squeeze(nanmean(nanmean(norm_resp_tc_run(:,good_ind,2,ind_off_run),2),4)))
    hold on
    title('Running')
    subplot(2,2,2)
    plot(tt, squeeze(nanmean(nanmean(norm_resp_tc_norun(:,good_ind,2,ind_off_norun),2),4)))
    hold on
    title('Stationary')
    subplot(2,1,2)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(nanmean(nanmean(norm_resp_tc_run(resp_win,good_ind,2,ind_off_run),1),2),4), nanstd(nanmean(nanmean(norm_resp_tc_run(resp_win,good_ind,2,ind_off_run),1),4),[],2)./sqrt(length(good_ind)), 'or');
    hold on;
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(nanmean(nanmean(norm_resp_tc_norun(resp_win,good_ind,2,ind_off_norun),1),2),4), nanstd(nanmean(nanmean(norm_resp_tc_norun(resp_win,good_ind,2,ind_off_norun),1),4),[],2)./sqrt(length(good_ind)), 'ok');
end

