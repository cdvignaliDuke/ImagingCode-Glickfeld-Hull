
lickTimes = cell2mat_padded(input.lickometerTimesUs);
lickTimes(lickTimes==0) = [];
c_vals = cell2mat(input.counterValues);
c_times = cell2mat(input.counterTimesUs);
cTargetOn = cell2mat(input.cTargetOn);
lickFrames = zeros(1,length(lickTimes));
for this_lick = 1:length(lickFrames)
    lickFrames(1,this_lick) = find(c_times<lickTimes(this_lick),1,'last');
end

save(fullfile(jake_out,img_fn, [img_fn '_events.mat']), 'c_vals', 'lickTimes', 'c_times', 'cTargetOn', 'lickFrames');

pre_cue = round(1500/30);
post_cue = round(4000/30);
cueAligned = zeros(nTrials,length([-pre_cue:post_cue]));

for trial_num = 1:nTrials
    cueAligned(trial_num,:) = tc_avg(1,[(cTargetOn(trial_num)-pre_cue) : (cTargetOn(trial_num)+post_cue)]);
end

cueAlignedOrig = cueAligned;
for trial_num = 1:nTrials
    cueAligned(trial_num,:) = cueAligned(trial_num,:)-mean(cueAligned(trial_num, pre_cue-36:pre_cue-6)); %take df/f of cuealigned
end
cueAlignedMean = mean(cueAligned);
cueAlignedSEM = std(cueAligned,[],1)/sqrt(size(cueAligned,1));
x_axis = [-pre_cue:post_cue]*33;

figure;
shadedErrorBar(x_axis,cueAlignedMean, cueAlignedSEM);
ylabel('DF/F');
xlabel('TIME FROM CUE (MS)');
title('BACKGROUND ROI: IMG1030');
vline(600, 'g')


figure;
x_axis_full = [1:length(tc_avg)]*33;
plot(x_axis_full, tc_avg(1,:),'k');
vline([lickFrames]*33,'b');
ylabel('RAW F');
xlabel('TIME (MS)');
title('LICKING DURING  THE ITI');

format long g
figure;
plot(x_axis_full, tc_avg(1,:),'k');
ylabel('RAW F');
xlabel('TIME (MS)');
title('POCKEL CELL MODULATION');


