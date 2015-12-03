%% plot pupil dynamics
%normalize to max area (for whole experiment)
rad_max = sqrt(max(Area_temp,[],1)/(pi))*calib;
rad_mat_target_norm = rad_mat_target./rad_max;
rad_mat_down_norm = rad_mat_down./rad_max;

b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
b2Ix = find(cell2mat(input.tBlock2TrialNumber));
reactTimes = cell2mat(input.reactTimesMs);
figure;
min_x = min(nanmean(rad_mat_target_norm(1:15,:),1),[],2)*0.9;
max_x = max(nanmean(rad_mat_target_norm(1:15,:),1),[],2)*1.1;
subplot(2,3,3)
scatter(nanmean(rad_mat_target_norm(1:15,b1Ix),1), reactTimes(:,b1Ix),'ok');
hold on
successIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
scatter(nanmean(rad_mat_target_norm(1:15,b2Ix),1), reactTimes(:,b2Ix),'og');
ylim([0 600])
xlim([min_x max_x])
title('Pupil size before target')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
subplot(2,3,1)
scatter(nanmean(rad_mat_down_norm(1:15,b1Ix),1), reactTimes(:,b1Ix),'ok');
hold on
scatter(nanmean(rad_mat_down_norm(1:15,b2Ix),1), reactTimes(:,b2Ix),'og');
ylim([0 600])
xlim([min_x max_x])
title('Pupil size before press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
subplot(2,3,2)
scatter(nanmean(rad_mat_down_norm(16:30,b1Ix),1), reactTimes(:,b1Ix),'ok');
hold on
scatter(nanmean(rad_mat_down_norm(16:30,b2Ix),1), reactTimes(:,b2Ix),'og');
ylim([0 600])
xlim([min_x max_x])
title('Pupil size after press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')

success1Ix = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'success')));
success2Ix = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
[n, edges] = histcounts([nanmean(rad_mat_target_norm(1:15,success1Ix),1) nanmean(rad_mat_target_norm(1:15,success2Ix),1)],8);
[n_b1, bin_b1] = histc(nanmean(rad_mat_target_norm(1:15,success1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_target_norm(1:15,success2Ix),1),edges);
subplot(2,3,6)
for i = 1:length(edges)
    ind = find(bin_b1 == i-1);
    errorbarxy(nanmean(nanmean(rad_mat_target_norm(1:15,success1Ix(ind)),1),2), mean(reactTimes(:,success1Ix(ind)),2), std(nanmean(rad_mat_target_norm(1:15,success1Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success1Ix(ind)),[],2)./sqrt(length(ind)), {'ok', 'k', 'k'});
    hold on;
    ind = find(bin_b2 == i-1);
    errorbarxy(nanmean(nanmean(rad_mat_target_norm(1:15,success2Ix(ind)),1),2), mean(reactTimes(:,success2Ix(ind)),2), std(nanmean(rad_mat_target_norm(1:15,success2Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success2Ix(ind)),[],2)./sqrt(length(ind)), {'og', 'g', 'g'});
    hold on;
end
ylim([0 600])
xlim([min_x max_x])
ylabel('React time (ms)')
xlabel('Normalized pupil size')
title('Average before target')
[n, edges] = histcounts([nanmean(rad_mat_down_norm(1:15,success1Ix),1) nanmean(rad_mat_down_norm(1:15,success2Ix),1)],8);
[n_b1, bin_b1] = histc(nanmean(rad_mat_down_norm(1:15,success1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_down_norm(1:15,success2Ix),1),edges);
subplot(2,3,4)
for i = 1:length(edges)
    ind = find(bin_b1 == i-1);
    errorbarxy(nanmean(nanmean(rad_mat_down_norm(1:15,success1Ix(ind)),1),2), mean(reactTimes(:,success1Ix(ind)),2), std(nanmean(rad_mat_down_norm(1:15,success1Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success1Ix(ind)),[],2)./sqrt(length(ind)), {'ok', 'k', 'k'})
    hold on;
    ind = find(bin_b2 == i-1);
    errorbarxy(nanmean(nanmean(rad_mat_down_norm(1:15,success2Ix(ind)),1),2), mean(reactTimes(:,success2Ix(ind)),2), std(nanmean(rad_mat_down_norm(1:15,success2Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success2Ix(ind)),[],2)./sqrt(length(ind)), {'og', 'g', 'g'})
    hold on;
end
ylim([0 600])
xlim([min_x max_x])
title('Average before press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
[n, edges] = histcounts([nanmean(rad_mat_down_norm(16:30,success1Ix),1) nanmean(rad_mat_down_norm(16:30,success2Ix),1)],8);
[n_b1, bin_b1] = histc(nanmean(rad_mat_down_norm(16:30,success1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_down_norm(16:30,success2Ix),1),edges);
subplot(2,3,5)
for i = 1:length(edges)
    ind = find(bin_b1 == i-1);
    errorbarxy(nanmean(nanmean(rad_mat_down_norm(16:30,success1Ix(ind)),1),2), mean(reactTimes(:,success1Ix(ind)),2), std(nanmean(rad_mat_down_norm(16:30,success1Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success1Ix(ind)),[],2)./sqrt(length(ind)), {'ok', 'k', 'k'})
    hold on;
    ind = find(bin_b2 == i-1);
    errorbarxy(nanmean(nanmean(rad_mat_down_norm(16:30,success2Ix(ind)),1),2), mean(reactTimes(:,success2Ix(ind)),2), std(nanmean(rad_mat_down_norm(16:30,success2Ix(ind)),1),[],2)./sqrt(length(ind)), std(reactTimes(:,success2Ix(ind)),[],2)./sqrt(length(ind)), {'og', 'g', 'g'})
    hold on;
end
ylim([0 600])
xlim([min_x max_x])
title('Average after press')
ylabel('React time (ms)')
xlabel('Normalized pupil size')
suptitle([mouse ' ' date ' Pupil size and react time'])
print([fnout '_PupilvsReact.pdf'], '-dpdf');    

%hit rate vs pupil size- before press
[n, edges] = histcounts([nanmean(rad_mat_down_norm(1:15,:),1)],5);
smIx = strcmp(input.trialOutcomeCell,'success') + strcmp(input.trialOutcomeCell,'ignore');
successIx = strcmp(input.trialOutcomeCell,'success');
missedIx = strcmp(input.trialOutcomeCell,'ignore');
sm1Ix = intersect(find(~isnan(rad_mat_down_norm(1,:))), intersect(b1Ix, find(smIx)));
sm2Ix = intersect(find(~isnan(rad_mat_down_norm(1,:))), intersect(b2Ix, find(smIx)));
[n_b1, bin_b1] = histc(nanmean(rad_mat_down_norm(1:15,sm1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_down_norm(1:15,sm2Ix),1),edges);
figure;
subplot(3,2,1)
for i = 1:length(edges)
    ind1 = find(bin_b1 == i);
    if (sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2))>5
        [pct1, err1] = binofit(sum(successIx(sm1Ix(ind1)),2), sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        errorbarxy(nanmean(nanmean(rad_mat_down_norm(1:15,sm1Ix(ind1)),1),2), pct1, std(nanmean(rad_mat_down_norm(1:15,sm1Ix(ind1)),1),[],2)./sqrt(length(ind1)), pct1-err1(1), {'ok', 'k', 'k'})
        hold on;
    end
    ind2 = find(bin_b2 == i);
    if sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2)>5
        [pct2, err2] = binofit(sum(successIx(sm2Ix(ind2)),2), sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        errorbarxy(nanmean(nanmean(rad_mat_down_norm(1:15,sm2Ix(ind2)),1),2), pct2, std(nanmean(rad_mat_down_norm(1:15,sm2Ix(ind2)),1),[],2)./sqrt(length(ind2)), pct2-err2(1), {'og', 'g', 'g'})
        hold on
    end
end
xlim([0 1])
ylim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')

tGratingDirection = cell2mat(input.tGratingDirectionDeg);
Oris = unique(tGratingDirection);
nOri = length(Oris);
M = nOri+1;
R = fliplr(linspace(0,1,M)) .';
colmat = flipud(horzcat(R, zeros(size(R)), zeros(size(R))));
ori_rad_mat = NaN(nOri,length(edges));
ori_hit_mat = NaN(nOri,length(edges));
ori_n_mat = NaN(nOri,length(edges));
for iOri = 2:nOri
    ori_ind = find(tGratingDirection(sm1Ix) == Oris(iOri));
    for  i = 1:length(edges)
        ind1 = intersect(ori_ind, find(bin_b1 == i));
        if length(ind1>1)
            ori_rad_mat(iOri,i) = nanmean(nanmean(rad_mat_down_norm(1:15,sm1Ix(ind1)),1),2);
            ori_hit_mat(iOri,i) = sum(successIx(sm1Ix(ind1)),2)/(sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
            x = (sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
            if isempty(x)
                x = 0;
            end
            ori_n_mat(iOri,i) = x;
        else
            ori_n_mat(iOri,i) = NaN;
        end
    end
    subplot(3,2,3)
    plot(ori_rad_mat(iOri,:), ori_hit_mat(iOri,:), '-', 'Color', colmat(iOri,:), 'Marker', 'o')
    hold on
    subplot(3,2,4)
    plot(ori_rad_mat(iOri,:), ori_n_mat(iOri,:), '-',  'Color', colmat(iOri,:), 'Marker', 'o')
    hold on;
end
subplot(3,2,3)
ylim([0 1])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')
title('Visual')
legend(num2str(chop(Oris(2:end)',2)))
subplot(3,2,4)
ylim([0 max(max(ori_n_mat,[],1),[],2)])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Number of trials')

tVolume = chop(double(celleqel2mat_padded(input.tSoundTargetAmplitude)),2);
Vols = unique(tVolume);
nVol = length(Vols);
M = nVol+1;
R = fliplr(linspace(0,1,M)) .';
colmat = flipud(horzcat(R, zeros(size(R)), zeros(size(R))));
vol_rad_mat = zeros(nVol,length(edges));
vol_hit_mat = zeros(nVol,length(edges));
vol_n_mat = zeros(nVol,length(edges));
for iVol = 2:nVol
    vol_ind = find(tVolume(sm2Ix) == Vols(iVol));
    for  i = 1:length(edges)
        ind2 = intersect(vol_ind, find(bin_b2 == i));
        vol_rad_mat(iVol,i) = nanmean(nanmean(rad_mat_down_norm(1:15,sm2Ix(ind2)),1),2);
        vol_hit_mat(iVol,i) = sum(successIx(sm2Ix(ind2)),2)/(sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        x = (sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        if isempty(x)
            x = 0;
        end
        vol_n_mat(iVol,i) = x;
    end
    subplot(3,2,5)
    plot(vol_rad_mat(iVol,:), vol_hit_mat(iVol,:), '-', 'Color', colmat(iVol,:), 'Marker', 'o')
    hold on
    subplot(3,2,6)
    plot(vol_rad_mat(iVol,:), vol_n_mat(iVol,:), '-', 'Color', colmat(iVol,:), 'Marker', 'o')
    hold on;
end
subplot(3,2,5)
ylim([0 1])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')
title('Auditory')
legend(num2str(chop(Vols(2:end)',2)))
subplot(3,2,6)
ylim([0 max(max(vol_n_mat,[],1),[],2)])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Number of trials')


subplot(3,2,2)
for iVol = 2:nVol
    vol_ind = find(tVolume == Vols(iVol));
    plot(Vols(iVol),sum(successIx(vol_ind),2)/(sum(successIx(vol_ind),2) + sum(missedIx(vol_ind),2)), 'og');
    hold on
end
for iOri = 2:nOri
    ori_ind = find(tGratingDirection == Oris(iOri));
    plot(Oris(iOri)./Oris(end),sum(successIx(ori_ind),2)/(sum(successIx(ori_ind),2) + sum(missedIx(ori_ind),2)), 'ok');
    hold on
end
ylim([0 1])
xlim([0 1])
xlabel('Change (%)')
ylabel('Hit rate')
suptitle([mouse ' ' date ' Pupil size (before press) and Hit rate'])
print([fnout '_prepress_pupilvHR.pdf'], '-dpdf');

%hit rate vs pupil size- before target
[n, edges] = histcounts([nanmean(rad_mat_target_norm(1:15,:),1)],5);
smIx = strcmp(input.trialOutcomeCell,'success') + strcmp(input.trialOutcomeCell,'ignore');
successIx = strcmp(input.trialOutcomeCell,'success');
missedIx = strcmp(input.trialOutcomeCell,'ignore');
sm1Ix = intersect(find(~isnan(rad_mat_target_norm(1,:))), intersect(b1Ix, find(smIx)));
sm2Ix = intersect(find(~isnan(rad_mat_target_norm(1,:))), intersect(b2Ix, find(smIx)));
[n_b1, bin_b1] = histc(nanmean(rad_mat_target_norm(1:15,sm1Ix),1),edges);
[n_b2, bin_b2] = histc(nanmean(rad_mat_target_norm(1:15,sm2Ix),1),edges);
figure;
subplot(3,2,1)
for i = 1:length(edges)
    ind1 = find(bin_b1 == i);
    if (sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2))>5
        [pct1, err1] = binofit(sum(successIx(sm1Ix(ind1)),2), sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        errorbarxy(nanmean(nanmean(rad_mat_target_norm(1:15,sm1Ix(ind1)),1),2), pct1, std(nanmean(rad_mat_target_norm(1:15,sm1Ix(ind1)),1),[],2)./sqrt(length(ind1)), pct1-err1(1), {'ok', 'k', 'k'})
        hold on;
    end
    ind2 = find(bin_b2 == i);
    if sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2)>5
        [pct2, err2] = binofit(sum(successIx(sm2Ix(ind2)),2), sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        errorbarxy(nanmean(nanmean(rad_mat_target_norm(1:15,sm2Ix(ind2)),1),2), pct2, std(nanmean(rad_mat_target_norm(1:15,sm2Ix(ind2)),1),[],2)./sqrt(length(ind2)), pct2-err2(1), {'og', 'g', 'g'})
        hold on
    end
end
xlim([0 1])
ylim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')

tGratingDirection = cell2mat(input.tGratingDirectionDeg);
Oris = unique(tGratingDirection);
nOri = length(Oris);
M = nOri+1;
R = fliplr(linspace(0,1,M)) .';
colmat = flipud(horzcat(R, zeros(size(R)), zeros(size(R))));
ori_rad_mat = zeros(nOri,length(edges));
ori_hit_mat = zeros(nOri,length(edges));
ori_n_mat = zeros(nOri,length(edges));
for iOri = 2:nOri
    ori_ind = find(tGratingDirection(sm1Ix) == Oris(iOri));
    for  i = 1:length(edges)
        ind1 = intersect(ori_ind, find(bin_b1 == i));
        ori_rad_mat(iOri,i) = nanmean(nanmean(rad_mat_target_norm(1:15,sm1Ix(ind1)),1),2);
        ori_hit_mat(iOri,i) = sum(successIx(sm1Ix(ind1)),2)/(sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        x = (sum(successIx(sm1Ix(ind1)),2) + sum(missedIx(sm1Ix(ind1)),2));
        if isempty(x)
            x = 0;
        end
        ori_n_mat(iOri,i) = x;
    end
    subplot(3,2,3)
    plot(ori_rad_mat(iOri,:), ori_hit_mat(iOri,:),'-', 'Color', colmat(iOri,:), 'Marker', 'o')
    hold on
    subplot(3,2,4)
    plot(ori_rad_mat(iOri,:), ori_n_mat(iOri,:),'-', 'Color', colmat(iOri,:), 'Marker', 'o')
    hold on;
end
subplot(3,2,3)
ylim([0 1])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')
title('Visual')
legend(num2str(chop(Oris(2:end)',2)))
subplot(3,2,4)
ylim([0 max(max(ori_n_mat,[],1),[],2)])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Number of trials')

M = nVol+1;
R = fliplr(linspace(0,1,M)) .';
colmat = flipud(horzcat(R, zeros(size(R)), zeros(size(R))));
vol_rad_mat = zeros(nVol,length(edges));
vol_hit_mat = zeros(nVol,length(edges));
vol_n_mat = zeros(nVol,length(edges));
for iVol = 2:nVol
    vol_ind = find(tVolume(sm2Ix) == Vols(iVol));
    for  i = 1:length(edges)
        ind2 = intersect(vol_ind, find(bin_b2 == i));
        vol_rad_mat(iVol,i) = nanmean(nanmean(rad_mat_target_norm(1:15,sm2Ix(ind2)),1),2);
        vol_hit_mat(iVol,i) = sum(successIx(sm2Ix(ind2)),2)/(sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        x = (sum(successIx(sm2Ix(ind2)),2) + sum(missedIx(sm2Ix(ind2)),2));
        if isempty(x)
            x = 0;
        end
        vol_n_mat(iVol,i) = x;
    end
    subplot(3,2,5)
    plot(vol_rad_mat(iVol,:), vol_hit_mat(iVol,:),'-', 'Color', colmat(iVol,:), 'Marker', 'o')
    hold on
    subplot(3,2,6)
    plot(vol_rad_mat(iVol,:), vol_n_mat(iVol,:),'-', 'Color', colmat(iVol,:), 'Marker', 'o')
    hold on;
end
subplot(3,2,5)
ylim([0 1])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Hit rate')
title('Auditory')
legend(num2str(chop(Vols(2:end)',2)))
subplot(3,2,6)
ylim([0 max(max(vol_n_mat,[],1),[],2)])
xlim([0 1])
xlabel('Normalized pupil radius')
ylabel('Number of trials')

subplot(3,2,2)
for iVol = 2:nVol
    vol_ind = find(tVolume == Vols(iVol));
    plot(Vols(iVol),sum(successIx(vol_ind),2)/(sum(successIx(vol_ind),2) + sum(missedIx(vol_ind),2)), 'og');
    hold on
end
for iOri = 2:nOri
    ori_ind = find(tGratingDirection == Oris(iOri));
    plot(Oris(iOri)./Oris(end),sum(successIx(ori_ind),2)/(sum(successIx(ori_ind),2) + sum(missedIx(ori_ind),2)), 'ok');
    hold on
end
ylim([0 1])
xlim([0 1])
xlabel('Change (%)')
ylabel('Hit rate')
suptitle([mouse ' ' date ' Pupil size (before target) and Hit rate'])
print([fnout '_pretarget_pupilvHR.pdf'], '-dpdf');


%distribution of pupil size
rad_all = (sqrt(Area_temp./pi))*calib;
cleverup = cell2mat(input.cLeverUp);
cleverdown = cell2mat(input.cLeverDown);
tr_frames = [];
for itrial = 1:ntrials
    tr_frames = [tr_frames double([cleverdown(itrial):cleverup(itrial)])];
end
min_x = min(rad_all,[],1)*0.9;
max_x = max(rad_all,[],1)*1.1;
figure;
subplot(2,2,1)
hist(rad_all,100)
title('All Frames')
xlim([min_x max_x])
subplot(2,2,2)
hist(rad_all(tr_frames,:),50)
title('Trial Frames')
iti_frames = find(ismember(1:size(rad_all,1), tr_frames)==0);
xlim([min_x max_x])
subplot(2,2,3)
hist(rad_all(iti_frames,:),100)
title('ITI Frames')
xlim([min_x max_x])
subplot(2,2,4)
hist(rad_all(tr_frames,:)./max(rad_all,[],1),50)
xlim([0 1])
title('Trial Frames- norm to max pupil size')
print([fnout '_pupil_size_hist.pdf'], '-dpdf');