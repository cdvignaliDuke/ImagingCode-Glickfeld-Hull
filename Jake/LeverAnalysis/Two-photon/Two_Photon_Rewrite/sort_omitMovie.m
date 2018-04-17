% omitBi_movie = []; omitLate_movie = [];
% omitEarly_movie = [];
load('_release_movies.mat'); load('_lick_stats.mat'); load('parse_behavior.mat');
load('_pvals.mat'); load('_cell_TCs.mat');
cc = struct(); cc.omitBi_movie = [];
cc.omitLate_movie = []; cc.omitEarly_movie =[];
cc.omitBi_movie_lick= []; cc.omitEarly_movie_lick = [];
cc.omitLate_movie_lick = []; cc.omitEarly_cell = [];
cc.omitLate_cell = []; cc.omitBi_cell = []; cc.prevBi = []; cc.prevEarly = []; cc.prevLate = [];
cc.omitBi_lever = []; cc.omitEarly_lever = []; 
cc.omitNeg_movie = [];
cc.omitLate_lever = []; cc.omitNeg_lick = [];
cc.omitNeg_cell = [];
zero_idx = round(size(omitReward_movie,3)/2);

[StartPress_duration, StopPress_duration]=analogLeverAna('Z:\home\andrew\Behavior\Data\data-i065-180110-1903.mat');
% 'Z:\home\andrew\Behavior\Data\data-i060-180104-1938.mat'
%'Z:\home\andrew\Behavior\Data\data-i063-180109-1847.mat'
% 'Z:\home\andrew\Behavior\Data\data-i065-180110-1903.mat'

omitRewardIndx = find(trial_outcome.omitRewardIndx);
prev_trial = [];

for j = 1:size(omitReward_movie,1)
    if trial_outcome.hasReward(omitRewardIndx(j) - 1) == 1
        omitR_prevR_TC = squeeze(omitReward_movie(j, :, :));
    else
        omitR_prevNR_TC = squeeze(omitReward_movie(j, :, :));
    end
end

fig = figure;
errorbar(tt2,mean(avg_OR_long(find(release_h),:),1),std(avg_OR_long(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'r')

for i = 1:size(omitReward_movie,2) % cell
    for j = 1:size(omitReward_movie,1) % trial
        
        v = squeeze(omitReward_movie(j,i,:));
        [maxtab, mintab] = peakdet(v, 0.01);
        max_peak = maxtab(maxtab(:,2) == max(maxtab(:,2)),:);
        max_peak_idx = max_peak(1,1); 
        max_peak = max_peak(1,2);
        temp = maxtab(maxtab(:,2) < max(maxtab(:,2)),:);
        second_peak = max(temp(:,2));
        second_peak_idx = maxtab(maxtab(:,2) == second_peak,1);
        
        third_temp =  maxtab (maxtab(:,2) < second_peak,:);
        third_peak = max(third_temp(:,2));
        third_temp_idx = maxtab(maxtab(:,2) == third_peak,1);

        if third_peak > max_peak*0.4 && (abs(second_peak_idx - max_peak_idx) <= 3 || (second_peak_idx - zero_idx) < -3)
            second_peak = third_peak;
            second_peak_idx = maxtab(maxtab(:,2) == second_peak,1);
        end
        % find bi phasic resp event
        if second_peak > max_peak*0.1 && ...
                abs(second_peak_idx - max_peak_idx) > 3 && ...
                 (second_peak_idx - zero_idx) > -3 && max_peak_idx - zero_idx > 10 %for 065
            
            cc.omitBi_movie = [cc.omitBi_movie; v'];
            cc.omitBi_movie_lick = [cc.omitBi_movie_lick; omitR_lick(j,:)];
            cc.omitBi_lever = [cc.omitBi_lever; StopPress_duration(j,:)];
            cc.omitBi_cell = [cc.omitBi_cell i];
            cc.prevBi = [cc.prevBi trial_outcome.hasReward(omitRewardIndx(j) - 1)];
        elseif max_peak_idx - zero_idx > 10
            
            cc.omitLate_movie = [cc.omitLate_movie; v'];
            cc.omitLate_movie_lick = [cc.omitLate_movie_lick; omitR_lick(j,:)];
            cc.omitLate_lever = [cc.omitLate_lever; StopPress_duration(j,:)];
            cc.omitLate_cell = [cc.omitLate_cell i];
            cc.prevLate = [cc.prevLate trial_outcome.hasReward(omitRewardIndx(j) - 1)];
        elseif max_peak_idx - zero_idx <= 10 && max_peak_idx - zero_idx > -15 && v(zero_idx + 10) < max_peak*0.5
            
            cc.omitEarly_movie = [cc.omitEarly_movie; v'];
            cc.omitEarly_movie_lick = [cc.omitEarly_movie_lick; omitR_lick(j,:)];
            cc.omitEarly_lever = [cc.omitEarly_lever; StopPress_duration(j,:)];
            cc.omitEarly_cell = [cc.omitEarly_cell i];
            cc.prevEarly = [cc.prevEarly trial_outcome.hasReward(omitRewardIndx(j) - 1)];
        end
        
        % find negative resp ~300 to 800 ms after release
        min_peak = mintab(mintab(:,2) == min(mintab(:,2)),:);
        min_peak_idx = min_peak(1,1); 
        min_peak = min_peak(1,2); 
        if min_peak_idx - zero_idx > 9 && min_peak < mean(v(1:(zero_idx-3)))
            cc.omitNeg_movie = [cc.omitNeg_movie; v'];
            cc.omitNeg_lick = [cc.omitNeg_lick; omitR_lick(j,:)];
            cc.omitNeg_cell = [cc.omitNeg_cell i];
        end
    end
end
maxy = 0.3; miny = -0.1;
fig = figure; p = panel(); p.pack(2,1);
p(1,1).select()
title(['Negative Response', ' nTimes = ', num2str(size(cc.omitNeg_movie,1))]);
errorbar((-15:15)*33, mean(cc.omitNeg_movie,1), std(cc.omitNeg_movie,1)/sqrt(size(cc.omitNeg_movie,1)));
p(2,1).select()
bar((-15:15)*33, mean(cc.omitNeg_lick,1)); ylim([miny*1.1 maxy*1.1]); xlim([-500 500])
ylim([0 0.5]); xlim([-500 500]);
title(['nCells = ', num2str(size(unique(cc.omitNeg_cell),2)), ' out of total Cells = ', num2str(size(omitReward_movie,2))]);
saveas(fig, ['_omitReward_negative.fig']);
print(['_omitReward_negative.eps'], '-depsc');
print(['_omitReward_negative.pdf'], '-dpdf');

bi_early_cell = intersect(unique(cc.omitBi_cell), unique(cc.omitEarly_cell))
biEarly_percent = length(bi_early_cell)/size(omitReward_movie,2)

late_early_cell = intersect(unique(cc.omitLate_cell), unique(cc.omitEarly_cell))
lateEarly_percent = length(late_early_cell)/size(omitReward_movie,2)

bi_prevReward_percent = sum(cc.prevBi) / length(cc.prevBi)*100;
early_prevReward_percent = sum(cc.prevEarly) / length(cc.prevEarly)*100;
late_prevReward_percent = sum(cc.prevLate) / length(cc.prevLate)*100;

% maxy = max([max(max(cc.omitBi_movie)) max(max(cc.omitEarly_movie)) max(max(cc.omitLate_movie))]);
% miny = min([min(min(cc.omitBi_movie)) min(min(cc.omitEarly_movie)) min(min(cc.omitLate_movie))]);
maxy = 0.3; miny = -0.1;
fig = figure; p = panel();
p.pack(3,3); p.de.margin = 0;

p(1,1).select()
title(['Bi phase', ' nTimes = ', num2str(size(cc.omitBi_movie,1)), 'previous trial is rewarded % = ', num2str(bi_prevReward_percent)]);
errorbar((-15:15)*33, mean(cc.omitBi_movie,1), std(cc.omitBi_movie,1)/sqrt(size(cc.omitBi_movie,1)));
ylim([miny*1.1 maxy*1.1]); xlim([-600 600])
% plot commands
p(2,1).select()
plot(-600:1:600, mean(cc.omitBi_lever,1));
ylim([-5 100]); xlim([-600 600])
p(3,1).select()
bar((-15:15)*33, mean(cc.omitBi_movie_lick,1));
ylim([0 0.5]); xlim([-600 600])

p(1,2).select()
title(['Early', ' nTimes = ', num2str(size(cc.omitEarly_movie,1)), 'previous trial is rewarded % = ', num2str(early_prevReward_percent)]);
errorbar((-15:15)*33, mean(cc.omitEarly_movie,1), std(cc.omitEarly_movie,1)/sqrt(size(cc.omitEarly_movie,1)));
ylim([miny*1.1 maxy*1.1]); xlim([-600 600])
% plot commands
p(2,2).select()
plot(-600:1:600, mean(cc.omitEarly_lever,1));
ylim([-5 100]); xlim([-600 600])
p(3,2).select()
bar((-15:15)*33, mean(cc.omitEarly_movie_lick,1));
ylim([0 0.5]); xlim([-600 600])

p(1,3).select()
title(['Late', ' nTimes = ', num2str(size(cc.omitLate_movie,1)), 'previous trial is rewarded % = ', num2str(late_prevReward_percent)]);
errorbar((-15:15)*33, mean(cc.omitLate_movie,1), std(cc.omitLate_movie,1)/sqrt(size(cc.omitLate_movie,1)));
ylim([miny*1.1 maxy*1.1]); xlim([-600 600])
% plot commands
p(2,3).select()
plot(-600:1:600, mean(cc.omitLate_lever,1));
ylim([-5 100]); xlim([-600 600])
p(3,3).select()
bar((-15:15)*33, mean(cc.omitLate_movie_lick,1));
ylim([0 0.5]); xlim([-600 600])

xlabel(['% both bi-phase and early = ', num2str(biEarly_percent*100), ' % both late and early = ', num2str(lateEarly_percent*100), '%']);
saveas(fig, ['_Early_Late_lickRate.fig']);
print(['_Early_Late_lickRate.eps'], '-depsc');
print(['_Early_Late_lickRate.pdf'], '-dpdf');

% plot commands
subplot(2,2,1);errorbar((-15:15)*33, mean(cc.omitBi_movie,1), std(cc.omitBi_movie,1)/sqrt(size(cc.omitBi_movie,1)));
hold on; bar((-15:15)*33, mean(cc.omitBi_movie_lick,1)/5); ylim([miny*1.1 maxy*1.1]); xlim([-600 600])
title(['Bi phase', ' nTrials = ', num2str(size(cc.omitBi_movie,1))]);
subplot(2,2,2);errorbar((-15:15)*33, mean(cc.omitLate_movie,1), std(cc.omitLate_movie,1)/sqrt(size(cc.omitLate_movie,1)));
hold on; bar((-15:15)*33, mean(cc.omitLate_movie_lick,1)/5); ylim([miny*1.1 maxy*1.1]); xlim([-600 600])
title('Late response')
subplot(2,2,3);errorbar((-15:15)*33, mean(cc.omitEarly_movie,1), std(cc.omitEarly_movie,1)/sqrt(size(cc.omitEarly_movie,1)));
hold on; bar((-15:15)*33, mean(cc.omitEarly_movie_lick,1)/5); ylim([miny*1.1 maxy*1.1]); xlim([-600 600])
title('Early response')
supertitle('dFF and lick rate for different type of response')

saveas(fig, ['_Early_Late_lickRate.fig']);
print(['_Early_Late_lickRate.eps'], '-depsc');
print(['_Early_Late_lickRate.pdf'], '-dpdf');
