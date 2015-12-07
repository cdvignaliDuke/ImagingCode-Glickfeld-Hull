%% Peak DF/F

nDay = size(elr_al,2);
x = 1:1:nDay;
baseline_avg = zeros(1,nDay);
peak_idx = zeros(1,nDay);
peak_val = zeros(1,nDay);
peak_amp = zeros(1,nDay);
frame_diff = zeros(1,nDay);
threshold_dfoverf = zeros(1,nDay);
threshold_idx = zeros(1,nDay);
threshold_val = zeros(1,nDay);

for i = 1:nDay
    baseline_avg(i) = mean(elr_v(10:40,i));
    [peak_val,peak_idx] = max(elr_v,[],1);
    peak_amp(i) = peak_val(i) - baseline_avg(i);
    frame_diff(i) = peak_idx(i) - 75;
end

for i = 1:nDay
    threshold_dfoverf(i) = baseline_avg(i) + 2*(std(elr_v(10:40,i)));
    threshold_idx(i) = find(elr_v(60:90,i) > threshold_dfoverf(i), 1);
end

%figure; scatter(x,peak_amp, 'filled', 'MarkerFaceColor', 'k');
figure; plot(threshold_idx, 'k');

%%
figure; hold on
scatter(x,peak_amp_AL, 'filled', 'MarkerFaceColor', 'b');
scatter(x,peak_amp_RL, 'filled', 'MarkerFaceColor', 'm');
scatter(x,peak_amp_LM, 'filled', 'MarkerFaceColor', 'c');
scatter(x,peak_amp_CS, 'filled', 'MarkerFaceColor', 'r');
scatter(x,peak_amp_CV, 'filled', 'MarkerFaceColor', 'g');

figure; hold on
scatter(x,threshold_idx_al, 'filled', 'MarkerFaceColor', 'b');
scatter(x,threshold_idx_rl, 'filled', 'MarkerFaceColor', 'm');
scatter(x,threshold_idx_lm, 'filled', 'MarkerFaceColor', 'c');
scatter(x,threshold_idx_cs, 'filled', 'MarkerFaceColor', 'r');
scatter(x,threshold_idx_cv, 'filled', 'MarkerFaceColor', 'g');
scatter(x,threshold_idx_v, 'filled', 'MarkerFaceColor', 'k');


figure; hold on
plot(x,peak_amp_al, '-b', 'linewidth', 1);
plot(x,peak_amp_rl, '-k', 'linewidth', 1);
plot(x,peak_amp_lm, '-c', 'linewidth', 1);
plot(x,peak_amp_cs, '-r', 'linewidth', 1);
plot(x,peak_amp_cv, '-g', 'linewidth', 1);
plot(x,peak_amp_v, '-m', 'linewidth', 1);

figure; hold on
plot(x,threshold_idx_al, '-b', 'linewidth', 1);
plot(x,threshold_idx_rl, '-k', 'linewidth', 1);
plot(x,threshold_idx_lm, '-c', 'linewidth', 1);
plot(x,threshold_idx_cs, '-r', 'linewidth', 1);
plot(x,threshold_idx_cv, '-g', 'linewidth', 1);
plot(x,threshold_idx_v, '-m', 'linewidth', 1);





