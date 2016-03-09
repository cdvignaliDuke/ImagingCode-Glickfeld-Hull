%Summary Statistic
colors = {'r', 'r', 'b', 'b', 'm', 'm', 'g', 'g', 'k', 'k', 'c', 'c', 'y', 'y', 'r', 'r', 'b', 'b'};
days = {'150718_img27', '150719_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '150518_img24', '150519_img24', '150518_img25', '150517_img25'};
DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryFolder\';
summary_fidget = {}; 
summary_fail = {};
for kk = 1:length(days)
    curr_file_fidget = strcat(DATA_DIR, days{kk}, '_fidget');
    summary_fidget{kk} = load(curr_file_fidget);
    curr_file_fail = strcat(DATA_DIR, days{kk}, '_fail');
    temp2 = load(curr_file_fail);
    summary_fail{kk} = temp2;
end 
%only plots ROIs in LS
summary_fidget_mat = []
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{1}.fidget_roi(:,2:3,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{2}.fidget_roi(:,2:3,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{3}.fidget_roi(:,1,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{4}.fidget_roi(:,:,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{5}.fidget_roi(:,2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{6}.fidget_roi(:,2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{7}.fidget_roi(:,2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{8}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{9}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{10}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{11}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{12}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{13}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{14}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{15}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{16}.fidget_roi(:,2:3,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{17}.fidget_roi(:,1:2,7:9),1),2))')
summary_fidget_mat = cat(1, summary_fidget_mat, squeeze(mean(mean(summary_fidget{18}.fidget_roi(:,1:2,7:9),1),2))')

summary_fail_mat = []
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{1}.fail_roi(:,2:3,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{2}.fail_roi(:,2:3,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{3}.fail_roi(:,1,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{4}.fail_roi(:,:,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{5}.fail_roi(:,2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{6}.fail_roi(:,2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{7}.fail_roi(:,2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{8}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{9}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{10}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{11}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{12}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{13}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{14}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{15}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{16}.fail_roi(:,2:3,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{17}.fail_roi(:,1:2,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{18}.fail_roi(:,1:2,7:9),1),2))')

plot_fidget_mat = mean(summary_fidget_mat,2)
plot_fail_mat = mean(summary_fail_mat,2)
plot_fidget_sm = std(summary_fidget_mat,[],2)/sqrt(size(summary_fidget_mat,2))
plot_fail_sm = std(summary_fail_mat,[],2)/sqrt(size(summary_fail_mat,2))

avg_summary_succ_mat = mean(mean(summary_fidget_mat))
avg_summary_fail_mat = mean(mean(summary_fail_mat))

%%
%PLOT SCATTER 
figure;
for i = 1:length(days);
   if i > 14
       plot(plot_fidget_mat(i), plot_fail_mat(i), ['x' colors{i}]); hold on;
   else
       plot(plot_fidget_mat(i), plot_fail_mat(i), ['o' colors{i}], 'MarkerFaceColor', colors{i}); hold on;
   end
end
legend(days{1:length(days)})
for i = 1:length(days);
    errorbarxy(plot_fidget_mat(i)', plot_fail_mat(i)', plot_fidget_sm(i)', plot_fail_sm(i)')%,...
    hold on   % 'Color',colors{i}); hold on; 
end
ylim([-.03 .25])
xlim([-.03 .25])
x = -.1:.1:1;
y = x;
hold on; plot(x,y,'k')
hline(0,'--k')
vline(0,'--k')
xlabel('df/f fidget condition');
ylabel('df/f fail condition');
title(['fidget vs fail summary']);

%% Plot no-lever control scatterplot


%Summary Statistic
colors = {'r', 'b', 'm', 'g'};
days = {'160208_img35', '160209_img36', '151222_img32', '151019_img30'};
DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryFolder\';
summary_succ = {}; 
summary_fail = {};
for kk = 1:length(days)
    curr_file_fidget = strcat(DATA_DIR, days{kk}, '_success');
    summary_succ{kk} = load(curr_file_fidget);
    curr_file_fail = strcat(DATA_DIR, days{kk}, '_fail');
    temp2 = load(curr_file_fail);
    summary_fail{kk} = temp2;
end 
%only plots ROIs in LS
summary_succ_mat = []
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{1}.success_roi(:,:,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{2}.success_roi(:,:,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{3}.success_roi(:,1:3,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{4}.success_roi(:,2:3,7:9),1),2))')

summary_fail_mat = []
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{1}.fail_roi(:,:,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{2}.fail_roi(:,:,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{3}.fail_roi(:,1:3,7:9),1),2))')
summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{4}.fail_roi(:,2:3,7:9),1),2))')

plot_fidget_mat = mean(summary_succ_mat,2)
plot_fail_mat = mean(summary_fail_mat,2)
plot_fidget_sm = std(summary_succ_mat,[],2)/sqrt(size(summary_succ_mat,2))
plot_fail_sm = std(summary_fail_mat,[],2)/sqrt(size(summary_fail_mat,2))

avg_summary_succ_mat = mean(mean(summary_succ_mat))
avg_summary_fail_mat = mean(mean(summary_fail_mat))

%%
%PLOT SCATTER 
figure;
for i = 1:length(days);
   if i > 14
       plot(plot_fidget_mat(i), plot_fail_mat(i), ['x' colors{i}]); hold on;
   else
       plot(plot_fidget_mat(i), plot_fail_mat(i), ['o' colors{i}], 'MarkerFaceColor', colors{i}); hold on;
   end
end
legend(days{1:length(days)})
for i = 1:length(days);
    errorbarxy(plot_fidget_mat(i)', plot_fail_mat(i)', plot_fidget_sm(i)', plot_fail_sm(i)')%,...
    hold on   % 'Color',colors{i}); hold on; 
end
ylim([-.03 .25])
xlim([-.03 .25])
x = -.1:.1:1;
y = x;
hold on; plot(x,y,'k')
hline(0,'--k')
vline(0,'--k')
xlabel('df/f reward condition');
ylabel('df/f no reward condition');
title(['No-Lever Control reward vs no reward']);

