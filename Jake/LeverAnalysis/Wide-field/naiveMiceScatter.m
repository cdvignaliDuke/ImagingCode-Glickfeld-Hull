%Summary Statistic for naive animals 
colors = {'r', 'b'};
days = {'160515_img48', '160516_img47'}; %'150718_img27', '150719_img27', '150718_img27',
DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryFolder\';
%DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryNoShift\';
summary_succ = {}; 
summary_fail = {};
for kk = 1:length(days)
    curr_file_succ = strcat(DATA_DIR, days{kk}, '_success');
    summary_succ{kk} = load(curr_file_succ);
    curr_file_fail = strcat(DATA_DIR, days{kk}, '_fail');
    temp2 = load(curr_file_fail);
    summary_fail{kk} = temp2;
end

%working on selecting peak response
ROIcell = {[1:4], [1:2]};
summary_succ_mat = [];
summary_succ_mat_ROI = {};
for i = 1:length(days)
    avgResp = squeeze(mean(summary_succ{i}.success_roi));
    maxResp = [];
    summary_succ_mat_temp = [];
    for ii = ROIcell{i}
        maxResp = find(avgResp(ii,:) == max(avgResp(ii,7:10)));
        peakWindow = [(maxResp-1):(maxResp+1)];
    end
    summary_succ_mat_temp = squeeze(mean(summary_succ{i}.success_roi(:,ROIcell{i},peakWindow),1));
    if length(ROIcell{i})==1
        summary_succ_mat_temp = reshape(summary_succ_mat_temp,1,3); %if there is only one ROI then squeeze will automatically reshape summary_succ_mat
    end
    summary_succ_mat_temp_avg = mean(summary_succ_mat_temp, 1);
    summary_succ_mat_ROI = {summary_succ_mat_ROI, summary_succ_mat_temp};
    summary_succ_mat = cat(1, summary_succ_mat, mean(summary_succ_mat_temp_avg,1));
end
summary_succ_mat_ROI(1) = [];

summary_fail_mat = [];
summary_fail_mat_ROI = {};
for i = 1:length(days)
    avgResp = squeeze(mean(summary_fail{i}.fail_roi));
    maxResp = [];
    summary_fail_mat_temp = [];
    for ii = ROIcell{i}
        maxResp = find(avgResp(ii,:) == max(avgResp(ii,7:10)));
        peakWindow = [(maxResp-1):(maxResp+1)];
    end
    summary_fail_mat_temp = squeeze(mean(summary_fail{i}.fail_roi(:,ROIcell{i},peakWindow),1));
    if length(ROIcell{i})==1
        summary_fail_mat_temp = reshape(summary_fail_mat_temp,1,3); %if there is only one ROI then squeeze will automatically reshape summary_fail_mat
    end
    summary_fail_mat_temp_avg = mean(summary_fail_mat_temp, 1);
    summary_fail_mat_ROI = {summary_fail_mat_ROI, summary_fail_mat_temp};
    summary_fail_mat = cat(1, summary_fail_mat, mean(summary_fail_mat_temp_avg,1));
end
summary_fail_mat_ROI(1) = [];

plot_succ_mat = mean(summary_succ_mat,2);
plot_fail_mat = mean(summary_fail_mat,2);
plot_succ_sm = std(summary_succ_mat,[],2)/sqrt(size(summary_succ_mat,2));
plot_fail_sm = std(summary_fail_mat,[],2)/sqrt(size(summary_fail_mat,2));

avg_summary_succ_mat = mean(mean(summary_succ_mat));
avg_summary_fail_mat = mean(mean(summary_fail_mat));

%%
%PLOT SCATTER 
figure;
for i = 1:length(days);
   if i < 5
       plot(plot_succ_mat(i), plot_fail_mat(i), ['o' colors{i}], 'MarkerFaceColor', colors{i}); hold on;
   elseif i > 18
       plot(plot_succ_mat(i), plot_fail_mat(i), ['o' colors{i}]); hold on;
   else
       plot(plot_succ_mat(i), plot_fail_mat(i), ['o' colors{i}], 'MarkerFaceColor', colors{i}); hold on;
   end
end
legend(days{1:length(days)})
for i = 1:length(days);
    errorbarxy(plot_succ_mat(i)', plot_fail_mat(i)', plot_succ_sm(i)', plot_fail_sm(i)')%,...
    hold on   % 'Color',colors{i}); hold on; 
end
ylim([-.03 .25])
xlim([-.03 .25])
x = -.1:.1:1;
y = x;
hold on; plot(x,y,'k')
hline(0,'--k')
vline(0,'--k')
xlabel('df/f success condition');
ylabel('df/f fail condition');
title(['success vs fail summary Shift']);







