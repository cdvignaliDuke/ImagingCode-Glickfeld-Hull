%Summary Statistic
colors = {'r', 'r', 'b', 'b', 'm', 'm', 'g', 'g', 'k', 'k', 'c', 'c', 'y', 'y', 'r', 'r', 'b', 'b'};
days = {'150718_img27', '150719_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '150518_img24', '150519_img24', '150518_img25', '150517_img25'};
DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryFolder\';
summary_succ = {}; 
summary_fail = {}; 
for kk = 1:length(days)
    curr_file_succ = strcat(DATA_DIR, days{kk}, '_success');
    summary_succ{kk} = load(curr_file_succ);
    
    curr_file_fail = strcat(DATA_DIR, days{kk}, '_fail');
    
    temp2 = load(curr_file_fail);
    summary_fail{kk} = temp2;
    
end 
%only plots ROIs in LS
summary_succ_mat = []
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{1}.success_roi(:,2:3,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{2}.success_roi(:,2:3,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{3}.success_roi(:,1,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{4}.success_roi(:,:,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{5}.success_roi(:,2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{6}.success_roi(:,2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{7}.success_roi(:,2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{8}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{9}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{10}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{11}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{12}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{13}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{14}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{15}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{16}.success_roi(:,2:3,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{17}.success_roi(:,1:2,7:9),1),2))')
summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{18}.success_roi(:,1:2,7:9),1),2))')

fail_succ_mat = []
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{1}.fail_roi(:,2:3,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{2}.fail_roi(:,2:3,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{3}.fail_roi(:,1,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{4}.fail_roi(:,:,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{5}.fail_roi(:,2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{6}.fail_roi(:,2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{7}.fail_roi(:,2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{8}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{9}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{10}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{11}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{12}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{13}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{14}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{15}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{16}.fail_roi(:,2:3,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{17}.fail_roi(:,1:2,7:9),1),2))')
fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{18}.fail_roi(:,1:2,7:9),1),2))')

%plots all lobules
% summary_succ_mat = []
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{1}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{2}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{3}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{4}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{5}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{6}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{7}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{8}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{9}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{10}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{11}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{12}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{13}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{14}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{15}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{16}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{17}.success_roi(:,:,7:9),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{18}.success_roi(:,:,7:9),1),2))')
% 
% fail_succ_mat = []
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{1}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{2}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{3}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{4}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{5}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{6}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{7}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{8}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{9}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{10}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{11}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{12}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{13}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{14}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{15}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{16}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{17}.fail_roi(:,:,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{18}.fail_roi(:,:,7:9),1),2))')

plot_succ_mat = mean(summary_succ_mat,2)
plot_fail_mat = mean(fail_succ_mat,2)
plot_succ_sm = std(summary_succ_mat,[],2)/sqrt(size(summary_succ_mat,2))
plot_fail_sm = std(fail_succ_mat,[],2)/sqrt(size(fail_succ_mat,2))

avg_summary_succ_mat = mean(mean(summary_succ_mat))
avg_summary_fail_mat = mean(mean(fail_succ_mat))

%%
%PLOT SCATTER 
figure;
for i = 1:length(days);
   if i > 14
       plot(plot_succ_mat(i), plot_fail_mat(i), ['x' colors{i}]); hold on;
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
title(['success vs fail summary']);

%%
% fail_avg = [];
% succ_avg = []; 
% fail_std = [];
% succ_std = []; 
% fail_sm = []; 
% succ_sm = []; 
% 
% for kk = 1:length(days)
% succ_avg(1,kk) = mean(mean(summary_succ{kk}.success_roi(:,7:9),1),2);
% fail_avg(1,kk) = mean(mean(summary_fail{kk}.fail_roi(:,7:9),1),2);
% end
% for kk = 1:length(days)
% succ_sm(1,kk) = std(mean(summary_succ{kk}.success_roi(:,7:9),2),[],1)./sqrt(size(summary_succ{kk}.success_roi,1));
% fail_sm(1,kk) = std(mean(summary_fail{kk}.fail_roi(:,7:9),2),[],1)./sqrt(size(summary_fail{kk}.fail_roi,1));
% end
