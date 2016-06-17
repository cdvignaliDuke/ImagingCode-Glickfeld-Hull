%% Behavioral criteria
clear
fidgetMax = 0.20;
corrMin = 0.50;
lapseMax = 0.10;
daysUnfiltered = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
ROIcell = {[1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:6], [2,3,5], [1], [1:2]};
ANALYSIS_DIR ='Z:\Analysis\LeverAnalysis\';
curr_cd = cd; 
days = {};
for kk= 1:length(daysUnfiltered);
    destySucc = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', daysUnfiltered{kk}, '_success');
    destyFail = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', daysUnfiltered{kk}, '_fail');
    destyFidget = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', daysUnfiltered{kk}, '_fidget');
    destyTooFast = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', daysUnfiltered{kk}, '_tooFast');
    destyLapse = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', daysUnfiltered{kk}, '_lapse');
    load(destySucc);
    load(destyFail);
    load(destyFidget);
    load(destyTooFast);
    if exist(strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', daysUnfiltered{kk}, '_lapse.mat'))==2;
        load(destyLapse);
        lapse = size(lapse_roi,1);
    else
        lapse=0;
    end
    corr = size(success_roi,1);
    early = size(fail_roi,1);
    fidget = size(fidget_roi,1);
    tooFast = size(tooFast_roi,1);
    total = corr + early + fidget + tooFast + lapse; 
    if fidget/total < fidgetMax
        if corr > early
            if lapse/total < lapseMax
                days = [days, daysUnfiltered{kk}];
            end
        end
    end
end
days

%Summary Statistic
colors = {'r', 'c', 'b', 'm', 'r', 'b', 'm', 'g', 'k', 'c', 'y',  'r', 'r', 'b', 'b', 'r', 'b', 'm', 'g', 'k', 'c', 'y',};
%days = {'150518_img24', '150519_img24', '150518_img25', '150517_img25', '150716_img27', '150718_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%ROIcell = {[1:2], [2:3], [1:2], [1:2], [1,2], [1,2,4], [1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:6], [2,3,5], [4], [1:2]};
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
for kk = 1:length(daysUnfiltered)
    temp = strcat('i', daysUnfiltered{kk});
    days_roi_matcher.(temp)= ROIcell{kk};
end

%working on selecting peak response
summary_succ_mat = [];
summary_succ_mat_ROI = {};
for i = 1:length(days)
    avgResp = squeeze(mean(summary_succ{i}.success_roi));
    maxResp = [];
    summary_succ_mat_temp = [];
    for ii = days_roi_matcher.(strcat('i', days{i}))
        maxResp = find(avgResp(ii,:) == max(avgResp(ii,7:10)));
        peakWindow = [(maxResp-1):(maxResp+1)];
    end
    summary_succ_mat_temp = squeeze(mean(summary_succ{i}.success_roi(:,days_roi_matcher.(strcat('i', days{i})),peakWindow),1));
    if length(days_roi_matcher.(strcat('i', days{i})))==1
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
    for ii = days_roi_matcher.(strcat('i', days{i}))
        maxResp = find(avgResp(ii,:) == max(avgResp(ii,7:10)));
        peakWindow = [(maxResp-1):(maxResp+1)];
    end
    summary_fail_mat_temp = squeeze(mean(summary_fail{i}.fail_roi(:,days_roi_matcher.(strcat('i', days{i})),peakWindow),1));
    if length(days_roi_matcher.(strcat('i', days{i})))==1
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
       plot(plot_succ_mat(i), plot_fail_mat(i), ['x' colors{i}]); hold on;
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

%% FOR NON LS AREAS
% 
% colors = {'r', 'r', 'b', 'b', 'm', 'm', 'g', 'g', 'k', 'k', 'y', 'c', 'c'};
% days = {'150718_img27', '150719_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '150519_img24', '150518_img25', '150517_img25'};
% aa = [1,4];
% %only plots ROIs NOT in LS
% summary_succ_mat = []
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{1}.success_roi(:,1,7:9),1),2))') %150718_img27
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{2}.success_roi(:,1,7:9),1),2))') %150719_img27
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{3}.success_roi(:,2,7:9),1),2))') %150716_img28
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{4}.success_roi(:,3,7:9),1),2))') %150717_img28
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{5}.success_roi(:,1,7:9),1),2))') %151021_img29
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{6}.success_roi(:,1,7:9),1),2))') %151022_img29
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{7}.success_roi(:,1,7:9),1),2))') %151009_img30
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{8}.success_roi(:,3:4,7:9),1),2))') %151011_img30
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{9}.success_roi(:,3:4,7:9),1),2))') %32
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{10}.success_roi(:,4:5,7:9),1),2))') %32
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{11}.success_roi(:,aa,7:9),1),2))')%24
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{12}.success_roi(:,3,7:9),1),2))') %25
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{13}.success_roi(:,3:4,7:9),1),2))')
% 
% fail_succ_mat = []
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{1}.fail_roi(:,1,7:9),1),2))') %27
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{2}.fail_roi(:,1,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{3}.fail_roi(:,2,7:9),1),2))') %28
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{4}.fail_roi(:,3,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{5}.fail_roi(:,1,7:9),1),2))') %29
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{6}.fail_roi(:,1,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{7}.fail_roi(:,1,7:9),1),2))') %30
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{8}.fail_roi(:,3:4,7:9),1),2))')
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{9}.fail_roi(:,3:4,7:9),1),2))') %32
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{10}.fail_roi(:,4:5,7:9),1),2))') 
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{11}.fail_roi(:,aa,7:9),1),2))') %24
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{12}.fail_roi(:,3,7:9),1),2))') %25
% fail_succ_mat = cat(1, fail_succ_mat, squeeze(mean(mean(summary_fail{13}.fail_roi(:,3:4,7:9),1),2))')


%%

% %Summary Statistic separeting rand=1000 vs rand=4500
% colors = {'r', 'r', 'r', 'r','r', 'b', 'b', 'b', 'b', 'b', 'b', 'b'};
% %days = {'150718_img27', '150719_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '150518_img24', '150519_img24', '150518_img25', '150517_img25'};
% days = {'160131_img36', '160131_img35', '151212_img32', '160315_img38', '160320_img41', '160129_img36', '160129_img35', '151009_img30', '151011_img30', '151211_img32', '160314_img38', '160319_img41'};
% 
% DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryFolder\';
% %DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryNoShift\';
% summary_succ = {}; 
% summary_fail = {};
% for kk = 1:length(days)
%     curr_file_succ = strcat(DATA_DIR, days{kk}, '_success');
%     summary_succ{kk} = load(curr_file_succ);
%     curr_file_fail = strcat(DATA_DIR, days{kk}, '_fail');
%     temp2 = load(curr_file_fail);
%     summary_fail{kk} = temp2;
% end 
% 
% %OneThousand = [10, 12, 14, '160315_img38', '160320_img41'];
% summary_succ_mat = []
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{1}.success_roi(:,1:2,6:8),1),2))') %10
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{2}.success_roi(:,1:2,6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{3}.success_roi(:,1:2,6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{4}.success_roi(:,1:4,6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{5}.success_roi(:,1:4,6:8),1),2))')
% %FortyFiveHundred = [7, 8, 9, 11, 13, '160314_img38', '160319_img41'];
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{6}.success_roi(:,2,6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{7}.success_roi(:,1:2,6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{8}.success_roi(:,1:2,6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{9}.success_roi(:,1:2,6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{10}.success_roi(:,1:2,6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{11}.success_roi(:,[1,2,3,5],6:8),1),2))')
% summary_succ_mat = cat(1, summary_succ_mat, squeeze(mean(mean(summary_succ{12}.success_roi(:,1:2,6:8),1),2))')
% 
% 
% %OneThousand = [10, 12, 14, '160315_img38', '160320_img41'];
% summary_fail_mat =[];
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{1}.fail_roi(:,1:2,6:8),1),2))') %10
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{2}.fail_roi(:,1:2,6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{3}.fail_roi(:,1:2,6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{4}.fail_roi(:,1:4,6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{5}.fail_roi(:,1:4,6:8),1),2))')
% %FortyFiveHundred = [7, 8, 9, 11, 13, '160314_img38', '160319_img41'];
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{6}.fail_roi(:,2,6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{7}.fail_roi(:,1:2,6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{8}.fail_roi(:,1:2,6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{9}.fail_roi(:,1:2,6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{10}.fail_roi(:,1:2,6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{11}.fail_roi(:,[1,2,3,5],6:8),1),2))')
% summary_fail_mat = cat(1, summary_fail_mat, squeeze(mean(mean(summary_fail{12}.fail_roi(:,1:2,6:8),1),2))')
% 
% plot_succ_mat = mean(summary_succ_mat,2)
% plot_fail_mat = mean(summary_fail_mat,2)
% plot_succ_sm = std(summary_succ_mat,[],2)/sqrt(size(summary_succ_mat,2))
% plot_fail_sm = std(summary_fail_mat,[],2)/sqrt(size(summary_fail_mat,2))
% 
% avg_summary_succ_mat = mean(mean(summary_succ_mat))
% avg_summary_fail_mat = mean(mean(summary_fail_mat))
% 
% %%
% %PLOT SCATTER 
% figure;
% for i = 1:length(days);
%     plot(plot_succ_mat(i), plot_fail_mat(i), ['o' colors{i}], 'MarkerFaceColor', colors{i}); hold on;
% end
% legend(days{1:length(days)})
% for i = 1:length(days);
%     errorbarxy(plot_succ_mat(i)', plot_fail_mat(i)', plot_succ_sm(i)', plot_fail_sm(i)')%,...
%     hold on   % 'Color',colors{i}); hold on; 
% end
% ylim([-.03 .25])
% xlim([-.03 .25])
% x = -.1:.1:1;
% y = x;
% hold on; plot(x,y,'k')
% hline(0,'--k')
% vline(0,'--k')
% xlabel('df/f success condition');
% ylabel('df/f fail condition');
% title(['success vs fail summary No Shift']);
% 

