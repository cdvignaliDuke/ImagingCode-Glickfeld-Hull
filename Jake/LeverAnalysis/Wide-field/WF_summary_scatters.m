%% bx criteria, days, ROIs, and path
clear
fidgetMax = 1.0;
corrMin = 0.60;
lapseMax = 1;
daysUnfiltered = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46'}; %'150718_img27', '150719_img27',
ROIcell = {[1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:4], [2,3,5], [1], [1:2], [3:4]};  %LS areas
ANALYSIS_DIR ='Z:\Analysis\WF Lever Analysis\';
curr_cd = cd;
days = {};

%% This forloop loads individual datasets, checks to see if they fit the bx criteria and edits the sessions included accordingly
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
    %check to see if lapse_roi exists. If there is only one lapse then it will be a different size than the others. 
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
    %Bx criteria for inclusion
    if fidget/total < fidgetMax
        if corr/total >= corrMin
            if lapse/total < lapseMax
                days = [days, daysUnfiltered{kk}];
            end
        end
    end
end
days  %report the days that meet bx criteria

% Non-LS areas
%daysUnfiltered = {'151009_img30', '151011_img30', '151212_img32', '160129_img36', '160314_img38', '160315_img38', '160606_img46'};
%days = {'151009_img30', '151011_img30', '151212_img32', '160129_img36', '160314_img38', '160315_img38', '160606_img46'};
%ROIcell = {[4,5], [4,5], [4,5], [3], [5], [4], [5]}; %for Non-LS areas
%colors = {'b', 'b', 'b', 'r', 'r', 'b', 'r'}; %blue=Vermis red=CrusI

%no-lever controls
days = {'160209_img36', '151222_img32', '151019_img30', '160725_img53'}; 
daysUnfiltered = {'160209_img36', '151222_img32', '151019_img30', '160725_img53'}; 
colors = {'b', 'g', 'm', kk'};
ROIcell = {[2], [3], [2:3], [1:2]};


%% (2)SUMMARY STATISTIC
%loads all the datasets again and stores them in summary_succ and summary_fail
%days = {'150518_img24', '150519_img24', '150518_img25', '150517_img25', '150716_img27', '150718_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%ROIcell = {[1:2], [2:3], [1:2], [1:2], [1,2], [1,2,4], [1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:6], [2,3,5], [4], [1:2]};

% UNCOMMENT AFTER RUNNING NOLEVER CONTROLS -------- colors = {'r', 'c', 'b', 'm', 'r', 'b', 'm', 'g', 'k', 'c', 'y',  'r', 'r', 'b', 'b', 'r', 'b', 'm', 'g', 'k', 'c', 'y',};
DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryFolder\';
summary_succ = {}; 
summary_fail = {};
summary_cue  = {};
for kk = 1:length(days)   %probably a much more efficient way to do this
    curr_file_succ = strcat(DATA_DIR, days{kk}, '_success');
    summary_succ{kk} = load(curr_file_succ);
    curr_file_fail = strcat(DATA_DIR, days{kk}, '_fail');
    summary_fail{kk} = load(curr_file_fail);
    curr_file_cue = strcat(DATA_DIR, days{kk}, '_cue');
    summary_cue{kk} = load(curr_file_cue);
end
for kk = 1:length(daysUnfiltered)
    temp = strcat('i', daysUnfiltered{kk});
    days_roi_matcher.(temp)= ROIcell{kk};
end

%shift all TCs of individual trials so they are baselined to 3 frame
for i = 1:length(days);
    for ii=1:length(summary_succ); %should be the same for succ and fail
        for iii = 1:size(summary_succ{i}.success_roi,2);
            for iiii = 1:size(summary_succ{i}.success_roi,1);
                shift = summary_succ{i}.success_roi(iiii,iii,[1:3]); 
                summary_succ{i}.success_roi(iiii,iii,:) = summary_succ{i}.success_roi(iiii,iii,:)-shift;
            end
            for iiii = 1:size(summary_fail{i}.fail_roi,1);
                shift = summary_fail{i}.fail_roi(iiii,iii,[1:3]);
                summary_fail{i}.fail_roi(iiii,iii,:) = summary_fail{i}.fail_roi(iiii,iii,:)-shift;
            end
            for iiii = 1:size(summary_cue{i}.cue_roi,1);
                shift = summary_cue{i}.cue_roi(iiii,iii,[1:3]);
                summary_cue{i}.cue_roi(iiii,iii,:) = summary_cue{i}.cue_roi(iiii,iii,:)-shift;
            end
        end
    end
end

%% FIND PEAK TIME AND VALUES USING LINEAR INTERPOLATION

XVector = [1:.1:size(success_roi,3)];
succ_time_to_peak = [];
fail_time_to_peak = [];
cue_time_to_peak = [];
peak_value_succ = [];
peak_value_fail = [];
peak_value_cue = [];
temp_lat_peak_succ = [];
temp_lat_peak_fail = [];
temp_lat_peak_cue = [];
tbyt_sm_lat_peak_fail = [];
tbyt_sm_lat_peak_succ = [];
tbyt_sm_lat_peak_cue = [];
tbyt_lat_peak_succ = [];
tbyt_lat_peak_fail = [];
tbyt_lat_peak_cue = [];
plot_succ_sm = [];
plot_fail_sm = [];
plot_cue_sm = [];
tbyt_peak_val_succ = [];
tbyt_peak_val_succ_sm = [];
tbyt_peak_val_cue_sm = [];
tbyt_peak_val_fail = [];
tbyt_peak_val_cue = [];
tbyt_peak_val_fail_sm = [];
tbyt_peak_val_cue_sm = [];

%figure;  %uncomment to plot all average TCs on a single figure. 
time_before = 5;
time_after = 5;
for kk = 1:length(days)  %for each animal..
    %does interpolation here
    success_roi = summary_succ{kk}.success_roi(:,[days_roi_matcher.(strcat('i', days{kk}))], :); %collects only the ROIs of interest. dims: 1=trial# 2=ROI# 3=Time
    success_roi_interp = nan(length(XVector), size(success_roi,1), size(success_roi,2)); %dims: 1=T2 2=trial# 3=ROI#
    for i = 1:size(success_roi,2); %for each roi...
        succ_temp = squeeze(success_roi(:,i,:))';  %dims 1=Time 2=trial number
        succ_temp = interp1(succ_temp, XVector);
        success_roi_interp(:,:,i) = succ_temp;  %dims 1=Time 2=trialNumber  3=roi
    end
    fail_roi = summary_fail{kk}.fail_roi(:,[days_roi_matcher.(strcat('i', days{kk}))], :); %collects only the ROIs of interest
    fail_roi_interp = nan(length(XVector), size(fail_roi,1), size(fail_roi,2));
    for i = 1:size(fail_roi,2);
        fail_temp = squeeze(fail_roi(:,i,:))';  %dims 1=Time 2=trial number
        fail_temp = interp1(fail_temp, XVector);
        fail_roi_interp(:,:,i) = fail_temp;   %dim 3=roi
    end
    cue_roi = summary_cue{kk}.cue_roi(:,[days_roi_matcher.(strcat('i', days{kk}))], :); %collects only the ROIs of interest
    cue_roi_interp = nan(length(XVector), size(cue_roi,1), size(cue_roi,2));
    for i = 1:size(cue_roi,2);
        cue_temp = squeeze(cue_roi(:,i,:))';  %dims 1=Time 2=trial number
        cue_temp = interp1(cue_temp, XVector);
        cue_roi_interp(:,:,i) = cue_temp;   %dim 3=roi
    end 
    
    %creates various averages of the interpolated data
    success_roi_interp_avg2 = squeeze(mean(success_roi_interp,2));   %only averages together the trials 
    fail_roi_interp_avg2    = squeeze(mean(fail_roi_interp,2)); %dims: 1=Time 2=roi#
    cue_roi_interp_avg2     = squeeze(mean(cue_roi_interp,2));
    success_roi_interp_avg  = squeeze(mean(success_roi_interp_avg2,2));  %dims: 1=Time
    fail_roi_interp_avg     = squeeze(mean(fail_roi_interp_avg2,2));
    cue_roi_interp_avg      = squeeze(mean(cue_roi_interp_avg2,2));
    success_roi_interp_avg3 = squeeze(mean(success_roi_interp,3));  %dims: 1=time 2=trial#
    fail_roi_interp_avg3    = squeeze(mean(fail_roi_interp,3));
    cue_roi_interp_avg3     = squeeze(mean(cue_roi_interp,3));
    
%     plot([-500:10:1000], success_roi_interp_avg,  'g'); hold on  %for showing all average curves on one plot
%     plot([-500:10:1000], fail_roi_interp_avg, 'r');
%     plot([-500:10:1000], cue_roi_interp_avg, 'b');
%     title(days(kk));

    %finding and storing peak times and values here for averaged trials using 100 ms window
    %Measuere of error is not a good one. 
    succ_time_to_peak = [succ_time_to_peak, find(success_roi_interp_avg==max(success_roi_interp_avg(51:110)))];  %finds peak within window of -1:5 frames relative to lever release frame (frame 6)
    fail_time_to_peak = [fail_time_to_peak, find(fail_roi_interp_avg==max(fail_roi_interp_avg(51:110)))];
    cue_time_to_peak  = [cue_time_to_peak, find(cue_roi_interp_avg==max(cue_roi_interp_avg(51:110)))];
    peak_value_succ = [peak_value_succ, mean(success_roi_interp_avg((succ_time_to_peak(kk)-time_before):(succ_time_to_peak(kk)+time_after)))];
    peak_value_fail = [peak_value_fail, mean(fail_roi_interp_avg((fail_time_to_peak(kk)-time_before):(fail_time_to_peak(kk)+time_after)))];
    peak_value_cue  = [peak_value_cue, mean(cue_roi_interp_avg((cue_time_to_peak(kk)-time_before):(cue_time_to_peak(kk)+time_after)))];
    plot_succ_sm = [plot_succ_sm, std(success_roi_interp_avg((succ_time_to_peak(kk)-time_before):(succ_time_to_peak(kk)+time_after)))/sqrt(time_before+time_after+1)]; %calculates the sem of the 11 points form the avg curve used to get the peak value
    plot_fail_sm = [plot_fail_sm, std(fail_roi_interp_avg((fail_time_to_peak(kk)-time_before):(fail_time_to_peak(kk)+time_after)))/sqrt(time_before+time_after+1)];
    plot_cue_sm  = [plot_cue_sm, std(cue_roi_interp_avg((cue_time_to_peak(kk)-time_before):(cue_time_to_peak(kk)+time_after)))/sqrt(time_before+time_after+1)];
    
    %find peak times for each trial   SHOULD ADAPT SO IT ALSO USES 100MS WINDOW
    %correct trials -------------------
    temp_lat_peak_succ = [];
    temp_tbyt_peak_val_succ = [];
    for i = 1:size(success_roi_interp_avg3,2);
        temp_val = find(success_roi_interp_avg3(:,i)==max(success_roi_interp_avg3([51:110],i)));  %finds the peak time of each trials TC
        if size(temp_val,1) > 1 %in case of finding multiple values it will select the first one
            temp_val = temp_val(find(temp_val >= 50,1, 'first'));
        end
        temp_lat_peak_succ = [temp_lat_peak_succ, temp_val]; %give peak time for each trial
        temp_tbyt_peak_val_succ = [temp_tbyt_peak_val_succ, success_roi_interp_avg3(temp_val,i)]; %collects all the trial by trial peak magnitudes
    end
    tbyt_peak_val_succ = [tbyt_peak_val_succ, mean(temp_tbyt_peak_val_succ)];  %collects the average peak value for all mice. Collected trial by trial. The peak value for that trial's peak
    tbyt_peak_val_succ_sm = [tbyt_peak_val_succ_sm, std(temp_tbyt_peak_val_succ)/sqrt(length(temp_tbyt_peak_val_succ))];
    tbyt_sm_lat_peak_succ = [tbyt_sm_lat_peak_succ, std((temp_lat_peak_succ-60)*10)/sqrt(length(temp_lat_peak_succ))]; %calculate SEM of peak times and avg the peak times
    tbyt_lat_peak_succ = [tbyt_lat_peak_succ, mean(temp_lat_peak_succ)];
    
    %early trials---------------------------
    temp_lat_peak_fail = [];
    temp_tbyt_peak_val_fail = [];
    for i = 1:size(fail_roi_interp_avg3,2);
        temp_val = find(fail_roi_interp_avg3(:,i)==max(fail_roi_interp_avg3([51:110],i)));
        if size(temp_val,1) > 1
            temp_val = temp_val(find(temp_val >= 50,1, 'first'));
        end
        temp_lat_peak_fail = [temp_lat_peak_fail, temp_val]; %give peak time for each trial
        temp_tbyt_peak_val_fail = [temp_tbyt_peak_val_fail, fail_roi_interp_avg3(temp_val,i)]; %collects all the trial by trial peaks 
    end
    tbyt_sm_lat_peak_fail = [tbyt_sm_lat_peak_fail, std((temp_lat_peak_fail-60)*10)/sqrt(length(temp_lat_peak_fail))]; %calculate SEM of peak times and avg the peak times
    tbyt_lat_peak_fail = [tbyt_lat_peak_fail, mean(temp_lat_peak_fail)];
    tbyt_peak_val_fail = [tbyt_peak_val_fail, mean(temp_tbyt_peak_val_fail)];  %collects the average peak value for all mice. Collected trial by trial. The peak value for that trial's peak
    tbyt_peak_val_fail_sm = [tbyt_peak_val_fail_sm, std(temp_tbyt_peak_val_fail)/sqrt(length(temp_tbyt_peak_val_fail))];
    
    %cue triggered----------------------
    temp_lat_peak_cue = [];
    temp_tbyt_peak_val_cue = [];
    for i = 1:size(cue_roi_interp_avg3,2);
        temp_val = find(cue_roi_interp_avg3(:,i)==max(cue_roi_interp_avg3([51:110],i)));
        if size(temp_val,1) > 1
            temp_val = temp_val(find(temp_val >= 50,1, 'first'));
        end
        temp_lat_peak_cue = [temp_lat_peak_cue, temp_val]; %give peak time for each trial
        temp_tbyt_peak_val_cue = [temp_tbyt_peak_val_cue, cue_roi_interp_avg3(temp_val,i)]; %collects all the trial by trial peaks 
    end
    tbyt_sm_lat_peak_cue = [tbyt_sm_lat_peak_cue, std((temp_lat_peak_cue-60)*10)/sqrt(length(temp_lat_peak_cue))]; %calculate SEM of peak times and avg the peak times
    tbyt_lat_peak_cue = [tbyt_lat_peak_cue, mean(temp_lat_peak_cue)];
    tbyt_peak_val_cue = [tbyt_peak_val_cue, mean(temp_tbyt_peak_val_cue)];  %collects the average peak value for all mice. Collected trial by trial. The peak value for that trial's peak
    tbyt_peak_val_cue_sm = [tbyt_peak_val_cue_sm, std(temp_tbyt_peak_val_cue)/sqrt(length(temp_tbyt_peak_val_cue))];
    
    
    %plot all trials vs avg for each animal on separate plots
    figure; 
    for i = 1:size(fail_roi_interp_avg3,2);
        subplot(1,3,2); plot([-500:10:1000], fail_roi_interp_avg3);
    end
    hold on; plot([-500:10:1000], mean(fail_roi_interp_avg3,2), 'LineWidth', 4); title([days(kk), 'fail']); vline(0);
    ylim([-0.2 0.5]);
    
    for i = 1:size(success_roi_interp_avg3,2);
        subplot(1,3,1); plot([-500:10:1000], success_roi_interp_avg3);
    end
    hold on; plot([-500:10:1000], mean(success_roi_interp_avg3,2), 'LineWidth', 4); title([days(kk), 'success']); vline(0);
    ylim([-0.2 0.5]);
    
    for i = 1:size(cue_roi_interp_avg3,2);
        subplot(1,3,3); plot([-500:10:1000], cue_roi_interp_avg3);
    end
    hold on; plot([-500:10:1000], mean(cue_roi_interp_avg3,2), 'LineWidth', 4); title([days(kk), 'cue']); vline(0);
    ylim([-0.2 0.5]);
end
% title('all seven animals that meet 60% correct');
% ylabel('df/f unshifted');
% xlabel('time from lever release (ms)');
avg_peak_succ = (mean(succ_time_to_peak)-60)*10;
avg_peak_fail = (mean(fail_time_to_peak)-60)*10;
avg_peak_cue = (mean(cue_time_to_peak)-60)*10;
succ_peak_std = std(succ_time_to_peak)*10;
fail_peak_std = std(fail_time_to_peak)*10;
cue_peak_std = std(cue_time_to_peak)*10;
%averages together all trials for a given animal then finds a peak time. Then avgs all the peak times across animals 
disp(['Average time to peak for correct trials = ' num2str(avg_peak_succ) ' +or- ' num2str(succ_peak_std)]);
disp(['Average time to peak for failed trials = ' num2str(avg_peak_fail) ' +or- ' num2str(fail_peak_std)]);
disp(['Average time to peak for cue triggered trials = ' num2str(avg_peak_cue) ' +or- ' num2str(cue_peak_std)]);
%find peak times from each trial after interpolation. Averages those peak times then averages all the avg peak times for all animals 
disp(['Average time to peak for correct trials = ' num2str(mean((tbyt_lat_peak_succ-60)*10)) ' +or- ' num2str(std((tbyt_lat_peak_succ-60)*10))]);
disp(['Average time to peak for failed trials = ' num2str(mean((tbyt_lat_peak_fail-60)*10)) ' +or- '  num2str(std((tbyt_lat_peak_fail-60)*10))]);
disp(['Average time to peak for cue triggered trials = ' num2str(mean((tbyt_lat_peak_cue-60)*10)) ' +or- '  num2str(std((tbyt_lat_peak_cue-60)*10))]);

%% METHOD 1
%plot scatter of correct vs early with shift and 100ms window around peak
succ_peak_val_mean = mean(peak_value_succ);
fail_peak_val_mean = mean(peak_value_fail);
min_scatter =  min([peak_value_succ,peak_value_fail]);
max_scatter =  max([peak_value_succ,peak_value_fail]);
colors = {'k', 'b', 'g', 'c', 'm', 'r', 'y'};
figure; hold on;
for i = 1:length(days);
plot(peak_value_succ(i),peak_value_fail(i), ['o' colors{i}], 'MarkerFaceColor', colors{i});  
end
xlim([-0.02 [max_scatter*1.1]]);
ylim([-0.02 [max_scatter*1.1]]);
x = -.1:.1:1;
y = x;
hold on; plot(x,y,'k')
hline(0,'--k')
vline(0,'--k')
xlabel('df/f success condition');
ylabel('df/f fail condition');
title({'corr vs early scatter Shift'; ['corrAvg=', num2str(succ_peak_val_mean), ' earlyAvg=', num2str(fail_peak_val_mean)]});
legend(days{1:length(days)});
for i = 1:length(days);
    errorbarxy(peak_value_succ(i), peak_value_fail(i), plot_succ_sm(i)', plot_fail_sm(i)');
    hold on;
end

%% METHOD 2
min_scatter =  min([tbyt_peak_val_succ, tbyt_peak_val_fail]);
max_scatter =  max([tbyt_peak_val_succ, tbyt_peak_val_fail]);
figure; hold on;
for i = 1:length(days);
plot(tbyt_peak_val_succ(i),tbyt_peak_val_fail(i), ['o' colors{i}], 'MarkerFaceColor', colors{i});  
end
xlim([-0.02 [max_scatter*1.1]]);
ylim([-0.02 [max_scatter*1.1]]);
x = -.1:.1:1;
y = x;
hold on; plot(x,y,'k')
hline(0,'--k')
vline(0,'--k')
xlabel('df/f success condition');
ylabel('df/f fail condition');
title('trial by trial corr vs early scatter Shift');
legend(days{1:length(days)});
for i = 1:length(days);
    errorbarxy(tbyt_peak_val_succ(i), tbyt_peak_val_fail(i), tbyt_peak_val_succ_sm(i)', tbyt_peak_val_fail_sm(i)');
    hold on;
end

%% LATENCY

%METHOD 1 scatter plot for latencies of average TCs. 
succ_time_to_peak = (succ_time_to_peak-60)*10;
fail_time_to_peak = (fail_time_to_peak-60)*10;
min_scatter =  min([succ_time_to_peak, fail_time_to_peak]);
max_scatter =  max([tbyt_lat_peak_succ, fail_time_to_peak]);
figure; hold on; 
for i = 1:length(days);
    plot(succ_time_to_peak(i),fail_time_to_peak(i), ['o' colors{i}], 'MarkerFaceColor', colors{i});  
end
xlim([-20 [max_scatter*1.1]]);
ylim([-20 [max_scatter*1.1]]);
x = -20:20:200;
y = x;
hold on; plot(x,y,'k');
hline(0,'--k');
vline(0,'--k');
xlabel('time to peak success condition (ms)');
ylabel('time to peak fail condition (ms)');
title('averaged TC corr vs early latency to peak scatter');
legend(days{1:length(days)});


%METHOD 2 Scatter plot by collecting latencies from trial by trial
tbyt_lat_peak_succ = (tbyt_lat_peak_succ-60)*10;
tbyt_lat_peak_fail = (tbyt_lat_peak_fail-60)*10;
min_scatter =  min([tbyt_lat_peak_succ, tbyt_lat_peak_fail]);
max_scatter =  max([tbyt_lat_peak_succ, tbyt_lat_peak_fail]);
figure; hold on;
for i = 1:length(days);
    plot(tbyt_lat_peak_succ(i),tbyt_lat_peak_fail(i), ['o' colors{i}], 'MarkerFaceColor', colors{i});  
end
xlim([-20 [max_scatter*1.1]]);
ylim([-20 [max_scatter*1.1]]);
x = -20:20:200;
y = x;
hold on; plot(x,y,'k');
hline(0,'--k');
vline(0,'--k');
xlabel('time to peak success condition (ms)');
ylabel('time to peak fail condition(ms)');
title('trial by trial corr vs early latency to peak scatter');
legend(days{1:length(days)});
for i = 1:length(days);
    errorbarxy(tbyt_lat_peak_succ(i), tbyt_lat_peak_fail(i), tbyt_sm_lat_peak_succ(i)', tbyt_sm_lat_peak_fail(i)');
    hold on;
end

%% CUE TRIGGERED MAGNITUDE
%METHOD1
succ_peak_val_mean = mean(peak_value_succ);
cue_peak_val_mean = mean(peak_value_cue);
min_scatter =  min([peak_value_succ,peak_value_cue]);
max_scatter =  max([peak_value_succ,peak_value_cue]);
colors = {'k', 'b', 'g', 'c', 'm', 'r', 'y'};
figure; hold on;
for i = 1:length(days);
    plot(peak_value_succ(i),peak_value_cue(i), ['o' colors{i}], 'MarkerFaceColor', colors{i});  
end
xlim([-0.02 [max_scatter*1.1]]);
ylim([-0.02 [max_scatter*1.1]]);
x = -.1:.1:1;
y = x;
hold on; plot(x,y,'k')
hline(0,'--k')
vline(0,'--k')
xlabel('df/f success condition');
ylabel('df/f cue condition');
title({'corr vs cue scatter Shift'; ['corrAvg=', num2str(succ_peak_val_mean), ' cueAvg=', num2str(cue_peak_val_mean)]});
legend(days{1:length(days)});
for i = 1:length(days);
    errorbarxy(peak_value_succ(i), peak_value_cue(i), plot_succ_sm(i)', plot_cue_sm(i)');
    hold on;
end

%METHOD2
min_scatter =  min([tbyt_peak_val_succ, tbyt_peak_val_cue]);
max_scatter =  max([tbyt_peak_val_succ, tbyt_peak_val_cue]);
figure; hold on;
for i = 1:length(days);
    plot(tbyt_peak_val_succ(i),tbyt_peak_val_cue(i), ['o' colors{i}], 'MarkerFaceColor', colors{i});  
end
xlim([-0.02 [max_scatter*1.1]]);
ylim([-0.02 [max_scatter*1.1]]);
x = -.1:.1:1;
y = x;
hold on; plot(x,y,'k')
hline(0,'--k')
vline(0,'--k')
xlabel('df/f success condition');
ylabel('df/f fail condition');
title('trial by trial corr vs cue scatter Shift');
legend(days{1:length(days)});
for i = 1:length(days);
    errorbarxy(tbyt_peak_val_succ(i), tbyt_peak_val_cue(i), tbyt_peak_val_cue_sm(i)', tbyt_peak_val_cue_sm(i)');
    hold on;
end

%% CUE TRIGGERED LATENCY
%METHOD 1 
cue_time_to_peak = (cue_time_to_peak-60)*10;
min_scatter =  min([succ_time_to_peak, cue_time_to_peak]);
max_scatter =  max([tbyt_lat_peak_succ, cue_time_to_peak]);
figure; hold on; 
for i = 1:length(days);
    plot(succ_time_to_peak(i),cue_time_to_peak(i), ['o' colors{i}], 'MarkerFaceColor', colors{i});  
end
xlim([-20 [max_scatter*1.1]]);
ylim([-20 [max_scatter*1.1]]);
x = -20:20:200;
y = x;
hold on; plot(x,y,'k');
hline(0,'--k');
vline(0,'--k');
xlabel('time to peak correct condition (ms)');
ylabel('time to peak cue condition (ms)');
title('averaged TC corr vs cue latency to peak scatter');
legend(days{1:length(days)});

%METHOD 2
tbyt_lat_peak_cue = (tbyt_lat_peak_cue-60)*10;
min_scatter =  min([tbyt_lat_peak_succ, tbyt_lat_peak_cue]);
max_scatter =  max([tbyt_lat_peak_succ, tbyt_lat_peak_cue]);
figure; hold on;
for i = 1:length(days);
    plot(tbyt_lat_peak_succ(i),tbyt_lat_peak_cue(i), ['o' colors{i}], 'MarkerFaceColor', colors{i});  
end
xlim([-20 [max_scatter*1.1]]);
ylim([-20 [max_scatter*1.1]]);
x = -20:20:200;
y = x;
hold on; plot(x,y,'k');
hline(0,'--k');
vline(0,'--k');
xlabel('time to peak success condition (ms)');
ylabel('time to peak cue condition (ms)');
title('trial by trial corr vs cue latency to peak scatter');
legend(days{1:length(days)});
for i = 1:length(days);
    errorbarxy(tbyt_lat_peak_succ(i), tbyt_lat_peak_cue(i), tbyt_sm_lat_peak_succ(i)', tbyt_sm_lat_peak_cue(i)');
    hold on;
end

