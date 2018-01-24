% use the variables save in WF_lever_plotting_TCs to compare the df/f for
% correct trials vs early trials with licking
%     -success_roi
%     -fail_roi
%     -licking data
%MODIFIED 12/13/17 so it plots all LS ROIs instead of averaging them.

clear;
%set directory and days to be analyzed
F_TC_dir    = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
lick_TC_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'; 
WF_plotting_lists_of_days;

lever_release = 6; 
lick_window = [lever_release-2:lever_release+2];
peak_window = [lever_release-2:lever_release+2];
fail_lick_mag = [];
fail_lick_sem = [];
list_of_days_fail_lick = [];

for session_num = [[1:2] [5:6] [8:length(days)]]
    %load the TCs
    load([F_TC_dir, days{session_num}, '_success']);
    load([F_TC_dir, days{session_num}, '_fail']);
    load([lick_TC_dir, days{session_num}, '_bx_outputs'], 'licking_data');
    lick_trace_succ = licking_data.lick_trace_succ;
    lick_trace_fail = licking_data.lick_trace_fail;
    
    %averaging across ROIs
%     if length(size(success_roi)) ==3
%         success_roi = squeeze(mean(success_roi(:,ROIcell{ii},:),2));
%         fail_roi = squeeze(mean(fail_roi(:,ROIcell{ii},:),2));
%     end
    
    %only keep the failed trials with licks
    fail_lick_inx =[];
    for trial_num = 1:size(fail_roi,1);
        if sum(lick_trace_fail(trial_num, lick_window)) >0;
            fail_lick_inx = [fail_lick_inx, trial_num];
        end
    end
    if isempty(fail_lick_inx)
        disp(['no fail trials with licks for ', days{session_num}])
        continue
    end
    fail_lick = lick_trace_fail(fail_lick_inx,:);
    
    %only keep the correct trials with licks
    
    
    %get mean and sem TCs of licking and f 
    fail_lick_mean= mean(fail_lick,1);
    fail_lick_tc = fail_roi(fail_lick_inx,:);
    fail_lick_tc_mean = mean(fail_lick_tc,1);
    if size(fail_lick_tc,1) >1
        fail_lick_tc_sem= std(fail_lick_tc,1)/sqrt(length(fail_lick_tc));
    else
        fail_lick_tc_sem = zeros(1,size(fail_lick_tc,2));
    end
    
    shift = mean(fail_lick_tc_mean(1:3));
    fail_lick_tc_mean = fail_lick_tc_mean-shift;
    
    figure;
    x_axis = (([1:size(success_roi,2)])-lever_release)*100;
    bar(x_axis, fail_lick_mean/10); hold on;
    plot(x_axis, fail_lick_tc_mean);
    errorbar(x_axis, fail_lick_tc_mean, fail_lick_tc_sem);
    xlim([x_axis(1), x_axis(end)]); 
    title(['Incorrect trials with at least 1 lick. n=', num2str(length(fail_lick_inx)), ' ', days{session_num}]);
    xlabel('time (ms) relative to lever release'); ylabel('df/f and  avg number of licks/10 per frame');
    
    fail_lick_mag = [fail_lick_mag, max(fail_lick_tc_mean(5:8))];
    fail_lick_sem = [fail_lick_sem, fail_lick_tc_sem( find(fail_lick_tc_mean(5:8)==max(fail_lick_tc_mean(5:8)))+4 )];
    list_of_days_fail_lick = [list_of_days_fail_lick; days{session_num}];
end

%need to store data from each mouse, find the peak, and make a scatterplot
save(['Z:\Analysis\WF Lever Analysis\licking_investigation\failed_trials_with_licks\fail_lick_tcs'], 'fail_lick_mag', 'fail_lick_sem', 'list_of_days_fail_lick')

%%code for plotting fail trials with licks vs correct trials
figure;
scatter( tbyt_peak_val_succ, fail_lick_mag); hold on; 
errorbar(tbyt_peak_val_succ, fail_lick_mag, fail_lick_sem, 'LineStyle', 'none');
for session_num = 1:length(tbyt_peak_val_succ)
    plot([ tbyt_peak_val_succ(session_num)-tbyt_peak_val_succ_sm(session_num),  tbyt_peak_val_succ(session_num)+tbyt_peak_val_succ_sm(session_num)], [fail_lick_mag(session_num), fail_lick_mag(session_num)]);
end
%errorbarxy(tbyt_peak_val_succ, fail_lick_mag, tbyt_peak_val_succ_sm, fail_lick_sem, 'LineStyle', 'none');
x = [0:.1:1]; y =x;
plot(x,y);
xlabel('peak df/f for correct trials');
ylabel('peak df/f for early trials with licks');
ylim([0 0.4]); xlim([0 0.4]);












