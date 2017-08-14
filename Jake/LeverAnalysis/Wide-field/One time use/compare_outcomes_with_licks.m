% use the variables save in WF_lever_plotting_TCs to compare the df/f for
% correct trials vs early trials with licking
%     -success_roi
%     -fail_roi
%     -licking data

clear;
%set directory and days to be analyzed
F_TC_dir    = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
lick_TC_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'; 

days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
ROIcell = {[2], [2], [1:3], [1,3], [1:4], [1:5], [1:2], [1:2], [1:2], [1:2], [3:4], [2,3,5], [1], [1:2], [3:4], [1,2,3], [3:5]}; 

lever_release = 6; 
lick_window = [lever_release-1:lever_release+2];
fail_lick_mag = [];
fail_lick_sem = [];
list_of_days_fail_lick = []

for ii = [[1:2] [5:6] [8:length(days)]]
    %load the TCs
    load([F_TC_dir, days{ii}, '_success']);
    load([F_TC_dir, days{ii}, '_fail']);
    load([lick_TC_dir, days{ii}, '_bx_outputs'], 'licking_data');
    lick_trace_succ = licking_data.lick_trace_succ;
    lick_trace_fail = licking_data.lick_trace_fail;
    if length(size(success_roi)) ==3
        success_roi = squeeze(mean(success_roi(:,ROIcell{ii},:),2));
        fail_roi = squeeze(mean(fail_roi(:,ROIcell{ii},:),2));
    end
    
    %isolate licking traces with at least one lick in lick window
    fail_lick_inx =[];
    for iii = 1:size(fail_roi,1);
        if sum(lick_trace_fail(iii, lick_window)) >0;
            fail_lick_inx = [fail_lick_inx, iii];
        end
    end
    if isempty(fail_lick_inx)
        disp(['no fail trials with licks for ', days{ii}])
        continue
    end
    fail_lick = lick_trace_fail(fail_lick_inx,:);
    
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
    title(['Incorrect trials with at least 1 lick. n=', num2str(length(fail_lick_inx)), ' ', days{ii}]);
    xlabel('time (ms) relative to lever release'); ylabel('df/f and  avg number of licks/10 per frame');
    
    fail_lick_mag = [fail_lick_mag, max(fail_lick_tc_mean(5:8))];
    fail_lick_sem = [fail_lick_sem, fail_lick_tc_sem( find(fail_lick_tc_mean(5:8)==max(fail_lick_tc_mean(5:8)))+4 )];
    list_of_days_fail_lick = [list_of_days_fail_lick; days{ii}];
end

%need to store data from each mouse, find the peak, and make a scatterplot
save(['Z:\Analysis\WF Lever Analysis\licking_investigation\failed_trials_with_licks\fail_lick_tcs'], 'fail_lick_mag', 'fail_lick_sem', 'list_of_days_fail_lick')

%%code for plotting fail trials with licks vs correct trials
figure;
scatter( tbyt_peak_val_succ, fail_lick_mag); hold on; 
errorbar(tbyt_peak_val_succ, fail_lick_mag, fail_lick_sem, 'LineStyle', 'none');
for ii = 1:length(tbyt_peak_val_succ)
    plot([ tbyt_peak_val_succ(ii)-tbyt_peak_val_succ_sm(ii),  tbyt_peak_val_succ(ii)+tbyt_peak_val_succ_sm(ii)], [fail_lick_mag(ii), fail_lick_mag(ii)]);
end
%errorbarxy(tbyt_peak_val_succ, fail_lick_mag, tbyt_peak_val_succ_sm, fail_lick_sem, 'LineStyle', 'none');
x = [0:.1:1]; y =x;
plot(x,y);
xlabel('peak df/f for correct trials');
ylabel('peak df/f for early trials with licks');
ylim([0 0.4]); xlim([0 0.4]);












