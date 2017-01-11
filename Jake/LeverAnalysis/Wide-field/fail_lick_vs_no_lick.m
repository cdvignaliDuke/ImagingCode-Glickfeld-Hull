clear;
%trying to compare failed trials with licking to failed trials without
%licking. 

%these are the days which I am including in my summary statistics so far. Need to update with newer datasets. 
%all datasets which were included in the scatterplot as of 10/13/16
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', }; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 

%naive mice 

%sessions to be analyzed
%days = {'160904_img55', '160905_img55', '160916_img61', '160918_img61', '160920_img61', '160921_img61', '161030_img62', '160904_img55'};
%days = {'161031_img68','161101_img68', '161030_img69', '161030_img70', '161101_img69', '161101_img70'};

bx_source     = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
TC_dir = ['Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\'];
colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 1,0,1; 1,1,0; 0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar     
old_cd = cd; %save old cd so I can restore it later

%% SECTION TWO 
for ii = 1:length(days);
    b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    load(bx_out_dir);
    load([TC_dir, days{ii}, '_fail']);
    days(ii)
    lick_trace_fail = licking_data.lick_trace_fail; 
    total_licks_per_trial = sum(lick_trace_fail,2);
    find(total_licks_per_trial);
    no_lick_fail_trials = [];
    lick_fail_trials = [];
    for iii = 1:size(lick_trace_fail,1)
        if total_licks_per_trial(iii) ==0
            no_lick_fail_trials = [no_lick_fail_trials, iii];
        end
        if total_licks_per_trial(iii) >4
            lick_fail_trials = [lick_fail_trials, iii];
        end
    end
    figure; hold on;
    for iii = 1:size(fail_roi,2)
        plot([-5:10], squeeze(mean(fail_roi(no_lick_fail_trials, iii, :),1)), 'Color', 'r');
        plot([-5:10], squeeze(mean(fail_roi(lick_fail_trials, iii, :),1)), 'Color', 'g');
    end
    title([days(ii), 'failed trials with licks n=', num2str(length(lick_fail_trials)), ' (green). Without licks n=', num2str(length(no_lick_fail_trials)), '(red)']);
    xlabel('frame number relative to lever release');
    ylabel('df/f');
    savefig(['\\crash\data\public\presentations\SfN\licking_investigation\failed_trials_licks_vs_no _licks\', days{ii}]);
end




