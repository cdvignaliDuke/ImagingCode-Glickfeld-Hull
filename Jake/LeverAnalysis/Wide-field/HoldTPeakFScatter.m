%colors = {'r', 'r', 'b', 'b', 'r', 'r', 'b', 'b', 'm', 'm', 'g', 'g', 'k', 'k', 'c', 'c', 'y', 'y', 'r', 'r', 'b', 'b'};
%days = {'150518_img24', '150519_img24', '150518_img25', '150517_img25', '150716_img27', '150718_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%ROIcell = {[1:2], [2:3], [1:2], [1:2], [1,2], [1,2,4], [1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:6], [2,3,5], [4], [1:2]};
colors = {'r'};
days = {'160314_img38'};
ROIcell = {[3:6]};
DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryFolder\';
BX_DIR = 'Z:\Analysis\LeverAnalysis\BxAndAnalysisOutputs\BxOutputs\';
%DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryNoShift\';
summary_succ = {}; 
summary_fail = {};
summary_bx = {};

for kk = 1:length(days)
    curr_file_succ = strcat(DATA_DIR, days{kk}, '_success');
    summary_succ{kk} = load(curr_file_succ);
    curr_file_fail = strcat(DATA_DIR, days{kk}, '_fail');
    temp2 = load(curr_file_fail);
    summary_fail{kk} = temp2;
    curr_file_bx = strcat(BX_DIR, days{kk}, '_bx_outputs');
    load(curr_file_bx); 
    
    %obtain peak windows, convert to single value of F for each trial. Plot
    %that F vs its corresponding hold time
    summary_succ_mat = [];
    summary_succ_mat = summary_succ{kk}.success_roi(:,ROIcell{kk},:); 
    summary_succ_mat = squeeze(mean(summary_succ_mat,2));
    succ_F = [];
    for ii = 1:size(summary_succ_mat,1)  %finding peak value for every single trial
        maxResp = find(summary_succ_mat(ii,:) == max(summary_succ_mat(ii,7:10)));
        peakWindow = [(maxResp-1):(maxResp+1)];
        succ_F = [succ_F, mean(summary_succ_mat(ii,peakWindow))];
    end
    
    summary_fail_mat = [];
    summary_fail_mat = summary_fail{kk}.fail_roi(:,ROIcell{kk},:); 
    summary_fail_mat = squeeze(mean(summary_fail_mat,2));
    fail_F = [];
    for ii = 1:size(summary_fail_mat,1)  %finding peak value for every single trial
        maxResp = find(summary_fail_mat(ii,:) == max(summary_fail_mat(ii,7:10)));
        peakWindow = [(maxResp-1):(maxResp+1)];
        fail_F = [fail_F, mean(summary_fail_mat(ii,peakWindow))];
    end
    figure; 
    plot(trial_outcome.succ_hold_dur, succ_F, ['o' colors{kk}]); 
    figure; 
    plot(trial_outcome.fail_hold_dur, fail_F, ['o' colors{kk}]); 
end

