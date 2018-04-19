clear
file_info_WF;

for sub = 1:length(date)
    behav_dir = 'Z:\home\andrew\Behavior\Data\';
    file_info_WF;
    
    data_dir = fullfile('Z:\home\jake\Analysis\WF Lever Analysis\PCA_ICA_output\');
    subfold = fullfile(data_dir,[date{sub}, '\']);
    WF_dir = fullfile('Z:','home','jake','Analysis','WF Lever Analysis','Meta-kmeans_output_dir', date{sub}, '\');
    
    load([subfold, 'ROI_TCs.mat']);
    
%     if exist([subfold, 'ROI_xy.mat'],'file') == 2
%             load([subfold, 'ROI_xy.mat']);
%         elseif exist([WF_dir, date{sub}, 'reg_out.mat'], 'file') == 2
%             load([WF_dir, date{sub}, 'reg_out.mat']);
%     end
%     mask_ori = zeros(1002,1004);
%     
%     mask_ori(ROI_x, ROI_y) = reshape(mask_final, sz(1), sz(2));
%     mask_final = reshape(mask_ori, 1002, 1004);
%     sz = [1002, 1004];
%     save([subfold, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
    
    date = date{sub};
    mouse = date(end-4:end);
    dataName   = dir([behav_dir, '*', '9', date(end-1:end), '-', date(1:6), '*']);
    load([behav_dir, dataName(end).name]);
    tc_dir  = subfold;
    dest = tc_dir;
    dest_sub = dest;
    % dest_sub = Z:\home\jake\Analysis\2P Analysis\Ziye_2P_figure
    
    if exist(dest_sub)
        load([tc_dir, 'ROI_TCs.mat']);
        %             load([tc_dir, '_npSubTCs.mat']);
        
        
    else
        error('subject %s_%s and subNum is %d input not found\n', date, mouse, sub);
    end
    %
%     getTC_events;
%     TC_quantification;
%      close all
%     
     getTC_Spike
%     spike_quantification
      close all
      clearvars -except sub 
    
    
end