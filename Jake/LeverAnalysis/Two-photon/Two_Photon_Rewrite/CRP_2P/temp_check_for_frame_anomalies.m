% search for and open frame_info for all datasets
%1) lever paper
%2) CRP data
%plot diff(counterTimes) and diff(counterValues)
%look for anomolies 

%% 1) lever data
clear

file_info;
for sub = 1:36                         % [29 ,30]%[8, 11, 13, 15:20] %36
    for rID = 1:2
        mouse = mouseID{sub}
        dateID = date;
        data_dir = fullfile('Z:\Analysis\2P Analysis\Lever');
        
        output_dest = fullfile('\\crash.dhe.duke.edu\data\home\ziye\2P_Analysis\2P_Analysis',[dateID{sub}, '_', runID{rID}, '_', mouse],'\');
        subfold = fullfile('\\crash.dhe.duke.edu\data\home\ziye\2P_Analysis\2P_Analysis',[dateID{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
        behav_dir = '\\crash.dhe.duke.edu\data\home\andrew\Behavior\Data\';
        tc_dir  = fullfile('\\crash.dhe.duke.edu\data\home\ziye\2P_Analysis\2P_Analysis',[dateID{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
       
        dest = tc_dir;
        dest_sub = dest;
        if exist(dest_sub)
            dataName   = dir([behav_dir, '*', behavID{sub}, '-', dateID{sub}(1:6), '*']);
            load([behav_dir, dataName(end).name]);
        end
    end
    
    %plot diff variables
    x_sz = length(cell2mat(cellfun(@int64,input.counterTimesUs,'UniformOutput',0)));
    figure('rend', 'painters', 'pos', [50 100 (796*1.5) (400*1.5)]); 
    subplot(2,1,1);
    %plot(diff(frame_info.times));
    plot(diff(cell2mat(cellfun(@int64,input.counterTimesUs,'UniformOutput',0))/1000));
    title('diff in counter times');
    ylim([0 80]); xlim([-5000 5000+x_sz]);
    
    subplot(2,1,2);
    %plot(diff(frame_info.counter_by_time));
    plot(diff(cell2mat(cellfun(@int64,input.counterValues,'UniformOutput',0))));
    title('diff in counter values');
    ylim([0 3]); xlim([-5000 5000+x_sz]);
    suptitle([dateID{sub}, ' ', mouse]);

end

%%  2) Check counter times and values for CRP data

% 
% 
% clear
% %file_info_CRP2;
% file_info_CRP_all;
% usFacs = 100;
% behav_dir = 'Z:\Data\2P_imaging\behavior\';
% crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration';
% doRerun = 0; useWarp = 0; doPCA = 0; useGPU = 0; checkImg = 0;
% for session_subset = 2
%     if session_subset ==1
%         stable_int = stable_int_1;
%         all_nPCA = all_nPCA_1;
%         all_PCuse = all_PCuse_1;
%         all_PCuse2 = all_PCuse2_1;
%         days_subset = days_1;
%         nPCA_to_use = nPCA_to_use_1;
%     elseif session_subset ==2
%         stable_int = stable_int_post;
%         all_nPCA = all_nPCA_post;
%         all_PCuse = all_PCuse_post;
%         all_PCuse2 = all_PCuse2_post;
%         days_subset = days_post;
%         nPCA_to_use = nPCA_to_use_post;
%     end
% end
% for sub = 1:15 %  %[2, 4, 5, 8, 9, 10, 12, 15]  
%      %collect session/mouse information    
%     session = days_subset{sub};
%     session_date = days_subset{sub}(1:6);
%     if session(end-2) == 'g'
%         mouse_num = ['9', session(end-1:end)];
%         mouse_ID = ['img', session(end-1:end)]; 
%     elseif session(end-2) =='0'
%         mouse_num = session(end-2:end);
%         mouse_ID = ['img', session(end-2:end)]; 
%     end
%     
%     %set pathnames
%     data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
%     for rID = 1:2
%         if exist([data_dir, mouse_ID, '_000_', runID{rID}, '.sbx'], 'file') == 2
%             break
%         end
%     end
%     out_dir =  fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID], '\');
%     data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
%     if ~exist(out_dir)
%         mkdir(out_dir);
%     end
%     
%     %set directories
%     subMat = dir([behav_dir, '*', mouse_num, '-', session_date, '*']); %sets directory as the behavior file
%     
%     %load the behavior data
%     load([behav_dir, subMat.name]);
%     
% %     %load frame info data
% %     load([out_dir, 'parse_behavior.mat'], 'frame_info');
% 
%     %plot diff variables
%     figure('rend', 'painters', 'pos', [50 100 (796*1.5) (400*1.5)]); 
%     subplot(2,1,1);
%     %plot(diff(frame_info.times));
%     plot(diff(cell2mat(input.counterTimesUs)/1000));
%     title('diff in counter times');
%     ylim([0 80]); xlim([-5000 5000+length(cell2mat(input.counterValues))]);
%     subplot(2,1,2);
%     %plot(diff(frame_info.counter_by_time));
%     plot(diff(cell2mat(input.counterValues)));
%     title('diff in counter values');
%     ylim([0 3]); xlim([-5000 5000+length(cell2mat(input.counterValues))]);
%     suptitle([session_date, ' ', mouse_num]);
%     
%     
% end






