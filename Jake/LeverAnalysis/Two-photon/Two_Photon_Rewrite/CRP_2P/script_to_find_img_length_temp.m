
clear
%file_info_CRP2;
file_info_CRP_all;
usFacs = 100;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration';
doRerun = 0; useWarp = 0; doPCA = 0; useGPU = 0; checkImg = 0;
for session_subset = 2
    if session_subset ==1
        stable_int = stable_int_1;
        all_nPCA = all_nPCA_1;
        all_PCuse = all_PCuse_1;
        all_PCuse2 = all_PCuse2_1;
        days_subset = days_1;
        nPCA_to_use = nPCA_to_use_1;
    elseif session_subset ==2
        stable_int = stable_int_post;
        all_nPCA = all_nPCA_post;
        all_PCuse = all_PCuse_post;
        all_PCuse2 = all_PCuse2_post;
        days_subset = days_post;
        nPCA_to_use = nPCA_to_use_post;
    end

for sub = 1:15 % [13,14]; %[2, 4, 5, 8, 9, 10, 12, 15]     %[1:2, 4:(length(stable_int)-2)]; %[2,5, 9, 22, 23, 31, 3, 6, 10, 25, 27, 32]
    %collect session/mouse information    
    session = days_subset{sub};
    session_date = days_subset{sub}(1:6);
    if session(end-2) == 'g'
        mouse_num = ['9', session(end-1:end)];
        mouse_ID = ['img', session(end-1:end)]; 
    elseif session(end-2) =='0'
        mouse_num = session(end-2:end);
        mouse_ID = ['img', session(end-2:end)]; 
    end
    
    %set pathnames
    data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
    for rID = 1:2
        if exist([data_dir, mouse_ID, '_000_', runID{rID}, '.sbx'], 'file') == 2
            break
        end
    end
    out_dir =  fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID], '\');
    data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
    if ~exist(out_dir)
        mkdir(out_dir);
    end
    
    config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
    laser_power_fn = dir(fullfile(data_dir,['*' runID{rID} '_realtime.mat']));
    [img_nframes, ~, ~] = get_nFrames_temp( data_dir, runID{rID}, [], []);
    
     save([out_dir, 'img_reg.mat'], 'img_nframes', '-v7.3', '-append');
end
end