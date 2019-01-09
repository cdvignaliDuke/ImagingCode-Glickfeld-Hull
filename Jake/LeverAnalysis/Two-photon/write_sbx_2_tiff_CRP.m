% make tiff for suite 2P from sbx file
clear
file_info_CRP_all;
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration';
for session_subset = 1%[1:2]
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
    elseif session_subset ==3
        stable_int = stable_int_UR;
        all_nPCA = all_nPCA_UR;
        all_PCuse = all_PCuse_UR;
        all_PCuse2 = all_PCuse2_UR;
        days_subset = days_UR;
        nPCA_to_use = nPCA_to_use_UR;
    elseif session_subset ==4
        stable_int = stable_int_1000;
        all_nPCA = all_nPCA_1000;
        all_PCuse = all_PCuse_1000;
        all_PCuse2 = all_PCuse2_1000;
        days_subset = days_1000;
        nPCA_to_use = nPCA_to_use_1000;
    end
    
    for sub = [8:9,14]%[1:15] %  %[2, 4, 5, 8, 9, 10, 12, 15]     %[1:2, 4:(length(stable_int)-2)]; %[2,5, 9, 22, 23, 31, 3, 6, 10, 25, 27, 32]
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
        session
        
        %set pathnames
        data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
        for rID = 1:2
            if exist([data_dir, mouse_ID, '_000_', runID{rID}, '.sbx'], 'file') == 2
                break
            end
        end
        if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
            rID=2;
        end
        data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
        
        %% Load imaging data and trim movie if needed
        config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
        laser_power_fn = dir(fullfile(data_dir,['*' runID{rID} '_realtime.mat']));
        [img, ~, ~] = loadFile(data_dir, runID{rID}, [], []);
        [img_mat_file, laser_power_vec_ttl] = get_laser_power_data(data_dir, config_fn, laser_power_fn);
        if ~isempty(laser_power_vec_ttl)
            laser_on_ind_conserv = conservative_laser_on(laser_power_vec_ttl);
            laser_on_ind = find(laser_power_vec_ttl);
            img =  img(:,:,laser_on_ind_conserv(find(laser_on_ind_conserv<=size(img,3))));
        end
        %min_x = 151; max_x = 600;
        %img1 = img(:,[min_x:max_x],[1:5000]);
        
        %% write tiff
        %set directory
        suite_tif_dir = ['Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\suite_2P_inputs\', mouse_num, '\20', session_date(1:2),'-', session_date(3:4),'-', session_date(5:6), '\'];
        if ~exist(suite_tif_dir)
            mkdir(suite_tif_dir);
        end
        
        %write tiffs in 5000 frame files
        n_frames = 5000;
        start_frame = 1;
        for file_num = 1:5
            img1 = img(:,:,[start_frame:start_frame+n_frames-1]);
            writetiff(img1, [suite_tif_dir, '\', mouse_ID, '_', num2str(file_num)]);
            start_frame = start_frame+n_frames;
        end

    end
end
