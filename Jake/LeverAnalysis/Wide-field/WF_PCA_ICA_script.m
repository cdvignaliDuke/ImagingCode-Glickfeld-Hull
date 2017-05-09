clear
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
bx_source          = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir     = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
PCA_output_dir_base= ['Z:\Analysis\WF Lever Analysis\PCA_output_dir\'];
old_cd = cd; %save old cd so I can restore it later

%% Select stable interval to use in motion correction later
for ii = [16, 17]  %160722_img53 41:120     %160606_img46 1:80    %160904_img55 1:80
    if ~exist([PCA_output_dir_base, days{ii}, '\', days{ii}, '_stable_movie_avg.mat'])
        stable_frame_int = input(['please open first movie file for ' days{ii} ' in imagej and select a frame interval in the form of a vector with no/few motion artifacts']);
        meta_data_dir = [image_dest_base days{ii} '\' days{ii} '_meta_data'];
        load(meta_data_dir, 'meta_data2');
        stable_movie = [];
        for iii = stable_frame_int
            img_subset = imread(meta_data2{1}(1).Filename,iii);
            stable_movie = cat(3, stable_movie, img_subset);
        end
        stable_movie_avg = mean(stable_movie, 3);
        save([PCA_output_dir_base, days{ii}, '\', days{ii}, '_stable_movie_avg'], 'stable_movie_avg');
    end
end
%% write new tiffs 
    %Method for selecting frames for PCA based on lever activity.
    %release_frames collects frames from -5:2 frames relative to lever
    %release in order to capture the initial, transient increase in F.
    %lick_frames collects frames from 3:23 frames relative to lever release
    %in order to collect the activity during licking. press_frames collects
    %frames -5:4 frames relative to lever press
for ii = [16, 17] 
    %pathnames and loading data
    days{ii}
    image_source = [image_source_base, days{ii} '\' days{ii} '_MMStack.ome.tif'];
    meta_data_dir = [image_dest_base days{ii} '\' days{ii} '_meta_data'];
    PCA_output_dir = [PCA_output_dir_base days{ii} '\'];
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    bx_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    load(meta_data_dir, 'meta_data2');
    load(bx_out_dir, 'frame_info', 'trial_outcome'); %lever, frame_info, data_tc,licking_data,lickTimes,tc_dfoverf

    %Check for existing black out mask of full frame. If none then manually generate one in imageJ
    if exist([PCA_output_dir, days{ii}, '_movies_full_frame_mask.tif'], 'file') == 2;
        full_frame_mask = readtiff([PCA_output_dir, days{ii}, '_movies_full_frame_mask.tif']);
    else
        disp(['PCA analysis paused. Please check to make sure a mask exist for ', days{ii}, ' if not then manually generate one in imageJ now']);
        pause
    end
    
    %only include indeces for trials which occured after the first frame
    frame_counter = frame_info.counter;
    first_frame_trial_num = frame_info.f_frame_trial_num;
    last_frame_trial_num = frame_info.l_frame_trial_num;
    num_frames_with_bx = length(frame_info.times);
    
    early_inx = trial_outcome.early_inx;
    corr_inx  = trial_outcome.corr_inx;
    early_inx([1:frame_info.f_frame_trial_num, frame_info.l_frame_trial_num:end])=0; %remove unimaged trials
    corr_inx([1:frame_info.f_frame_trial_num, frame_info.l_frame_trial_num:end])=0;
    early_inx = early_inx(find(early_inx >0));
    corr_inx  = corr_inx(find(corr_inx >0));
    corr_hold_dur = trial_outcome.succ_hold_dur;
    early_hold_dur = trial_outcome.fail_hold_dur;
    assert(length(corr_inx) == length(corr_hold_dur));
    assert(length(early_inx) == length(early_hold_dur));
    %find frame numbers for the last frame of the iti
    iti_dur_ms = bx_data.itiTimeMs; %round(bx_data.itiTimeMs/mode(diff(frame_info.times))); %determine the duration of the iti in frames
    if iti_dur_ms < 2500
        error('iti time does not exceed the recommended minimum duration in order to isolate an iti interval without licking');
    end
    trial_start_times_ms = bx_data.tThisTrialStartTimeMs(first_frame_trial_num+1:last_frame_trial_num-1); %exclude trials which were not completely imaged 
    iti_end_inx = cell2mat(trial_start_times_ms) + double(iti_dur_ms);
    iti_end_inx = round(iti_end_inx - frame_info.f_frame_MWorks_time); %zero indexed times relative to first frame time
    
    %convert indeces to frame numbers. convert hold durations to times relative to releases. 
    early_press = early_inx-early_hold_dur;
    corr_press  = corr_inx-corr_hold_dur;
    early_inx    = frame_counter(early_inx);
    corr_inx     = frame_counter(corr_inx);
    early_press  = frame_counter(early_press);
    corr_press   = frame_counter(corr_press);
    
    early_inx = early_inx(find(early_inx > 10 & early_inx < num_frames_with_bx-20)); % must have 10 frames before and after the lever event 
    corr_inx  = corr_inx(find(corr_inx > 10 & corr_inx < num_frames_with_bx-20)); % must have 10 frames before and after the lever event 
    combined_inx = sort([early_inx, corr_inx]);
    combined_inx = frame_counter(combined_inx);
    iti_end_inx = frame_info.counter(iti_end_inx); %convert iti end times to frame numbers
    
    %download the entire movie so specific frames can be indexed
    disp('starting tiff file download');
    full_movie = download_full_tiff_movie(meta_data2);
    
    %motion correction of the entire movie
    stable_movie_avg = load([PCA_output_dir_base, days{ii}, '\', days{ii}, '_stable_movie_avg']);
    stable_movie_avg = stable_movie_avg.stable_movie_avg;
    [out full_movie_reg] = stackRegister(full_movie, stable_movie_avg);  %does stackRegister alter the size of each frame? 
    clear full_movie; %the motion correction does not alter the size of the movie
    
    %apply the mask
    full_movie_size = size(full_movie_reg);
    full_frame_mask = reshape(full_frame_mask, 1, full_movie_size(1)*full_movie_size(2) );
    full_movie_reg = reshape(full_movie_reg, full_movie_size(1)*full_movie_size(2), full_movie_size(3) );
    zeros_in_mask = find(full_frame_mask == 0);
    full_movie_reg(zeros_in_mask, :) = 0;
    full_movie_reg = reshape(full_movie_reg, full_movie_size(1), full_movie_size(2), full_movie_size(3));
    full_frame_mask = reshape(full_frame_mask, full_movie_size(1), full_movie_size(2));
        
    %find coordinates to crop image based on the min and max x/y values of non-zero values
    x_max = max(find(sum(full_frame_mask,1)));
    x_min = min(find(sum(full_frame_mask,1)));
    y_max = max(find(sum(full_frame_mask,2)));
    y_min = min(find(sum(full_frame_mask,2)));
    full_movie_reg = full_movie_reg([y_min:y_max], [x_min:x_max], :);
    save([PCA_output_dir_base, days{ii}, '\', days{ii}, '_cropped_image_coords'], 'x_min', 'x_max', 'y_min', 'y_max');
    
    %extract specific frame intervals from the full, motion-corrected movie
    corr_release_frames = frame_int_extractor(full_movie_reg, corr_inx, -5, 2);
    early_release_frames = frame_int_extractor(full_movie_reg, early_inx, -5, 2);
    combined_release_frames = cat(3, corr_release_frames, early_release_frames);
    lick_frames = frame_int_extractor(full_movie_reg, corr_inx, 3, 23);
    iti_frames = frame_int_extractor(full_movie_reg, iti_end_inx, -9, -1);
    
    %isolate and extract frames for full trial
    early_full_inx = [];
    for iii = 1:length(early_inx)
        early_full_inx = [early_full_inx, [(early_press(iii)-5):(early_inx(iii)+20)]];
    end
    corr_full_inx = [];
    for iii = 1:length(corr_inx)
        corr_full_inx = [corr_full_inx, [corr_press(iii)-5:corr_inx(iii)+20]];
    end
    early_full_frames = full_movie_reg(:,:,early_full_inx);
    corr_full_frames  = full_movie_reg(:,:,corr_full_inx);
    
    %generate average images for each category
    corr_release_frames_avg = mean(corr_release_frames,3);
    early_release_frames_avg = mean(early_release_frames,3);
    combined_release_frames_avg = mean(combined_release_frames,3);
    lick_frames_avg = mean(lick_frames,3);
    iti_frames_avg = mean(iti_frames,3);
    corr_full_frames_avg = mean(corr_full_frames,3);
    early_full_frames_avg = mean(early_full_frames,3);
    
%     %write tiffs of all frames for each category
%     writetiff(corr_release_frames, [PCA_output_dir, days{ii}, '_corr_release_movie']);
%     writetiff(early_release_frames, [PCA_output_dir, days{ii}, '_early_release_movie']);
%     writetiff(combined_release_frames, [PCA_output_dir, days{ii}, '_combined_release_movie']);
%     writetiff(lick_frames, [PCA_output_dir, days{ii}, '_licking_movie']);
%     writetiff(iti_frames, [PCA_output_dir, days{ii}, '_iti_movie']);
    
    %save matfiles of all frames for each category
    disp('saving movie files')
    save([PCA_output_dir, days{ii}, '_corr_release_movie'], 'corr_release_frames');
    save([PCA_output_dir, days{ii}, '_early_release_movie'], 'early_release_frames');
    save([PCA_output_dir, days{ii}, '_combined_release_movie'], 'combined_release_frames');
%     save([PCA_output_dir, days{ii}, '_licking_movie'], 'lick_frames');
%     save([PCA_output_dir, days{ii}, '_iti_movie'], 'iti_frames');
    save([PCA_output_dir, days{ii}, '_corr_full_movie'], 'corr_full_frames');
    save([PCA_output_dir, days{ii}, '_early_full_movie'], 'early_full_frames');
    
%     %write tiffs for the average images of each category
%     writetiff(corr_release_frames_avg, [PCA_output_dir, days{ii}, '_corr_release_movie_avg']);
%     writetiff(early_release_frames_avg, [PCA_output_dir, days{ii}, '_early_release_movie_avg']);
%     writetiff(combined_release_frames_avg, [PCA_output_dir, days{ii}, '_combined_release_movie_avg']);
%     writetiff(lick_frames_avg, [PCA_output_dir, days{ii}, '_licking_movie_avg']);
%     writetiff(iti_frames_avg, [PCA_output_dir, days{ii}, '_iti_movie_avg']);
    
    %save matfiles for the average images of each category
    save([PCA_output_dir, days{ii}, '_corr_release_movie_avg'], 'corr_release_frames_avg');
    save([PCA_output_dir, days{ii}, '_early_release_movie_avg'], 'early_release_frames_avg');
    save([PCA_output_dir, days{ii}, '_combined_release_movie_avg'], 'combined_release_frames_avg');
%     save([PCA_output_dir, days{ii}, '_licking_movie_avg'], 'lick_frames_avg');
%     save([PCA_output_dir, days{ii}, '_iti_movie_avg'], 'iti_frames_avg');
    save([PCA_output_dir, days{ii}, '_corr_full_movie_avg'], 'corr_full_frames_avg');
    save([PCA_output_dir, days{ii}, '_early_full_movie_avg'], 'early_full_frames_avg');
end

%% PCA

% for ii = [15,17]
%     %pauses inserted to allow for the viewing and selection of PCs
%     days{ii}
%     PCA_output_dir = [PCA_output_dir_base days{ii} '\'];
% 
%     global rt
%     early_release_frames_avg = readtiff([PCA_output_dir, days{ii}, '_early_release_movie_avg.tif']);
%     rt = readtiff([PCA_output_dir, days{ii}, '_early_release_movie.tif']);
%     [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA2([PCA_output_dir, days{ii}, '_early_release_movie.tif'], [], 100,1, [PCA_output_dir], []);
%     [PCuse] = CellsortChoosePCs(early_release_frames_avg, mixedfilters);
% 
%     pause  %I cant remember why I put these pauses here
% 
%     corr_release_frames_avg = readtiff([PCA_output_dir, days{ii}, '_corr_release_movie_avg.tif']);
%     rt = readtiff([PCA_output_dir, days{ii}, '_corr_release_movie.tif']);
%     [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA2([PCA_output_dir, days{ii}, '_corr_release_movie.tif'], [], 100,1, [PCA_output_dir], []);
%     [PCuse] = CellsortChoosePCs(corr_release_frames_avg, mixedfilters);
% 
%     pause
% 
%     combined_release_frames_avg = readtiff([PCA_output_dir, days{ii}, '_combined_release_movie_avg.tif']);
%     rt = readtiff([PCA_output_dir, days{ii}, '_combined_release_movie.tif']);
%     [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA2([PCA_output_dir, days{ii}, '_combined_release_movie.tif'], [], 100,1, [PCA_output_dir], []);
%     [PCuse] = CellsortChoosePCs(combined_release_frames_avg, mixedfilters);
% 
%     pause
% 
%     lick_frames_avg = readtiff([PCA_output_dir, days{ii}, '_licking_movie_avg.tif']);
%     rt = readtiff([PCA_output_dir, days{ii}, '_licking_movie.tif']);
%     [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA2([PCA_output_dir, days{ii}, '_licking_movie.tif'], [], 100,1, [PCA_output_dir], []);
%     [PCuse] = CellsortChoosePCs(lick_frames_avg, mixedfilters);
% 
%     pause
%     
%     iti_frames_avg = readtiff([PCA_output_dir, days{ii}, '_iti_movie_avg.tif']);
%     rt = readtiff([PCA_output_dir, days{ii}, '_iti_movie.tif']);
%     [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA2([PCA_output_dir, days{ii}, '_iti_movie.tif'], [], 100,1, [PCA_output_dir], []);
%     [PCuse] = CellsortChoosePCs(iti_frames_avg, mixedfilters);
% end
% 
% % function [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm, PCuse] = run_PCA_all_categories()
% % for iii = 1:length()
% %     rt = readtiff([PCA_output_dir, days{ii}, '\', days{ii}, '_' movie_tag '.tif']);
% % end
% % end

%% ICA
% 
% for ii = [15,17]
%     days{ii}
%     image_source = [image_source_base, days{ii} '\' days{ii} '_MMStack.ome.tif'];
%     load([PCA_output_dir_base days{ii}, '\' days{ii} '_corr_release_movie_1,1312_20-Jan-2017']);
%     %PCuse = []; 
%     PCuse = [2:16];
%     mu = 0.05;
%     nIC = 12;
%     ica_A_guess = []; 
%     termtol = 0.00001;
%     maxrounds = 1000; 
%     
%     [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
%     mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds);
% 
%     figure; 
%     ICuse = 1:nIC;
%     dt = 100;
%     CellsortICAplot('contour', ica_filters, ica_sig, movm, [], dt, 1, 1, ICuse);
% end
% % scale
% %  better masks
% % Compare lever release areas to licking areas and see if they are seperable 


