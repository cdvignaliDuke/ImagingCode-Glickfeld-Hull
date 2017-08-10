% Make frame by frame figures of the average df/f or F within a rectangular
% ROI averaged across trials for correct and early lever releases. 
%Adapted from HAD_cmp_success_fail

%generates frame by frame movies of an ROI. Heatmaps. 
clear
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)
bx_source     = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
PCA_output_dir_base = ['Z:\Analysis\WF Lever Analysis\PCA_output_dir\'];
kmeans_output_dir_base = ['Z:\Analysis\WF Lever Analysis\kmeans_output_dir\'];
frame_by_frame_outputs = ['Z:\Analysis\WF Lever Analysis\Frame_by_frame_output\'];
days = {'151211_img32'};

for kk=1:length(days)
    % ---- load behavior file
    b_data = get_bx_data(bx_source, days{kk});  %find the correct behavior file and loads it.
    
    % --- load frame times
    frame_info_dest = [image_dest_base '\' days{kk} '\' days{kk} '_frame_times.mat'];
    load(frame_info_dest);
    ifi = (frame_times(end)-frame_times(1))/length(frame_times);
    Sampeling_rate = round(1000/ifi);
    
    % --- load bxOutputs and metadata
    load([bx_outputs_dir, days{kk}, '_bx_outputs']);
    
    %Determine the ROI and crop the movie accordingly
    load([image_dest_base, days{kk}, '\', days{kk}, '_meta_data']);
    image_source = [image_source_base, days{kk}];
    [ROI_x ROI_y]  = get_movie_ROI(meta_data2, frame_info);
    [img, sz]  = crop_movie_by_ROI(image_source, meta_data2, ROI_x, ROI_y, frame_info);
    
    %save cropped movie and ROI related variables
    if exist([frame_by_frame_outputs, days{kk}], 'file') == 0;
        mkdir([frame_by_frame_outputs, days{kk}])
    end
    save([frame_by_frame_outputs, days{kk}, '\', days{kk}, '_ROI'], 'img');
    save([frame_by_frame_outputs, days{kk}, '\', days{kk}, 'ROI_related_variables'], 'sz', 'ROI_x', 'ROI_y');
    
    %reshape img for the df/f calculation
    img = reshape(img, [size(img,1)*size(img,2), size(img,3)]);
    img_dfoverf = tbyt_dfoverf_full_frame(b_data, [bx_outputs_dir, days{kk}, '_bx_outputs'], img);
    
    % PLOTTING --------------------------------
    %do simple movie analysis
    func = @median;  % func = @mean;   %func = @std;
    pre_frames = 5;
    post_frames = 10;  
    rm_baseline_plot = 0; % 1 for removing baseline when plotiing
    ts = (-pre_frames:post_frames)*1000/round(Sampeling_rate);
    
    %Correct trials: identify and plot mean correct trial movie
    use_ev_success = trial_outcome.success_time;
    success_movie = trigger_movie_by_event(img_dfoverf, frame_info, use_ev_success, pre_frames, post_frames);
    avg_success_move = squeeze(func(success_movie,1));
    fig1 = figure; plot_movie(avg_success_move,sz,rm_baseline_plot, ts);
    title('Correct');
    clim([-0.15 0.151]);
    disp(['Correct: n = ' num2str(size(success_movie,1))]);
   
    %Failed trials: identify and plot mean correct trial movie
    use_ev_fail = trial_outcome.early_time;
    fail_movie = trigger_movie_by_event(img_dfoverf, frame_info, use_ev_fail, pre_frames, post_frames);
    avg_fail_move = squeeze(func(fail_movie,1));
    fig2 = figure; plot_movie(avg_fail_move,sz,rm_baseline_plot, ts);
    title('Early');
    disp(['Early: n = ' num2str(size(fail_movie,1))]);
    clim([-0.15 0.151]);
    
    % ----- success  - early!!
    diff_succ_fail = avg_success_move - avg_fail_move;
    fig3 = figure; plot_movie(diff_succ_fail,sz,rm_baseline_plot, ts);    %JH ALTERED   Added "fig3"
    title('Correct - Early');
    figure(fig3);                             
    clim([-0.15 0.151]);
%  
%     % -- Trigger movie off lever press IF subsequent hold time is >1000ms
%     if lever.release(1) < lever.press(1)
%         holdDuration = NaN(1, size((lever.release),2)-1);
%         for i = 1:(length(holdDuration))
%             holdDuration(i) = lever.release(i+1)-lever.press(i);
%         end
%     else
%         holdDuration = NaN(1, size((lever.release),2));
%         for i = 1:(length(holdDuration))
%             holdDuration(i) = lever.release(i)-lever.press(i);
%         end
%     end
%     
%     longHoldInd = zeros(size(holdDuration));
%     for i = 1:size(holdDuration, 2);
%         if holdDuration(i) > 1000
%             longHoldInd(i) = 1;
%         end
%     end
%     
%     use_ev_press = NaN(1, sum(longHoldInd));
%     aa = 1;
%     for i = 1:length(longHoldInd)
%         if longHoldInd(i) == 1;
%             use_ev_press(aa) = lever.press(i);
%             aa=aa+1;
%         end
%     end
%     use_ev_press = round(use_ev_press);
%     %----- uncomment to use only events w/o lever press after release
%     %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
%     %         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
%     %------ uncomment to use event only w/o lever press before release time
%     %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
%     %         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
%     %
%     press_movie = trigger_movie_by_event(img_dfoverf, frame_info, ...
%         use_ev_press, pre_frames, post_frames);
%     avg_press_move = squeeze(func(press_movie,1));
%     fig4 = figure; plot_movie(avg_press_move,sz,rm_baseline_plot, ts);
%     title('lever press');
%     clim([-0.5 0.5]);
%     disp(['press: n = ' num2str(size(press_movie,1))]);
end
